library(tidyverse)
library(wql)
library(patchwork)
setwd("Documents/salp_ecomorphology/")

#Load data and label control rows
presens <- read.csv("respirometry_Kona21.tsv", header = T, sep = "\t", stringsAsFactors = F)
presens$Species[which(is.na(presens$Species))] <- "Control"
presens$Specimen[which(is.na(presens$Specimen))] <- "Control"

#Filter by plastic and intact treatments
presens <- presens[which(presens$Container=="Plastic" & presens$Treatment == "Intact"),]

#Correct for O2 saturation limit wit temperature
presens$Temperature...C.[which(is.na(presens$Temperature...C.))]<-29 #### MAKING SHIT UP ####
presens <- mutate(presens, sat_O2 = O2..mg.L./oxySol(presens$Temperature...C., 30.1, 1))

#Subset controls and specimen measures
controls <- presens[which(presens$Specimen=="Control"),]
specimens <- presens[which(presens$Specimen != "Control"),]

#Join Specimens+COntrols and subtract for specific O2 measurements
join_presens <- full_join(specimens, controls[-c(2:10,14,16:18)], by=c("Experiment","Time.point..min."), suffix=c("_animal","_control"))
names(join_presens)
join_presens <- mutate(join_presens, O2..mg.L.specific = O2..mg.L._animal-O2..mg.L._control, sat_O2.specific = sat_O2_animal - sat_O2_control)

#Extract T0 points and calculate relative O2 values
t_zeros <- join_presens[which(join_presens$Time.point..min.==0),c("Specimen","Species","O2..mg.L._animal", "O2..mg.L._control","O2..mg.L.specific", "sat_O2_control", "sat_O2_animal", "sat_O2.specific")]
dif_presens <- full_join(join_presens, t_zeros, by=c("Species","Specimen"), suffix=c("","T0"))
dif_presens <- mutate(dif_presens, dif_O2.animal = O2..mg.L._animal-O2..mg.L._animalT0, dif_O2.control = O2..mg.L._control-O2..mg.L._controlT0, dif_O2.specific = O2..mg.L.specific-O2..mg.L.specificT0, dif_sat_O2.animal = sat_O2_animal-sat_O2_animalT0, dif_sat_O2.control = sat_O2_control-sat_O2_controlT0, dif_sat_O2.specific = sat_O2.specific-sat_O2.specificT0)

#Correct for zooid length and number
norm_presens <- mutate(dif_presens, corr_O2.animal = dif_O2.animal/((((Zooid.length..mm./2)^2)*Zooid.length..mm.)*Number.of.zooids), corr_O2.control = dif_O2.control, corr_O2.specific = dif_O2.specific/((((Zooid.length..mm./2)^2)*Zooid.length..mm.)*Number.of.zooids), corr_sat_O2.animal = dif_sat_O2.animal/((((Zooid.length..mm./2)^2)*Zooid.length..mm.)*Number.of.zooids), corr_sat_O2.control = dif_sat_O2.control, corr_sat_O2.specific = dif_sat_O2.specific/((((Zooid.length..mm./2)^2)*Zooid.length..mm.)*Number.of.zooids))

#Remove unorthodox measuements
norm_presens <- norm_presens[which(norm_presens$Specimen != "D4-CP-B-1"),] #remove weird C.polae from MGCL2 experiment
norm_presens <- norm_presens[which(norm_presens$corr_O2.specific < 0),] #something has to have gone wrong in these datapoints

#Factorize specimen and species concepts
norm_presens$Specimen <- as.factor(norm_presens$Specimen)
norm_presens$Species <- as.factor(norm_presens$Species)


### PLOTS ###
#Raw O2 plots
rawControl <- ggplot(norm_presens, aes(x=Time.point..min., y=O2..mg.L._control)) + geom_point(aes(col=Species)) + ylab("O2 (mg/L) controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw() + theme(legend.position = "none")
rawAnimal <- ggplot(norm_presens, aes(x=Time.point..min., y=O2..mg.L._animal)) + geom_point(aes(col=Species)) + ylab("O2 (mg/L) animals") + geom_line(aes(col = Species, group=Specimen)) + theme_bw() + theme(legend.position = "none")
rawSpecific <- ggplot(norm_presens, aes(x=Time.point..min., y=O2..mg.L.specific)) + geom_point(aes(col=Species)) + ylab("O2 (mg/L) animals-controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()

pdf("RawO2.pdf", height=4, width=14)
wrap_plots(rawControl,rawAnimal,rawSpecific)
dev.off()

#Raw SAT O2 plots
rawSatControl <- ggplot(norm_presens, aes(x=Time.point..min., y=sat_O2_control)) + geom_point(aes(col=Species)) + ylab("Sat O2 (%) Controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw() + theme(legend.position = "none")
rawSatAnimal <- ggplot(norm_presens, aes(x=Time.point..min., y=sat_O2_animal)) + geom_point(aes(col=Species)) + ylab("Sat O2 (%) (mg/L) Animals") + geom_line(aes(col = Species, group=Specimen)) + theme_bw() + theme(legend.position = "none")
rawSatSpecific <- ggplot(norm_presens, aes(x=Time.point..min., y=sat_O2.specific)) + geom_point(aes(col=Species)) + ylab("Sat O2 (%) Animals-Controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()

pdf("RawSatO2.pdf", height=4, width=14)
wrap_plots(rawSatControl,rawSatAnimal,rawSatSpecific)
dev.off()

#Relative O2 from T-0 plots
RelControl <- ggplot(norm_presens, aes(x=Time.point..min., y=dif_O2.control)) + geom_point(aes(col=Species)) + ylab("O2 (mg/L) from T0 - Controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw() + theme(legend.position = "none")
RelAnimal <- ggplot(norm_presens, aes(x=Time.point..min., y=dif_O2.animal)) + geom_point(aes(col=Species))+ ylab("O2 (mg/L) from T0 - Animals") + geom_line(aes(col = Species, group=Specimen)) + theme_bw() + theme(legend.position = "none")
RelSpecific <- ggplot(norm_presens, aes(x=Time.point..min., y=dif_O2.specific)) + geom_point(aes(col=Species))+ ylab("O2 (mg/L) from T0 - Animals-Controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()

pdf("RelO2.pdf", height=4, width=14)
wrap_plots(RelControl,RelAnimal,RelSpecific)
dev.off()

#Relative SAT O2 T-0 plots
RelSatControl <- ggplot(norm_presens, aes(x=Time.point..min., y=dif_sat_O2.control)) + geom_point(aes(col=Species)) + ylab("Sat O2 (%) from T0 - Controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw() + theme(legend.position = "none")
RelSatAnimal <- ggplot(norm_presens, aes(x=Time.point..min., y=dif_sat_O2.animal)) + geom_point(aes(col=Species))+ ylab("Sat O2 (%) from T0 - Animals") + geom_line(aes(col = Species, group=Specimen)) + theme_bw() + theme(legend.position = "none")
RelSatSpecific <- ggplot(norm_presens, aes(x=Time.point..min., y=dif_sat_O2.specific)) + geom_point(aes(col=Species))+ ylab("Sat O2 (%) from T0 - Animals-Controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()

wrap_plots(RelSatControl,RelSatAnimal,RelSatSpecific)

#Zooid N & Size corrected O2 plots
CorrControl <- ggplot(norm_presens, aes(x=Time.point..min., y=corr_O2.control)) + geom_point(aes(col=Species)) + ylab("Corrected O2 (mg/L) from T0 - Controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw() + theme(legend.position = "none")
CorrAnimal <- ggplot(norm_presens, aes(x=Time.point..min., y=corr_O2.animal)) + geom_point(aes(col=Species))+ ylab("Corrected O2 (mg/L) from T0 - Animals") + geom_line(aes(col = Species, group=Specimen)) + theme_bw() + theme(legend.position = "none")
CorrSpecific <- ggplot(norm_presens, aes(x=Time.point..min., y=corr_O2.specific)) + geom_point(aes(col=Species)) + ylab("Corrected O2 (mg/L) from T0 - Animals-Controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()

wrap_plots(CorrControl,CorrAnimal,CorrSpecific)

#Zooid N & Size corrected O2 % saturation
CorrSatControl <- ggplot(norm_presens, aes(x=Time.point..min., y=corr_sat_O2.control)) + geom_point(aes(col=Species)) + ylab("Corrected Sat O2 (%) from T0 - Controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw() + theme(legend.position = "none")
CorrSatAnimal <- ggplot(norm_presens, aes(x=Time.point..min., y=corr_sat_O2.animal)) + geom_point(aes(col=Species))+ ylab("Corrected Sat O2 (%) from T0 - Animals") + geom_line(aes(col = Species, group=Specimen)) + theme_bw() + theme(legend.position = "none")
CorrSatSpecific <- ggplot(norm_presens, aes(x=Time.point..min., y=corr_sat_O2.specific)) + geom_point(aes(col=Species)) + ylab("Corrected Sat O2 (%) from T0 - Animals-Controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()

wrap_plots(CorrSatControl,CorrSatAnimal,CorrSatSpecific)

#Estimate slopes
slopes <- as.data.frame(matrix(ncol=3,nrow=length(unique(norm_presens$Specimen))))
names(slopes) <- c("Species","Specimen","Slope")
for(i in 1:length(unique(norm_presens$Specimen))){
  spm_i <- unique(norm_presens$Specimen)[i]
  series_i <- norm_presens[which(norm_presens$Specimen==spm_i),c("Specimen","Species","Time.point..min.","corr_sat_O2.specific")]
  lm_i <- lm(Time.point..min.~corr_sat_O2.specific, series_i)
  print(series_i$Specimen[1] %>% as.character());print(series_i$Species[1] %>% as.character());print(lm_i$coefficients)
  slopes[i,1] <- series_i$Species[1] %>% as.character()
  slopes[i,2] <- series_i$Specimen[1] %>% as.character()
  slopes[i,3] <- lm_i$coefficients[2] %>% as.numeric()
}

ggplot(slopes,aes(x=Species,y=Slope))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90))
## GROUP BY Specimen ##

