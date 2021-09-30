library(tidyverse)
library(wql)
library(patchwork)
setwd("~/Documents/salp_ecomorphology/")

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

pdf("RelSatO2.pdf", height=4, width=14)
wrap_plots(RelSatControl,RelSatAnimal,RelSatSpecific)
dev.off()

#Zooid N & Size corrected O2 plots
CorrControl <- ggplot(norm_presens, aes(x=Time.point..min., y=corr_O2.control)) + geom_point(aes(col=Species)) + ylab("Corrected O2 (mg/L) from T0 - Controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw() + theme(legend.position = "none")
CorrAnimal <- ggplot(norm_presens, aes(x=Time.point..min., y=corr_O2.animal)) + geom_point(aes(col=Species))+ ylab("Corrected O2 (mg/L) from T0 - Animals") + geom_line(aes(col = Species, group=Specimen)) + theme_bw() + theme(legend.position = "none")
CorrSpecific <- ggplot(norm_presens, aes(x=Time.point..min., y=corr_O2.specific)) + geom_point(aes(col=Species)) + ylab("Corrected O2 (mg/L) from T0 - Animals-Controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()

pdf("CorrO2.pdf", height=4, width=14)
wrap_plots(CorrControl,CorrAnimal,CorrSpecific)
dev.off()

#Zooid N & Size corrected O2 % saturation
CorrSatControl <- ggplot(norm_presens, aes(x=Time.point..min., y=corr_sat_O2.control)) + geom_point(aes(col=Species)) + ylab("Corrected Sat O2 (%) from T0 - Controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw() + theme(legend.position = "none")
CorrSatAnimal <- ggplot(norm_presens, aes(x=Time.point..min., y=corr_sat_O2.animal)) + geom_point(aes(col=Species))+ ylab("Corrected Sat O2 (%) from T0 - Animals") + geom_line(aes(col = Species, group=Specimen)) + theme_bw() + theme(legend.position = "none")
CorrSatSpecific <- ggplot(norm_presens, aes(x=Time.point..min., y=corr_sat_O2.specific)) + geom_point(aes(col=Species)) + ylab("Corrected Sat O2 (%) from T0 - Animals-Controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()

pdf("CorrSatO2.pdf", height=4, width=14)
wrap_plots(CorrSatControl,CorrSatAnimal,CorrSatSpecific)
dev.off()

#Estimate slopes
slopes <- as.data.frame(matrix(ncol=8,nrow=length(unique(norm_presens$Specimen))))
names(slopes) <- c("Species","Specimen","Zooid.length..mm.","Number.of.zooids","Slope_O2","Slope_CorrO2sat","Timespan","Temperature...C.")
for(i in 1:length(unique(norm_presens$Specimen))){
  spm_i <- unique(norm_presens$Specimen)[i]
  series_i <- norm_presens[which(norm_presens$Specimen==spm_i),c("Specimen","Species","Zooid.length..mm.","Number.of.zooids","Time.point..min.","Temperature...C.","sat_O2.specific","corr_sat_O2.specific")]
  lm_i <- lm(sat_O2.specific~Time.point..min., series_i)
  lm_j <- lm(corr_sat_O2.specific~Time.point..min., series_i)
  print(series_i$Specimen[1] %>% as.character());print(series_i$Species[1] %>% as.character());print(lm_i$coefficients);print(lm_j$coefficients)
  time.span_i <- max(series_i$Time.point..min.)-min(series_i$Time.point..min.)
  slopes[i,1] <- series_i$Species[1] %>% as.character()
  slopes[i,2] <- series_i$Specimen[1] %>% as.character()
  slopes[i,3] <- series_i$Zooid.length..mm.[1] %>% as.character()
  slopes[i,4] <- series_i$Number.of.zooids[1] %>% as.character()
  slopes[i,5] <- lm_i$coefficients[2] %>% as.numeric()
  slopes[i,6] <- lm_j$coefficients[2] %>% as.numeric()
  slopes[i,7] <- time.span_i
  slopes[i,8] <- series_i$Temperature...C. %>% mean()
}

pdf("slopes_SatO2.pdf", height=6, width=10)
ggplot(slopes,aes(x=Species,y=-Slope_O2))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90))
dev.off()

pdf("slopes_CorrSatO2.pdf", height=6, width=10)
ggplot(slopes,aes(x=Species,y=-Slope_CorrO2sat))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90))
dev.off()

## GROUP BY Specimen ##

### COST OF TRANSPORT ###
swim <- data.frame(Species=c("Pegea confoederata","Iasis (Weelia) cylindrica","Cyclosalpa affinis"),Speed.cm.s=c(1.68,1.23,1.72))
COT <- left_join(slopes, swim, by="Species") 
COT$Zooid.length..mm. <- as.numeric(COT$Zooid.length..mm.)
COT$Number.of.zooids <- as.numeric(COT$Number.of.zooids)
COT %>% mutate(Speed.body.s = Zooid.length..mm.*Speed.cm.s/10) -> COT
COT %>% mutate(absO2slope = abs(Slope_O2)*0.225*oxySol(Temperature...C., 30.1, 1)) -> COT
#Define COT.abs as mgO2/(Zooid.vol*cm_moved)  v.v. COT.rel as mgO2/(Zooid.vol*bodylength)
COT %>% mutate(COT.abs = absO2slope*Speed.cm.s*60/((((Zooid.length..mm./2)^2)*Zooid.length..mm.)*Number.of.zooids), COT.rel = abs(Slope_O2)*0.225*Speed.body.s*60/((((Zooid.length..mm./2)^2)*Zooid.length..mm.)*Number.of.zooids))->COT
COT <- COT[which(!is.na(COT$COT.abs)),]

pdf("COT_abs.pdf", height=6, width=10)
ggplot(COT,aes(x=Species,y=COT.abs))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90))
dev.off()

pdf("COT_rel.pdf", height=6, width=10)
ggplot(COT,aes(x=Species,y=COT.rel))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90))
dev.off()

#Estimate Carbon content for each species and recalculate carbon-based COT
mm_to_carbon <- read.csv("salpcarbon.tsv", header=T, stringsAsFactors = F, sep='\t')
mm_to_carbon$Species[mm_to_carbon$Species=="Pegea socia"] <- "Pegea confoederata" ##PROXY P. socia for P. confoederata
mm_to_carbon$Species[mm_to_carbon$Species=="Soestia (Iasis) zonaria"] <- "Iasis (Weelia) cylindrica" ##PROXY Soestia for Weelia
COT <- left_join(COT, mm_to_carbon[,c(1,4,5)], by="Species")
COT <- mutate(COT, mgC = Regression_b*Zooid.length..mm.^Regression_alpha)
COT %>% mutate(COT.abs = absO2slope*Speed.cm.s*60/((((Zooid.length..mm./2)^2)*Zooid.length..mm.)*Number.of.zooids), COT.rel = abs(Slope_O2)*0.225*Speed.body.s*60/((((Zooid.length..mm./2)^2)*Zooid.length..mm.)*Number.of.zooids))->COT