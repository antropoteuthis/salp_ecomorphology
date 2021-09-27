library(tidyverse)
setwd("Documents/salp_ecomorphology/")

presens <- read.csv("respirometry_Kona21.tsv", header = T, sep = "\t", stringsAsFactors = F)
presens$Species[which(is.na(presens$Species))] <- "Control"
presens$Specimen[which(is.na(presens$Specimen))] <- "Control"

presens <- presens[which(presens$Container=="Plastic" & presens$Treatment == "Intact"),]
controls <- presens[which(presens$Specimen=="Control"),]
specimens <- presens[which(presens$Specimen != "Control"),]

####
join_presens <- full_join(specimens, controls[-c(2:10,14,16:18)], by=c("Experiment","Time.point..min."), suffix=c("_animal","_control"))
names(join_presens)
join_presens <- mutate(join_presens, O2..mg.L.specific = O2..mg.L._animal-O2..mg.L._control)

t_zeros <- join_presens[which(join_presens$Time.point..min.==0),c("Specimen","Species","O2..mg.L._animal", "O2..mg.L._control","O2..mg.L.specific")]
dif_presens <- full_join(join_presens, t_zeros, by=c("Species","Specimen"), suffix=c("","T0"))
dif_presens <- mutate(dif_presens, dif_O2.animal = O2..mg.L._animal-O2..mg.L._animalT0, dif_O2.control = O2..mg.L._control-O2..mg.L._controlT0, dif_O2.specific = O2..mg.L.specific-O2..mg.L.specificT0)
dif_presens$Temperature...C.[which(is.na(dif_presens$Temperature...C.))]<-29 #### MAKING SHIT UP ####
norm_presens <- mutate(dif_presens, corr_O2.animal = Temperature...C.*dif_O2.animal/((Zooid.length..mm.^3)*Number.of.zooids), corr_O2.control = Temperature...C.*dif_O2.control, corr_O2.specific = Temperature...C.*dif_O2.specific/((Zooid.length..mm.^3)*Number.of.zooids))

norm_presens <- norm_presens[which(norm_presens$Specimen != "D4-CP-B-1"),] #remove weird C.polae from MGCL2 experiment
norm_presens <- norm_presens[which(norm_presens$corr_O2.specific < 0),] #something has to have gone wrong in these datapoints
norm_presens$Specimen <- as.factor(norm_presens$Specimen)
norm_presens$Species <- as.factor(norm_presens$Species)

ggplot(norm_presens, aes(x=Time.point..min., y=O2..mg.L._control)) + geom_point(aes(col=Species)) + ylab("O2 (mg/L) controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()
ggplot(norm_presens, aes(x=Time.point..min., y=O2..mg.L._animal)) + geom_point(aes(col=Species)) + ylab("O2 (mg/L) animals") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()
ggplot(norm_presens, aes(x=Time.point..min., y=O2..mg.L.specific)) + geom_point(aes(col=Species)) + ylab("O2 (mg/L) animals-controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()
ggplot(norm_presens, aes(x=Time.point..min., y=dif_O2.control)) + geom_point(aes(col=Species)) + ylab("O2 (mg/L) from T0 - Controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()
ggplot(norm_presens, aes(x=Time.point..min., y=dif_O2.animal)) + geom_point(aes(col=Species))+ ylab("O2 (mg/L) from T0 - Animals") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()
ggplot(norm_presens, aes(x=Time.point..min., y=dif_O2.specific)) + geom_point(aes(col=Species))+ ylab("O2 (mg/L) from T0 - Animals-Controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()
ggplot(norm_presens, aes(x=Time.point..min., y=corr_O2.control)) + geom_point(aes(col=Species)) + ylab("Corrected O2 (mg/L) from T0 - Controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()
ggplot(norm_presens, aes(x=Time.point..min., y=corr_O2.animal)) + geom_point(aes(col=Species))+ ylab("Corrected O2 (mg/L) from T0 - Animals") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()
ggplot(norm_presens, aes(x=Time.point..min., y=corr_O2.specific)) + geom_point(aes(col=Species)) + ylab("Corrected O2 (mg/L) from T0 - Controls") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()

slopes <- as.data.frame(matrix(ncol=3,nrow=length(unique(norm_presens$Specimen))))
names(slopes) <- c("Species","Specimen","Slope")
for(i in 1:length(unique(norm_presens$Specimen))){
  spm_i <- unique(norm_presens$Specimen)[i]
  series_i <- norm_presens[which(norm_presens$Specimen==spm_i),c("Specimen","Species","Time.point..min.","corr_O2.specific")]
  lm_i <- lm(Time.point..min.~corr_O2.specific, series_i)
  print(series_i$Specimen[1] %>% as.character());print(series_i$Species[1] %>% as.character());print(lm_i$coefficients)
  slopes[i,1] <- series_i$Species[1] %>% as.character()
  slopes[i,2] <- series_i$Specimen[1] %>% as.character()
  slopes[i,3] <- lm_i$coefficients[2] %>% as.numeric()
}

ggplot(slopes,aes(x=Species,y=Slope))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90))
## GROUP BY Specimen ##

