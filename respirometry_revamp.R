library(tidyverse)
library(wql)
library(patchwork)
require(data.table)
library(mgcv)
setwd("~/Documents/salp_ecomorphology/")

#Load data and label control rows
presens <- read.csv("respirometry_Kona22.tsv", header = T, sep = "\t", stringsAsFactors = F)
presens$Species[which(is.na(presens$Species))] <- "Control"
presens$Specimen[which(is.na(presens$Specimen))] <- "Control"

#Filter by plastic and intact treatments
presens <- presens[which(presens$Container=="Plastic" & presens$Treatment == "Intact"),]

#Filter by blastozooids
presens <- presens[which(presens$Stage != "Oozoid" & presens$Stage != "Oozoid+Stolon" | is.na(presens$Stage)),]

#Fill in contained volume
presens$Container.volume..ml.[which(is.na(presens$Container.volume..ml.) & presens$Sensor.ID %in% c("P1","P2","P4","P5"))] <- 170
presens$Container.volume..ml.[which(is.na(presens$Container.volume..ml.) & presens$Sensor.ID == "P3 control")] <- 72

#Estimate absolute oxygen mg
ggplot(presens,aes(x=Specimen,y=O2..mg.L.))+geom_point(aes(col=Time.point..min.))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
presens <- mutate(presens, abs_O2.mg. = O2..mg.L.*Container.volume..ml./1000)
ggplot(presens,aes(x=Specimen,y=abs_O2.mg.))+geom_point(aes(col=Container.volume..ml. %>% log()))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Correct for O2 saturation limit wit temperature
#presens$Temperature...C.[which(is.na(presens$Temperature...C.))]<-29 #### MAKING SHIT UP ####
#presens <- mutate(presens, sat_O2 = O2..mg.L./oxySol(presens$Temperature...C., 30.1, 1))

#Subset controls and specimen measures
controls <- presens[which(presens$Specimen=="Control"),]
specimens <- presens[which(presens$Specimen != "Control"),]

#Keep only controls with the largest jar volume and in case of a tie, those that report the highest O2 concentration (purportedly those with least organic matter decay).
controls <- controls %>% 
  group_by(Experiment, Time.point..min.) %>% 
  filter(Container.volume..ml. == max(Container.volume..ml.) & O2..mg.L. == max(O2..mg.L.)) %>% .[,-which(names(.)=="Sensor.ID" | names(.)=="Measurement")] 

controls=unique(controls)

#Extract T0 points and calculate O2 consumed
# t_zeros_control <- controls[which(controls$Time.point..min.==0),c("Experiment","abs_O2.mg.")]
# dif_controls <- full_join(controls, t_zeros_control, by=c("Experiment"), suffix=c("","T0"))
# dif_controls <- mutate(dif_controls, dif_O2.mg. = abs_O2.mg.-abs_O2.mg.T0)
# 
# t_zeros_specimens <- specimens[which(specimens$Time.point..min.==0),c("Specimen","Experiment","abs_O2.mg.")]
# dif_specimens <- full_join(specimens, t_zeros_specimens, by=c("Specimen","Experiment"), suffix=c("","T0"))
# dif_specimens <- mutate(dif_specimens, dif_O2.mg. = abs_O2.mg.-abs_O2.mg.T0)

#Join Specimens+Controls and subtract correcting for the difference in volume for specific O2 measurements
join_presens <- full_join(specimens, controls[c("Date","Experiment","Time.point..min.","abs_O2.mg.","Container.volume..ml.")], 
                            by=c("Experiment","Time.point..min.","Date"), suffix=c("_animal","_control"))
join_presens <- mutate(join_presens, abs_O2.mg.specific = abs_O2.mg._animal-(abs_O2.mg._control*(Container.volume..ml._animal/Container.volume..ml._control)))
join_presens <- join_presens[which(!is.na(join_presens$Measurement)),]

ggplot(join_presens,aes(x=Specimen,y=abs_O2.mg.specific))+geom_point(aes(col=Container.volume..ml._animal %>% log()))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

join_presens %>% 
  mutate(exbool = as.character(Experiment)) %>% 
ggplot(aes(x=Time.point..min., y=abs_O2.mg.specific)) + 
  geom_point(aes(col= Species)) + ylab("O2 (mg)") + 
  geom_line(aes(col = Species, group=Specimen)) + 
  theme_bw()

#Colony volume imputation

imputed_vols <- data.frame(imputed_vol = round(join_presens[which(!is.na(join_presens$Colony.volume..ml.)),"Number.of.zooids"]*pi*join_presens[which(!is.na(join_presens$Colony.volume..ml.)),"Zooid.length..mm."]*((join_presens[which(!is.na(join_presens$Colony.volume..ml.)),"Zooid.length..mm."]/2)^2)*0.001,1),
                           real_vol = join_presens[which(!is.na(join_presens$Colony.volume..ml.)),"Colony.volume..ml."], 
                           Number.of.zooids = join_presens[which(!is.na(join_presens$Colony.volume..ml.)), "Number.of.zooids"],
                           Zooid.length..mm. = join_presens[which(!is.na(join_presens$Colony.volume..ml.)),"Zooid.length..mm."],
                           Species = join_presens[which(!is.na(join_presens$Colony.volume..ml.)),"Species"])

mutate(imputed_vols, estimate_vol=Number.of.zooids*(0.00015*pi*Zooid.length..mm.*((0.3*Zooid.length..mm.)^2 - ((0.2*Zooid.length..mm.)^2)))) %>%
  ggplot(aes(x=estimate_vol, y=real_vol)) + 
  geom_point(aes(col=Species))

fit3 = gam(real_vol~Zooid.length..mm.+Number.of.zooids, data = imputed_vols)

unknown=as.data.frame(join_presens[which(is.na(join_presens$Colony.volume..ml.)),c("Number.of.zooids","Zooid.length..mm.")])
join_presens$Colony.volume..ml.[which(is.na(join_presens$Colony.volume..ml.))] <- predict.gam(fit3, unknown) %>% as.vector()

norm_presens <- join_presens

#Remove unorthodox measuements
#norm_presens <- join_presens[which(join_presens$Specimen != "D4-CP-B-1"),] #remove weird C.polae from MGCL2 experiment
#norm_presens <- norm_presens[which(norm_presens$Specimen != "D25-Pso-B-1"),] #weird Pegea socia o2 increase
#norm_presens <- norm_presens[which(norm_presens$corr_O2.specific < 0),] #something has to have gone wrong in these datapoints

#Factorize specimen and species concepts
norm_presens$Specimen <- as.factor(norm_presens$Specimen)
norm_presens$Species <- as.factor(norm_presens$Species)
norm_presens$Experiment <- as.factor(norm_presens$Experiment)

### PLOTS ###
#Raw O2 plots
rawControl <- ggplot(norm_presens, aes(x=Time.point..min., y=abs_O2.mg._control)) + geom_point(aes(cex=log(Container.volume..ml._control))) + ylab("O2 (mg) Controls") + geom_line(aes(group=Specimen)) + theme_bw()
rawAnimal <- ggplot(norm_presens, aes(x=Time.point..min., y=abs_O2.mg._animal)) + geom_point(aes(col=Species)) + ylab("O2 (mg) Animals") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()+ theme(legend.position = "none")
rawDifference <- ggplot(norm_presens, aes(x=Time.point..min., y=abs_O2.mg.specific)) + geom_point(aes(col=Species)) + ylab("O2 (mg) Difference") + geom_line(aes(col = Species, group=Specimen)) + theme_bw()

pdf("RawO2.pdf", height=4, width=14)
wrap_plots(rawControl,rawAnimal,rawDifference)
dev.off()

#Estimate slopes
slopes <- as.data.frame(matrix(ncol=8,nrow=length(unique(norm_presens$Specimen))))
names(slopes) <- c("Species","Specimen","Colony.volume..ml.","Zooid.length..mm.","Number.of.zooids","Slope_O2","Timespan","Temperature...C.")
for(i in 1:length(unique(norm_presens$Specimen))){
  spm_i <- unique(norm_presens$Specimen)[i]
  series_i <- norm_presens[which(norm_presens$Specimen==spm_i),c("Specimen","Species","Zooid.length..mm.","Number.of.zooids","Colony.volume..ml.","Time.point..min.","Temperature...C.","abs_O2.mg.specific")]
  lm_i <- lm(abs_O2.mg.specific~Time.point..min., series_i)
  print(series_i$Specimen[1] %>% as.character());print(series_i$Species[1] %>% as.character());print(lm_i$coefficients)
  time.span_i <- max(series_i$Time.point..min.)-min(series_i$Time.point..min.)
  slopes[i,1] <- series_i$Species[1] %>% as.character()
  slopes[i,2] <- series_i$Specimen[1] %>% as.character()
  slopes[i,3] <- series_i$Colony.volume..ml.[1] %>% as.numeric()
  slopes[i,4] <- series_i$Zooid.length..mm.[1] %>% as.numeric()
  slopes[i,5] <- series_i$Number.of.zooids[1] %>% as.numeric()
  slopes[i,6] <- lm_i$coefficients[2] %>% as.numeric()
  slopes[i,7] <- time.span_i
  slopes[i,8] <- series_i$Temperature...C. %>% mean()
}

pdf("slopes_SatO2.pdf", height=6, width=10)
ggplot(slopes,aes(x=Species,y=-Slope_O2))+geom_boxplot()+theme_bw()+theme(axis.text.x = element_text(angle = 90))
dev.off()

slopes %>% mutate(Slope_normalized = Slope_O2/Colony.volume..ml.) -> slopes

pdf("slopes_Corrected.pdf", height=6, width=10)
ggplot(slopes,aes(x=Species,y=-Slope_normalized))+geom_boxplot()+ylab("-Slope / Specimen biovolume (ml)")+theme_bw()+theme(axis.text.x = element_text(angle = 90))
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
COT_pruned <- COT[which(!is.na(COT$COT.abs)),]

pdf("COT_abs.pdf", height=6, width=10)
ggplot(COT_pruned,aes(x=Species,y=COT.abs))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90))
dev.off()

pdf("COT_rel.pdf", height=6, width=10)
ggplot(COT_pruned,aes(x=Species,y=COT.rel))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90))
dev.off()

#Estimate Carbon content for each species and recalculate carbon-based COT
mm_to_carbon <- read.csv("Madin1981_salpcarbon.tsv", header=T, stringsAsFactors = F, sep='\t')
mm_to_carbon$Generation[mm_to_carbon$Species=="Iasis (Weelia) cylindrica"] <- "a" ##PROXY 
COT <- left_join(COT, mm_to_carbon[,c(1,4,5)], by="Species")
COT <- mutate(COT, mgC = Regression_b*Zooid.length..mm.^Regression_alpha)
COT <- mutate(COT, COT_C.abs = absO2slope*Speed.cm.s*60/mgC, COT_C.rel = abs(Slope_O2)*0.225*Speed.body.s*60/mgC)

pdf("COT_C_abs.pdf", height=6, width=10)
ggplot(COT,aes(x=Species,y=COT_C.abs))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90))
dev.off()

pdf("COT_C_rel.pdf", height=6, width=10)
ggplot(COT,aes(x=Species,y=COT_C.rel))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90))
dev.off()
