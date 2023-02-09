############ LOAD UP AND SET UP ##############

library(tidyverse)
library(wql)
library(ggrepel)
library(patchwork)
require(data.table)
library(mgcv)
setwd("~/salp_ecomorphology/")

#Load data and label: control rows, paired specimens
presens <- read.csv("respirometry_Kona21-22s.tsv", header = T, sep = "\t", stringsAsFactors = F)
presens$Species[which(is.na(presens$Species))] <- "Control"
presens$Specimen[which(is.na(presens$Specimen))] <- "Control"
presens <- mutate(presens, is.paired=ifelse(Measurement.notes=="Paired", "Yes", "No"))
as.numeric(as.factor(presens$Injection.time)) -> presens$Injection.time 

#Filter to get only plastic and non-mgcl2 blastozooid measurements
presens <- presens[which(presens$Container=="Plastic" & presens$Treatment != "MgCl2"),]
presens <- presens[which(presens$Stage=="Blastozooid" | is.na(presens$Stage)),]

#Estimate absolute oxygen mg
ggplot(presens,aes(x=Specimen,y=O2..mg.L.))+geom_point(aes(col=Time.point..min.))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

presens <- mutate(presens, abs_O2.mg. = O2..mg.L.*Container.volume..ml./1000)

ggplot(presens,aes(x=Specimen,y=abs_O2.mg.))+geom_point(aes(col=Container.volume..ml. %>% log()))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Correct for O2 saturation limit wit temperature
#presens$Temperature...C.[which(is.na(presens$Temperature...C.))]<-29 #### MAKING SHIT UP ####
#presens <- mutate(presens, sat_O2 = O2..mg.L./oxySol(presens$Temperature...C., 30.1, 1))

#Colony volume imputation for those specimens where we didn't measure the specimen biovolume in the field
imputed_vols <- data.frame(imputed_vol = round(presens[which(!is.na(presens$Colony.volume..ml.)),"Number.of.zooids"]*pi*presens[which(!is.na(presens$Colony.volume..ml.)),"Zooid.length..mm."]*((presens[which(!is.na(presens$Colony.volume..ml.)),"Zooid.length..mm."]/2)^2)*0.001,1),
                           real_vol = presens[which(!is.na(presens$Colony.volume..ml.)),"Colony.volume..ml."], 
                           Number.of.zooids = presens[which(!is.na(presens$Colony.volume..ml.)), "Number.of.zooids"],
                           Zooid.length..mm. = presens[which(!is.na(presens$Colony.volume..ml.)),"Zooid.length..mm."],
                           Species = presens[which(!is.na(presens$Colony.volume..ml.)),"Species"]) 

  #try out an imputation fit by hand using a 3D ellipsoid formula
mutate(imputed_vols, estimate_vol=Number.of.zooids*(ifelse(Species=="Salpa maxima",0.00005,0.00015)*pi*Zooid.length..mm.*((0.35*Zooid.length..mm.)^2 - ((0.25*Zooid.length..mm.)^2)))) %>%
  ggplot(aes(x=real_vol, y=estimate_vol)) + 
  geom_point(aes(col=Species))

  #try a GAM-based modeling for the imputation
fit3 = gam(real_vol~Zooid.length..mm.+Number.of.zooids, data = imputed_vols)
fit3$coefficients[1] <- 0
unknown=as.data.frame(presens[which(is.na(presens$Colony.volume..ml.) & presens$Species != "Control"),c("Number.of.zooids","Zooid.length..mm.")])
presens$Colony.volume..ml.[which(is.na(presens$Colony.volume..ml.) & presens$Species != "Control")] <- predict.gam(fit3, unknown) %>% as.vector()

#Sensitivity of container effects by species (needs more development)
ggplot(presens[which(presens$Time.point..min.>115),],aes(x=Sensor.ID,y=O2..mg.L.))+geom_point(aes(col=Species))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(presens[which(presens$Time.point..min.>115),],aes(x=Sensor.ID,y=abs_O2.mg.))+geom_point(aes(col=Species))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

########### ANALYSES ##############

#Subset controls and specimen measures for swimmers, anesthetized, and paired experiments

#### SOME CONTROLS' TREATMENTS DONT MATCH THE TREATMENT OF THEIR SPECIMENS!!!  ####

# controls_swim <- presens[which(presens$Specimen=="Control" & presens$Treatment == "Intact" & !str_detect(presens$Measurement.notes, "Paired")),] %>% 
#   group_by(Experiment, Time.point..min.) %>% 
#   filter(Container.volume..ml. == max(Container.volume..ml.) & O2..mg.L. == max(O2..mg.L.)) %>% .[,-which(names(.)=="Sensor.ID" | names(.)=="Measurement")] %>% unique()
# controls_ko <-  presens[which(presens$Specimen=="Control" & presens$Treatment == "Anesthetized" & !str_detect(presens$Measurement.notes, "Paired")),] %>% 
#   group_by(Experiment, Time.point..min.) %>% 
#   filter(Container.volume..ml. == max(Container.volume..ml.) & O2..mg.L. == max(O2..mg.L.)) %>% .[,-which(names(.)=="Sensor.ID" | names(.)=="Measurement")] %>% unique()
# controls_swim_paired <- presens[which(presens$Specimen=="Control" & presens$Treatment == "Intact" & str_detect(presens$Measurement.notes, "Paired")),] %>% 
#   group_by(Experiment, Time.point..min.) %>% 
#   filter(Container.volume..ml. == max(Container.volume..ml.) & O2..mg.L. == max(O2..mg.L.)) %>% .[,-which(names(.)=="Sensor.ID" | names(.)=="Measurement")] %>% unique()
# controls_ko_paired <- presens[which(presens$Specimen=="Control" & presens$Treatment == "Anesthetized" & str_detect(presens$Measurement.notes, "Paired")),] %>% 
#   group_by(Experiment, Time.point..min.) %>% 
#   filter(Container.volume..ml. == max(Container.volume..ml.) & O2..mg.L. == max(O2..mg.L.)) %>% .[,-which(names(.)=="Sensor.ID" | names(.)=="Measurement")] %>% unique()
# 
# controls <- presens[which(presens$Specimen=="Control"),] %>%  group_by(Experiment, Time.point..min.) %>%  filter(Container.volume..ml. == max(Container.volume..ml.) & O2..mg.L. == max(O2..mg.L.)) %>% .[,-which(names(.)=="Sensor.ID" | names(.)=="Measurement")] %>% unique()

  #Specimens
# specimens_swim <- presens[which(presens$Specimen != "Control" & presens$Treatment == "Intact" & !str_detect(presens$Measurement.notes, "Paired")),]
# specimens_ko <- presens[which(presens$Specimen != "Control" & presens$Treatment == "Anesthetized" & !str_detect(presens$Measurement.notes, "Paired")),]
# specimens_swim_paired <- presens[which(presens$Specimen != "Control" & presens$Treatment == "Intact" & str_detect(presens$Measurement.notes, "Paired")),]
# specimens_ko_paired <- presens[which(presens$Specimen != "Control" & presens$Treatment == "Anesthetized" & str_detect(presens$Measurement.notes, "Paired")),]
# 
# specimens <- presens[which(presens$Specimen != "Control"),]


#Extract T0 points and calculate O2 consumed
# t_zeros_control <- controls[which(controls$Time.point..min.==0),c("Experiment","abs_O2.mg.")]
# dif_controls <- full_join(controls, t_zeros_control, by=c("Experiment"), suffix=c("","T0"))
# dif_controls <- mutate(dif_controls, dif_O2.mg. = abs_O2.mg.-abs_O2.mg.T0)
# 
# t_zeros_specimens <- specimens[which(specimens$Time.point..min.==0),c("Specimen","Experiment","abs_O2.mg.")]
# dif_specimens <- full_join(specimens, t_zeros_specimens, by=c("Specimen","Experiment"), suffix=c("","T0"))
# dif_specimens <- mutate(dif_specimens, dif_O2.mg. = abs_O2.mg.-abs_O2.mg.T0)

#Join Specimens+Controls and subtract correcting for the difference in volume for specific O2 measurements

#Match each measurement with its relevant control
# presens <- mutate(presens, abs_O2.mg.specific = abs_O2.mg._animal-(abs_O2.mg._control*(Container.volume..ml._animal/Container.volume..ml._control)))
# presens$Injection.time[join_presens$Injection.time==""]<-NA
# presens$Start.time[join_presens$Start.time==""]<-NA
# ggplot(join_presens,aes(x=Specimen,y=abs_O2.mg.specific)) +
#   geom_point(aes(col=Container.volume..ml. %>% log())) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# join_presens %>% 
#   mutate(exbool = as.character(Experiment)) %>% 
#   group_by(Treatment) %>% 
#   ggplot(aes(x=Time.point..min., y=abs_O2.mg.specific)) + 
#   geom_point(aes(col= Species)) + ylab("O2 (mg)") + 
#   geom_line(aes(col = Species, group=Specimen)) + 
#   theme_bw() +
#   ylim(c(-5,0))+
#   facet_wrap(~Treatment)

# join_presens %>% filter(is.paired=="Yes") %>% 
#   mutate(exbool = as.character(Experiment)) %>% 
#   ggplot(aes(x=Time.point..min., y=abs_O2.mg.specific)) + 
#   geom_point(aes(col=Species, alpha=Treatment))+ ylab("O2 (mg)") + 
#   geom_line(aes(col = Species, group=Specimen, alpha=Treatment)) + 
#   theme_bw()
# 
#   #Swimmers
# join_swim <- full_join(specimens_swim, controls_swim[c("Date","Experiment","Time.point..min.","abs_O2.mg.","Container.volume..ml.")], 
#                             by=c("Experiment","Time.point..min.","Date"), suffix=c("_animal","_control"))
# join_swim <- mutate(join_swim, abs_O2.mg.specific = abs_O2.mg._animal-(abs_O2.mg._control*(Container.volume..ml._animal/Container.volume..ml._control)))
# join_swim <- join_swim[which(!is.na(join_swim$Measurement)),]
# 
# ggplot(join_swim,aes(x=Specimen,y=abs_O2.mg.specific))+geom_point(aes(col=Container.volume..ml._animal %>% log()))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# join_swim %>% 
#   mutate(exbool = as.character(Experiment)) %>% 
# ggplot(aes(x=Time.point..min., y=abs_O2.mg.specific)) + 
#   geom_point(aes(col= Species)) + ylab("O2 (mg)") + 
#   geom_line(aes(col = Species, group=Specimen)) + 
#   theme_bw()
# 
#   #Anesthetized
# join_ko <- full_join(specimens_ko, controls_ko[c("Date","Experiment","Time.point..min.","abs_O2.mg.","Container.volume..ml.")], 
#                      by=c("Experiment","Time.point..min.","Date"), suffix=c("_animal","_control"))
# join_ko <- mutate(join_ko, abs_O2.mg.specific = abs_O2.mg._animal-(abs_O2.mg._control*(Container.volume..ml._animal/Container.volume..ml._control)))
# join_ko <- join_ko[which(!is.na(join_ko$Measurement)),]
# 
# ggplot(join_ko,aes(x=Specimen,y=abs_O2.mg.specific))+geom_point(aes(col=Container.volume..ml._animal %>% log()))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# join_ko %>% 
#   mutate(exbool = as.character(Experiment)) %>% 
#   ggplot(aes(x=Time.point..min., y=abs_O2.mg.specific)) + 
#   geom_point(aes(col= Species)) + ylab("O2 (mg)") + 
#   geom_line(aes(col = Species, group=Specimen)) + 
#   theme_bw()
# 
#   #Paired
# join_swim_paired <- full_join(specimens_swim_paired, controls_swim_paired[c("Date","Experiment","Time.point..min.","abs_O2.mg.","Container.volume..ml.")], 
#                               by=c("Experiment","Time.point..min.","Date"), suffix=c("_swim_animal","_swim_control"))
# join_swim_paired <- mutate(join_swim_paired, abs_O2.mg._swim_specific = abs_O2.mg._swim_animal-(abs_O2.mg._swim_control*(Container.volume..ml._swim_animal/Container.volume..ml._swim_control)))
# join_swim_paired <- join_swim_paired[which(!is.na(join_swim_paired$Measurement)),]
# 
# ggplot(join_swim_paired,aes(x=Specimen,y=abs_O2.mg._swim_specific))+geom_point(aes(col=Container.volume..ml._swim_animal %>% log()))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# join_swim_paired %>% 
#   mutate(exbool = as.character(Experiment)) %>% 
#   ggplot(aes(x=Time.point..min., y=abs_O2.mg._swim_specific)) + 
#   geom_point(aes(col= Species)) + ylab("O2 (mg)") + 
#   geom_line(aes(col = Species, group=Specimen)) + 
#   theme_bw()
# 
# join_ko_paired <- full_join(specimens_ko_paired, controls_ko_paired[c("Date","Experiment","Time.point..min.","abs_O2.mg.","Container.volume..ml.")], 
#                               by=c("Experiment","Time.point..min.","Date"), suffix=c("_ko_animal","_ko_control"))
# join_ko_paired <- mutate(join_ko_paired, abs_O2.mg._ko_specific = abs_O2.mg._ko_animal-(abs_O2.mg._ko_control*(Container.volume..ml._ko_animal/Container.volume..ml._ko_control)))
# join_ko_paired <- join_ko_paired[which(!is.na(join_ko_paired$Measurement)),]
# 
# ggplot(join_ko_paired,aes(x=Specimen,y=abs_O2.mg._ko_specific))+geom_point(aes(col=Container.volume..ml._ko_animal %>% log()))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# join_ko_paired %>% 
#   mutate(exbool = as.character(Experiment)) %>% 
#   ggplot(aes(x=Time.point..min., y=abs_O2.mg._ko_specific)) + 
#   geom_point(aes(col= Species)) + ylab("O2 (mg)") + 
#   geom_line(aes(col = Species, group=Specimen)) + 
#   theme_bw()

#######

#Match each measurement with its relevant control
presens$abs_O2.mg._control <- presens$Measurement #placeholder column
for(i in 1:nrow(presens)){
  if(!is.na(presens$Relevant.control[i]) & presens$Specimen[i] != "Control"){
    presens$abs_O2.mg._control[i] <- presens$abs_O2.mg.[which(presens$Measurement == presens$Relevant.control[i])]
  }
}

#remove control rows
rename(as_tibble(presens) %>% filter(presens$Specimen != "Control"), abs_O2.mg._animal=abs_O2.mg.) %>% as.data.frame() -> presens

#Remove unorthodox measuements
#norm_presens <- join_presens[which(join_presens$Specimen != "D4-CP-B-1"),] #remove weird C.polae from MGCL2 experiment
#norm_presens <- norm_presens[which(norm_presens$Specimen != "D25-Pso-B-1"),] #weird Pegea socia o2 increase
#norm_presens <- norm_presens[which(norm_presens$corr_O2.specific < 0),] #something has to have gone wrong in these datapoints

#Factorize specimen and species concepts
factSPsp <- function(df){
  df$Species <- as.factor(df$Species)
  df$Specimen <- as.factor(df$Specimen)
  df$Experiment <- as.factor(df$Specimen)
}

factSPsp(presens)
# factSPsp(join_ko)
# factSPsp(join_ko_paired)
# factSPsp(join_swim)
# factSPsp(join_swim_paired)
# factSPsp(join_presens)

### PLOTS ###

# ggplot(join_swim, aes(x=Time.point..min., y=abs_O2.mg.specific)) +
#   geom_point(aes(col=Species)) +
#   ylab("O2 (mg) Difference Intact") +
#   geom_line(aes(col = Species, group=Specimen)) +
#   ylim(c(-2.1,0)) +
#   theme_bw()
# 
# ggplot(join_ko, aes(x=Time.point..min., y=abs_O2.mg.specific)) +
#   geom_point(aes(col=Species)) +
#   ylab("O2 (mg) Difference Intact") +
#   geom_line(aes(col = Species, group=Specimen)) +
#   theme_bw()

#Raw O2 plots for swim
rawControlswim <- ggplot(presens, aes(x=Time.point..min., y=abs_O2.mg._control)) +
  geom_point() +
  ylab("O2 (mg) Controls Intact") +
  theme_bw()
rawAnimalswim <- ggplot(presens, aes(x=Time.point..min., y=abs_O2.mg._animal)) +
  geom_point(aes(col=Species)) +
  ylab("O2 (mg) Animals Intact") +
  geom_line(aes(col = Species, group=Specimen)) +
  theme_bw()

# rawDifferenceSwim <- ggplot(presens %>% 
#                               filter(Treatment=="Intact"), aes(x=Time.point..min., y=abs_O2.mg.specific)) +
#   geom_point(aes(col=Species)) + 
#   ylab("O2 (mg) Difference Intact") + 
#   geom_line(aes(col = Species, group=Specimen)) + 
#   theme_bw()

pdf("Figures_respirometry/Kona2022/RawO2.pdf", height=4, width=14)
wrap_plots(rawControlswim,rawAnimalswim)
dev.off()

#Estimate slopes   ### ALERT WE ARE REMOVING -VE NUMBERS AND PUTTING ZEROES INSTEAD  ####

#SLOPE CALCULATIONS
slopes <- as.data.frame(matrix(ncol=12,nrow=length(unique(paste0(presens$Specimen,presens$Treatment)))))
names(slopes) <- c("Species","Specimen","Colony.volume..ml.","Zooid.length..mm.","Number.of.zooids","Slope_O2", "Slope_O2.control", "Slope_O2_dif", "Timespan","Temperature_range","Treatment","Paired")
for(i in 1:length(unique(paste0(presens$Specimen,presens$Treatment)))){
  spm_i <- unique(paste(presens$Specimen,presens$Treatment),sep=" ")[i] %>% 
    str_split(pattern=" ") %>%
    .[[1]] %>% .[1]
  treat_i <- unique(paste(presens$Specimen,presens$Treatment),sep=" ")[i] %>% 
    str_split(pattern=" ") %>%
    .[[1]] %>% .[2]
  series_i <- presens[which(presens$Specimen==spm_i & presens$Treatment==treat_i),c("Specimen","Species","Zooid.length..mm.","Number.of.zooids","Colony.volume..ml.","Time.point..min.","Temperature...C.","abs_O2.mg._animal","abs_O2.mg._control","Treatment","is.paired")]
  lm_animal <- lm(abs_O2.mg._animal~Time.point..min., series_i)
  lm_control <- lm(abs_O2.mg._control~Time.point..min., series_i)
  print(series_i$Specimen[1] %>% as.character());print(series_i$Species[1] %>% as.character());print(lm_animal$coefficients);print(lm_control$coefficients)
  time.span_i <- max(series_i$Time.point..min.)-min(series_i$Time.point..min.)
  #"Species","Specimen","Colony.volume..ml.","Zooid.length..mm.","Number.of.zooids",
  slopes[i,1] <- series_i$Species[1] %>% as.character()
  slopes[i,2] <- series_i$Specimen[1] %>% as.character()
  slopes[i,3] <- series_i$Colony.volume..ml.[1] %>% as.numeric()
  slopes[i,4] <- series_i$Zooid.length..mm.[1] %>% as.numeric()
  slopes[i,5] <- series_i$Number.of.zooids[1] %>% as.numeric()
  #"Slope_O2", "Slope_O2.control", "Slope_O2_dif"
  slopes[i,6] <- lm_animal$coefficients[2]
  slopes[i,7] <- lm_control$coefficients[2]
  slopes[i,8] <- lm_animal$coefficients[2] - lm_control$coefficients[2]
  #"Timespan","Temperature_range","Treatment","Paired"
  slopes[i,9] <- time.span_i
  slopes[i,10] <- if(!is.na(mean(series_i$Temperature...C., na.rm = T))){max(series_i$Temperature...C., na.rm=T) - min(series_i$Temperature...C., na.rm=T)} else NA
  print(series_i$Treatment)
  slopes[i,11] <- series_i$Treatment %>% unique()
  slopes[i,12] <- series_i$is.paired %>% unique()
}

#Normalize by colony volume, transform to PicoGrams of O2
slopes %>% mutate(Slope_normalized = 1000000*Slope_O2/Colony.volume..ml., Slope_O2_dif_normalized = 1000000*Slope_O2_dif/Colony.volume..ml.) -> slopes

#Get carbon estimates
#Estimate Carbon content for each species and recalculate carbon-based rawCOT
mm_to_carbon <- read.csv("Madin1981_salpcarbon.tsv", header=T, stringsAsFactors = F, sep='\t')
mm_to_carbon$Generation[mm_to_carbon$Species=="Iasis (Weelia) cylindrica"] <- "a" ##PROXY 
mm_to_carbon$Species[mm_to_carbon$Species=="Thalia democratica"] <- "Thalia sp."  #### amalgamation!

#Normalize by carbon
slopes %>% left_join(mm_to_carbon[which(mm_to_carbon$Generation=="a"),c(1,4,5)], by="Species") %>%
  mutate(mgC = Regression_b*Zooid.length..mm.^Regression_a) %>% 
  mutate(Slopes_mgC = 1000000*Slope_O2/mgC, Slopes_dif_mgC = 1000000*Slope_O2_dif/mgC) -> slopes

slopes$Species %>% factor(levels=c("Pegea sp.", "Helicosalpa virgula","Cyclosalpa affinis", 
                                      "Cyclosalpa bakeri", "Cyclosalpa quadriluminis","Cyclosalpa polae",
                                      "Cyclosalpa sewelli", "Ritteriella sp.", "Ritteriella amboinensis",
                                      "Brooksia rostrata", "Thalia sp.", "Metcalfina hexagona", "Salpa fusiformis", "Salpa aspera", "Salpa maxima", 
                                      "Iasis (Weelia) cylindrica", "Soestia zonaria")) -> slopes$Species

pdf("Figures_respirometry/Kona2022/slopes_raw.pdf", height=6, width=10)
slopes %>% filter(Treatment=="Intact") %>% 
  ggplot(aes(x=Species,y=-Slope_O2))->slope_plot
slope_plot+
  geom_violin(color="red",fill="red", alpha=0.7)+
  geom_point(color="red",fill="red", alpha=0.7)+
  geom_violin(data=slopes %>% filter(Treatment=="Anesthetized"), color="blue", fill="blue", alpha=0.7, width = 0.5)+
  geom_point(data=slopes %>% filter(Treatment=="Anesthetized"), color="blue",fill="red", alpha=0.7)+
  ylab("Gross respiration rate (pgO2/min)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))#+
  #facet_wrap(~Paired)
dev.off()

pdf("Figures_respirometry/Kona2022/slope_diffs.pdf", height=6, width=10)
slopes %>% filter(Treatment=="Intact") %>% 
  ggplot(aes(x=Species,y=-Slope_O2_dif))+
  geom_violin(color="red",fill="red", alpha=0.7)+
  geom_point(color="red",fill="red", alpha=0.7)+
  geom_violin(data=slopes %>% filter(Treatment=="Anesthetized"), color="blue", fill="blue", alpha=0.7, width = 0.5)+
  geom_point(data=slopes %>% filter(Treatment=="Anesthetized"), color="blue",fill="red", alpha=0.7)+
  ylab("Net respiration rate (mgO2/min)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))#+
#facet_wrap(~Paired)
dev.off()

pdf("Figures_respirometry/Kona2022/slopes_volume.pdf", height=6, width=10)
slopes %>% filter(Treatment=="Intact") %>% 
  ggplot(aes(x=Species,y=-Slope_normalized))+
  geom_violin(color="red",fill="red", alpha=0.7)+
  geom_point(color="red",fill="red", alpha=0.7)+
  geom_violin(data=slopes %>% filter(Treatment=="Anesthetized"), color="blue", fill="blue", alpha=0.7, width = 0.5)+
  geom_point(data=slopes %>% filter(Treatment=="Anesthetized"), color="blue",fill="red", alpha=0.7)+
  ylab("Raw respiration rate (pgO2/min) per specimen biovolume (ml)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))
dev.off()

pdf("Figures_respirometry/Kona2022/slope_diffs_volume.pdf", height=6, width=10)
slopes %>% filter(Treatment=="Intact") %>% 
  ggplot(aes(x=Species,y=-Slope_O2_dif_normalized))+
  geom_violin(color="red",fill="red", alpha=0.7)+
  geom_point(color="red",fill="red", alpha=0.7)+
  geom_violin(data=slopes %>% filter(Treatment=="Anesthetized"), color="blue", fill="blue", alpha=0.7, width = 0.5)+
  geom_point(data=slopes %>% filter(Treatment=="Anesthetized"), color="blue",fill="red", alpha=0.7)+
  ylab("Net respiration rate (pgO2/min) per specimen biovolume (ml)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))
dev.off()

pdf("Figures_respirometry/Kona2022/slopes_mgC.pdf", height=6, width=10)
slopes %>% filter(Treatment=="Intact", !is.na(Slopes_mgC)) %>% 
  ggplot(aes(x=Species,y=-Slopes_mgC))+
  geom_violin(color="red",fill="red", alpha=0.7)+
  geom_point(color="red",fill="red", alpha=0.7)+
  geom_violin(data=slopes %>% filter(Treatment=="Anesthetized", !is.na(Slopes_mgC)), color="blue", fill="blue", alpha=0.7, width = 0.5)+
  geom_point(data=slopes %>% filter(Treatment=="Anesthetized", !is.na(Slopes_mgC)), color="blue",fill="red", alpha=0.7)+
  ylab("Raw respiration rate (pgO2/min) / Specimen carbon (mg)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))
dev.off()

pdf("Figures_respirometry/Kona2022/slopes_mgC_dif.pdf", height=6, width=10)
slopes %>% filter(Treatment=="Intact", !is.na(Slopes_dif_mgC)) %>% 
  ggplot(aes(x=Species,y=-Slopes_dif_mgC))+
  geom_point(color="red",fill="red", alpha=0.7)+
  geom_point(data=slopes %>% filter(Treatment=="Anesthetized", !is.na(Slopes_dif_mgC)), color="blue",fill="red", alpha=0.7)+
  ylab("Net respiration rate (pgO2/min) / Specimen carbon (mg)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))
dev.off()

## GROUP BY Specimen ##

### COST OF TRANSPORT ###
#Get swimming speeds
salplit <- read.csv("salplit.tsv", sep='\t', stringsAsFactors = F)
swimraw <- salplit[,c(1,2,4)] %>% filter(Variable=="Mean swimming speed cms")
swim <- data.frame(Species=swimraw$Species, Speed.cm.s=as.numeric(swimraw$Value))

swim <- read.csv("SalpPreliminaryPass.tsv",sep='\t', stringsAsFactors = F)[,c(1,15)]
swim <- swim[!is.na(swim$Speed..cm.s.),]
names(swim)[2] <- "Speed.cm.s"

      #### CAUTION ##### WE ARE SWAPPING NAMES HERE
swim$Species[swim$Species=="Pegea socia"] <- "Pegea sp."
swim$Species[swim$Species=="Ritteriella retracta"] <- "Ritteriella sp."
swim$Species[swim$Species=="Rittereilla amboinensis"] <- "Ritteriella amboinensis"
swim$Species[swim$Species=="Thalia democratica"] <- "Thalia sp."

energetics <- mutate(slopes, Species=as.character(Species))
energetics$Species[energetics$Species=="Iasis (Weelia) cylindrica"] <- "Iasis cylindrica"
energetics <- left_join(energetics, swim, by="Species") 
energetics$Zooid.length..mm. <- as.numeric(energetics$Zooid.length..mm.)
energetics$Number.of.zooids <- as.numeric(energetics$Number.of.zooids)
energetics %>% mutate(Speed.body.s = 10*Speed.cm.s/Zooid.length..mm.) -> energetics #speed per body lengths
energetics %>% group_by(Species) %>% summarize(-Slope_O2_dif_normalized %>% mean(na.rm=T))

#Define COST of LIVING: COT.abs as mgO2/(Zooid.vol*cm_moved)  v.v. COT.rel as mgO2/(Zooid.vol*bodylength)
energetics %>% filter(energetics$Treatment=="Intact") %>% 
  mutate(COT.abs.ml = -Slope_O2_dif_normalized/(Speed.cm.s*60), 
         COT.rel.ml = -Slope_O2_dif_normalized/(Speed.body.s*60), 
         COT.abs.mgC = -Slope_O2_dif_normalized/(Speed.cm.s*60), 
         COT.rel.mgC = -Slope_O2_dif_normalized/(Speed.body.s*60)) -> COT_intact
COT_intact <- COT_intact[which(!is.na(COT_intact$COT.abs.ml)),]

COT_intact$Species %>% factor(levels=c("Pegea sp.", "Helicosalpa virgula","Cyclosalpa affinis", 
                                        "Cyclosalpa bakeri", "Cyclosalpa quadriluminis","Cyclosalpa polae",
                                        "Cyclosalpa sewelli", "Ritteriella sp.", "Ritteriella amboinensis",
                                        "Brooksia rostrata", "Thalia sp.", "Metcalfina hexagona", "Salpa fusiformis", "Salpa aspera", "Salpa maxima", 
                                        "Iasis cylindrica", "Soestia zonaria")) -> COT_intact$Species

#raw cost of living by volume per cm
pdf("Figures_respirometry/Kona2022/rawCOT_intact_abs.pdf", height=6, width=10)
COT_intact %>% filter(!is.na(COT.abs.ml)) %>% 
  ggplot(aes(x=Species,y=COT.abs.ml))+
  geom_boxplot(aes(fill=Speed.cm.s))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90))+ylab("Gross Cost of Transport (pgO2/Biovolume per cm moved)") #+ ylim(c(0, 7e-05))
dev.off()

#raw cost of living by volume per body length
pdf("Figures_respirometry/Kona2022/rawCOT_intact_rel.pdf", height=6, width=10)
COT_intact %>% filter(!is.na(COT.rel.ml)) %>% 
  ggplot(aes(x=Species,y=COT.rel.ml))+
  geom_boxplot(aes(fill=Speed.cm.s))+
  #scale_fill_gradient(low="white",high="red")+
  theme_bw()+theme(axis.text.x = element_text(angle = 90))+ylab("Gross Cost of Transport (pgO2/Biovolume per zooid length moved)")  #+ ylim(c(0, 6e-05))
dev.off()

#raw cost of living by carbon per cm
pdf("Figures_respirometry/Kona2022/COT_mgC_abs.pdf", height=6, width=10)
COT_intact %>% filter(!is.na(COT.abs.mgC)) %>% 
  ggplot(aes(x=Species,y=COT.abs.mgC))+
  geom_boxplot(aes(fill=Speed.cm.s))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90))+ylab("Gross Cost of Transport (pgO2/mgC per cm moved)") #+ ylim(c(0, 6e-05))
dev.off()

#raw cost of living by carbon per bodylength
pdf("Figures_respirometry/Kona2022/COT_mgC_rel.pdf", height=6, width=10)
COT_intact %>% filter(!is.na(COT.rel.mgC)) %>% 
  ggplot(aes(x=Species,y=COT.rel.mgC))+
  geom_boxplot(aes(fill=Speed.cm.s))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90))+ylab("Gross Cost of Transport (pgO2/mgC per zooid length moved)")
dev.off()

#COT for paired species 
COT <- energetics[,c("Species","Specimen","Colony.volume..ml.","Zooid.length..mm.","Number.of.zooids","Treatment","Paired","Slope_O2_dif_normalized","Slopes_dif_mgC","Speed.cm.s","Speed.body.s")]
  
##### WARNING ##### SKETCHY STUFF TO REMOVE NEGATIVES !!!!
COT$Slope_O2_dif_normalized[which(COT$Slope_O2_dif_normalized>0)] <- 0

COT %>% filter(Paired=="Yes") %>% 
  pivot_wider(names_from = Treatment, values_from = c(Slope_O2_dif_normalized, Slopes_dif_mgC)) %>% 
  mutate(COT.abs.ml = -(Slope_O2_dif_normalized_Intact-Slope_O2_dif_normalized_Anesthetized)/(Speed.cm.s*60), 
         COT.rel.ml = -(Slope_O2_dif_normalized_Intact-Slope_O2_dif_normalized_Anesthetized)/(Speed.body.s*60)) %>% 
  mutate(COT.abs.mgC = -(Slopes_dif_mgC_Intact-Slopes_dif_mgC_Anesthetized)/(Speed.cm.s*60), 
         COT.rel.mgC = -(Slopes_dif_mgC_Intact-Slopes_dif_mgC_Anesthetized)/(Speed.body.s*60)) -> COT_paired

#COT by biovolume per cm across species
COT_paired %>% filter(!is.na(COT.abs.ml)) %>% 
ggplot(aes(x=Species,y=(COT.abs.ml) ))+
  geom_boxplot(aes(col=Speed.cm.s, fill= Speed.cm.s))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Net cost of Transport (pgO2/ml per cm moved)")

#COT by biovolume per zooid length across species
COT_paired %>% filter(!is.na(COT.rel.ml)) %>% 
  ggplot(aes(x=Species,y=(COT.rel.ml) ))+
  geom_boxplot(aes(col=Speed.cm.s, fill= Speed.cm.s))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Net cost of Transport (pgO2/ml per zooid length moved)")

#COT by mgC per cm across species
COT_paired %>% filter(!is.na(COT.abs.mgC)) %>%
  ggplot(aes(x=Species,y=(COT.abs.mgC) ))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Net cost of Transport (pgO2/mgC per cm moved)")

#COT by mgC per zooid length across species
COT_paired %>% filter(!is.na(COT.rel.mgC)) %>% 
  ggplot(aes(x=Species,y=COT.rel.mgC ))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Net cost of Transport (pgO2/mgC per zooid length moved)")

###  COT paired by species   ###
COT %>% filter(!is.na(Species)) %>% 
  pivot_wider(names_from = Treatment, values_from = c(Slope_O2_dif_normalized, Slopes_dif_mgC)) %>% 
  group_by(Species) %>% 
  summarise_at(vars("Speed.cm.s", "Speed.body.s", "Slope_O2_dif_normalized_Intact", "Slope_O2_dif_normalized_Anesthetized", "Slopes_dif_mgC_Intact", "Slopes_dif_mgC_Anesthetized", "Zooid.length..mm."),
               function(x){median(x, na.rm=T)}) %>% as.data.frame() -> COT_species

##### WARNING ##### SKETCHY STUFF TO REMOVE NEGATIVES !!!!
#COT_species$Slope_O2_dif_normalized_Intact[which(COT_species$Slope_O2_dif_normalized_Intact>COT_species$Slope_O2_dif_normalized_Anesthetized)] <- COT_species$Slope_O2_dif_normalized_Anesthetized[which(COT_species$Slope_O2_dif_normalized_Intact>COT_species$Slope_O2_dif_normalized_Anesthetized)]

#Systemic corrections
COT_species[,c(1,4,5)] -> playDF
  ## Calculate maximum negative deviation between intact and anesthetized across species
max(playDF[which(playDF$Slope_O2_dif_normalized_Intact>
               playDF$Slope_O2_dif_normalized_Anesthetized),"Slope_O2_dif_normalized_Intact"]-
  playDF[which(playDF$Slope_O2_dif_normalized_Intact>
                 playDF$Slope_O2_dif_normalized_Anesthetized),"Slope_O2_dif_normalized_Anesthetized"]) -> Intact_systemic_correction
playDF$Slope_O2_dif_normalized_Intact <- playDF$Slope_O2_dif_normalized_Intact - Intact_systemic_correction
  ## Calculate maximum anesthetized slope
# max(playDF$Slope_O2_dif_normalized_Anesthetized, na.rm = T) -> Anesthetized_systemic_correction
# playDF$Slope_O2_dif_normalized_Anesthetized <- playDF$Slope_O2_dif_normalized_Anesthetized - Anesthetized_systemic_correction
# playDF$Slope_O2_dif_normalized_Intact <- playDF$Slope_O2_dif_normalized_Intact - Anesthetized_systemic_correction

playDF %>%  mutate(Percent=100*round((Slope_O2_dif_normalized_Intact-Slope_O2_dif_normalized_Anesthetized)/Slope_O2_dif_normalized_Intact,6))-> playDF
playDF
 COT_species[,c("Slope_O2_dif_normalized_Intact","Slope_O2_dif_normalized_Anesthetized")] <- playDF[,c("Slope_O2_dif_normalized_Intact","Slope_O2_dif_normalized_Anesthetized")]

 COT_species %>% mutate(COT.abs.ml = -(Slope_O2_dif_normalized_Intact-Slope_O2_dif_normalized_Anesthetized)/(Speed.cm.s*60), 
                        COT.rel.ml = -(Slope_O2_dif_normalized_Intact-Slope_O2_dif_normalized_Anesthetized)/(Speed.body.s*60)) %>% 
   mutate(COT.abs.mgC = -(Slopes_dif_mgC_Intact-Slopes_dif_mgC_Anesthetized)/(Speed.cm.s*60), 
          COT.rel.mgC = -(Slopes_dif_mgC_Intact-Slopes_dif_mgC_Anesthetized)/(Speed.body.s*60)) %>% 
   mutate(COT.p.ml = -100*(Slope_O2_dif_normalized_Anesthetized-Slope_O2_dif_normalized_Intact)/Slope_O2_dif_normalized_Intact, 
          COT.p.mgC = -100*(Slopes_dif_mgC_Anesthetized-Slopes_dif_mgC_Intact)/Slopes_dif_mgC_Intact) %>% 
   mutate(COL.abs.ml = -Slope_O2_dif_normalized_Intact/(Speed.cm.s*60), 
          COL.rel.ml = -Slope_O2_dif_normalized_Intact/(Speed.body.s*60), 
          COL.abs.mgC = -Slope_O2_dif_normalized_Intact/(Speed.cm.s*60), 
          COL.rel.mgC = -Slope_O2_dif_normalized_Intact/(Speed.body.s*60))-> COT_species 
 
#order species as factors by architecture
COT_species$Species %>% factor(levels=c("Pegea sp.", "Helicosalpa virgula","Cyclosalpa affinis", 
                                        "Cyclosalpa bakeri", "Cyclosalpa quadriluminis","Cyclosalpa polae",
                                        "Cyclosalpa sewelli", "Ritteriella sp.", "Ritteriella amboinensis",
                                        "Brooksia rostrata", "Thalia sp.", "Metcalfina hexagona", "Salpa fusiformis", "Salpa aspera", "Salpa maxima", 
                                        "Iasis cylindrica", "Soestia zonaria")) -> COT_species$Species

#COT_species <- filter(COT_species, Species != "Cyclosalpa bakeri")  ### REMOVED OUTLIER WONKY SPP !!!
COT_species[,c(1,15,9,13)]


#Net respiration rate by biovolume in Swimmers and KO
COT_species %>% filter(!is.na(COT.abs.ml)) %>% 
  ggplot(aes(x=Species))+
  geom_point(aes(y=(-Slope_O2_dif_normalized_Intact)),color="red", alpha=0.7)+
  geom_point(aes(y=(-Slope_O2_dif_normalized_Anesthetized)), color="blue",alpha=0.7)+
  ylab("Net respiration rate (pgO2/min) per specimen biovolume (ml)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))

COT_species %>% filter(!is.na(COT.abs.ml)) %>% 
  ggplot(aes(x=Species))+
  geom_bar(aes(y=5-Slope_O2_dif_normalized_Intact+Slope_O2_dif_normalized_Anesthetized), stat="identity")+
  ylab("Net respiration rate (pgO2/min) per specimen biovolume (ml)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))

#COT by biovolume per cm across absolute speeds
ggplot(COT_species %>% filter(!is.na(COT.abs.ml)), aes(x=Speed.cm.s,y=COT.abs.ml))+
  #geom_point(aes(color=Speed.body.s ),cex=3)+
  theme_bw()+
  geom_text(label=COT_species %>% filter(!is.na(COT.abs.ml)) %>% .$Species, hjust=0.4, vjust=-0.7)+
  ylab("Net cost of Transport (pgO2/ml per cm moved)")

#COT by biovolume per cm across relative speeds
ggplot(COT_species %>% filter(!is.na(COT.abs.ml)), aes(x=Speed.body.s,y=COT.abs.ml))+
  geom_point(aes(color=Speed.cm.s ),cex=3)+
  theme_bw()+
  geom_text(label=COT_species %>% filter(!is.na(COT.abs.ml)) %>% .$Species, hjust=0.4, vjust=-0.7)+
  ylab("Net cost of Transport (pgO2/ml per body length moved)")

# COT_species <- mutate(COT_species, percentSwim = -5000*(Slope_O2_dif_normalized_Intact - Slope_O2_dif_normalized_Anesthetized) )
# 
# ggplot(COT_species %>% filter(!is.na(percentSwim)), aes(x=Speed.cm.s,y=percentSwim))+
#   geom_point()+
#   geom_text(label=COT_species %>% filter(!is.na(percentSwim)) %>% .$Species)+
#   theme_bw()

#COT by biovolume per cm across species
ggplot(COT_species %>% filter(!is.na(COT.abs.ml)), aes(x=Species,y=COT.abs.ml))+
  geom_point(aes(color=Speed.cm.s))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Net cost of Transport (pgO2/mgC per cm moved)")

#COT by biovolume per zooid length across species
ggplot(COT_species %>% filter(!is.na(COT.rel.ml)), aes(x=Species,y=COT.rel.ml))+
  geom_point(aes(color=Speed.body.s))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Net cost of Transport (pgO2/ml per body length moved)")

#COT by biovolume per zooid length across absolute speeds
ggplot(COT_species %>% filter(!is.na(COT.rel.ml)), aes(x=Speed.cm.s,y=COT.rel.ml))+
  geom_point(aes(color=Speed.cm.s))+
  theme_bw()+
  geom_text(label=COT_species %>% filter(!is.na(COT.rel.ml)) %>% .$Species, hjust=0.4, vjust=-0.7)+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Net cost of Transport (pgO2/ml per body length moved)")

#COT by biovolume per zooid length across relative speeds
ggplot(COT_species %>% filter(!is.na(COT.rel.ml)), aes(x=Speed.body.s,y=COT.rel.ml))+
  geom_point(aes(color=Speed.body.s))+
  theme_bw()+
  geom_text(label=COT_species %>% filter(!is.na(COT.rel.ml)) %>% .$Species, hjust=0.4, vjust=-0.7)+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Net cost of Transport (pgO2/ml per body length moved)")

#COT by mgC per zooid length across species
ggplot(COT_species %>% filter(!is.na(COT.rel.mgC)), aes(x=Species,y=COT.rel.mgC))+
  geom_point(aes(size=Zooid.length..mm.,color=Speed.body.s))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Cost of Transport (pgO2/mgC per body length moved)") 

# % cost (by biovolume) invested in swimming across species order: COT.p, color: Speed_cm
ggplot(COT_species %>%  filter(!is.na(COT.p.ml)), aes(x = reorder(Species %>% as.character(), Speed.cm.s),y=COT.p.ml))+
  geom_point(aes(size=Zooid.length..mm.,color=log(Speed.cm.s)))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
 # ylim(c(0,100))+
  xlab("Species")+
  ylab("% Cost of Transport Volume-normalized")

# % cost (by biovolume) invested in swimming across species order: COT.p, , color: COT
ggplot(COT_species %>%  filter(!is.na(COT.p.ml)), aes(x = reorder(Species %>% as.character(), COL.abs.ml),y=COT.p.ml))+
  geom_point(aes(size=Zooid.length..mm.,color=COL.abs.ml))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  # ylim(c(0,100))+
  xlab("Species")+
  ylab("% Cost of Transport Volume-normalized")


# % cost (by biovolume) invested in swimming across absolute SPEED , color: Species
ggplot(COT_species %>%  filter(!is.na(COT.p.ml)), aes(x = Speed.cm.s,y=COT.p.ml))+
  geom_point(aes(size=Zooid.length..mm.,color=Species))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  #ylim(c(0,100))+
  ylab("% Cost of Transport Volume-normalized")

# % cost (by biovolume) invested in swimming across relative SPEED , color: Species
ggplot(COT_species %>%  filter(!is.na(COT.p.ml)), aes(x = Speed.body.s,y=COT.p.ml))+
  geom_point(aes(size=Zooid.length..mm.,color=Species))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  #ylim(c(0,100))+
  ylab("% Cost of Transport Volume-normalized")

# % cost (by biovolume) invested in swimming across COT_vol , color: Species
ggplot(COT_species %>%  filter(!is.na(COT.p.ml)), aes(x = COT.abs.ml,y=COT.p.ml))+
  geom_point(aes(size=Zooid.length..mm.,color=Species))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  #ylim(c(0,100))+
  ylab("% Cost of Transport Volume-normalized")

# % cost (by biovolume) invested in swimming across basal rate , color: Species
ggplot(COT_species %>%  filter(!is.na(COT.p.ml)), aes(x = COL.abs.ml-COT.abs.ml,y=COT.p.ml))+
  geom_point(aes(size=Zooid.length..mm.,color=Species))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  #ylim(c(0,100))+
  ylab("% Cost of Transport Volume-normalized")+xlab("Basal rate")

# % cost (by biovolume) invested in swimming across species , color: Relative speed
ggplot(COT_species %>% filter(!is.na(COT.p.mgC)), aes(x=Species,y=COT.p.mgC))+
  geom_point(aes(size=Speed.cm.s,color=Speed.body.s))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  #ylim(c(0,100))+
  ylab("% Cost of Transport Carbon-normalized") #

### plots with speed  ###

ggplot(COT_intact, aes(x=Speed.cm.s,y=COT.abs.ml))+
  geom_point()+
  theme_bw()+
 geom_smooth(method="lm")+
  ylab("Cost of Living (mgO2/ml per body length moved)")

ggplot(COT_species, aes(x=Speed.cm.s,y=COT.abs.ml))+
  geom_point()+
  geom_text(aes(label=Species), vjust=-0.8)+
  theme_bw()+
  geom_smooth(method="lm")+
  ylab("Cost of Transport (mgO2/mgC per body length moved)")

