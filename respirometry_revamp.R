############ LOAD UP AND SET UP ##############

# COREECT exponentially for body mass/vol ???????? ######

library(tidyverse)
library(wql)
library(ggrepel)
library(patchwork)
require(data.table)
library(mgcv)
setwd("~/salp_ecomorphology/")

#Load data and label: control rows, paired specimens
presens <- read.csv("~/salp_ecomorphology/respirometry_Kona21-23.tsv", header = T, sep = "\t", stringsAsFactors = F)
presens$Species[which(is.na(presens$Species))] <- "Control"
presens$Specimen[which(is.na(presens$Specimen))] <- "Control"
presens <- mutate(presens, is.paired=ifelse(grepl("Paired", Measurement.notes), "Yes", "No"))
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

# Make SM Table 2

summary_table <- presens %>%
  group_by(Specimen) %>%
  summarize(
    Experiment = first(Experiment),
    Date = first(Date),
    Species = first(Species),
    Activity.level = first(Activity.level),
    Zooid.length..mm. = first(Zooid.length..mm.),
    Number.of.zooids = first(Number.of.zooids),
    Colony.volume..ml. = first(Colony.volume..ml.),
    Container.volume..ml. = first(Container.volume..ml.),
    Organism.notes = first(Organism.notes),
    is.paired = first(is.paired),
    Number.of.measurements = n()
  )

# Print the summary table
summary_table[which(summary_table$Species == "Iasis (Weelia) cylindrica"),"Species"]<-"Iasis cylindrica"
summary_table[which(summary_table$Species == "Ritteriella retracta"),"Species"]<-"Ritteriella sp."
summary_table[which(summary_table$Species == "Thalia cicar"),"Species"]<-"Thalia sp."
print(summary_table)

#Species level summary table
summary_table %>% group_by(Species) %>% summarise(
  `Mean Number of zooids respirometry` = mean(Number.of.zooids, na.rm = T),
  `Mean zooid length (mm) respirometry` = mean(Zooid.length..mm., na.rm = T),
  `Mean Colony volume (ml)` = mean(`Colony.volume..ml.`, na.rm = T),
  `Number of Respirometry Specimens` = n(),
  `Number of Respirometry Measurements` = sum(Number.of.measurements, na.rm = T)
) -> resp_spp

full_join(species_table, resp_spp, by=c("Species")) -> whole_summary
whole_summary[18,2]<-"Linear"
write_csv(whole_summary, "SMTable3_Species.csv")

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

pdf("Figures_respirometry/Kona2021-2023/RawO2.pdf", height=4, width=14)
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

slopes %>% filter(Slope_O2_dif<0) -> slopes

#Normalize by colony volume, transform to PicoGrams of O2
    # slopes %>% mutate(Slope_normalized = 1000000*Slope_O2/(Colony.volume..ml.^0.75), Slope_O2_dif_normalized = 1000000*Slope_O2_dif/(Colony.volume..ml.^0.75)) -> slopes
slopes %>% mutate(Slope_normalized = 1000000*Slope_O2/(Colony.volume..ml.), Slope_O2_dif_normalized = 1000000*Slope_O2_dif/(Colony.volume..ml.)) -> slopes

#Get carbon estimates
#Estimate Carbon content for each species and recalculate carbon-based rawCOT
mm_to_carbon <- read.csv("~/salp_ecomorphology/Madin1981_salpcarbon.tsv", header=T, stringsAsFactors = F, sep='\t')
mm_to_carbon$Generation[mm_to_carbon$Species=="Iasis (Weelia) cylindrica"] <- "a" ##PROXY 
mm_to_carbon$Species[mm_to_carbon$Species=="Thalia democratica"] <- "Thalia sp."  #### amalgamation!

#Normalize by carbon
slopes %>% left_join(mm_to_carbon[which(mm_to_carbon$Generation=="a"),c(1,4,5)], by="Species") %>%
  mutate(mgC = Regression_b*Zooid.length..mm.^Regression_a) %>% 
  mutate(Slopes_mgC = 1000000*Slope_O2/mgC, Slopes_dif_mgC = 1000000*Slope_O2_dif/mgC) -> slopes

slopes$Species[which(slopes$Species == "Thalia cicar")] <- "Thalia sp."
slopes$Species[which(slopes$Species == "Ritteriella sp.")] <- "Ritteriella retracta"
slopes$Species %>% factor(levels=c("Pegea sp.", "Helicosalpa virgula","Cyclosalpa affinis", 
                                      "Cyclosalpa bakeri", "Cyclosalpa quadriluminis","Cyclosalpa polae",
                                      "Cyclosalpa sewelli", "Ritteriella retracta", "Ritteriella amboinensis",
                                      "Brooksia rostrata", "Thalia sp.", "Metcalfina hexagona", "Salpa fusiformis", "Salpa aspera", "Salpa maxima", 
                                      "Iasis (Weelia) cylindrica", "Ihlea punctata", "Soestia zonaria")) -> slopes$Species

# Specify the columns from slopes that you want
selected_columns <- c("Specimen", "Timespan", "Temperature_range", "Treatment", 
                      "Slope_O2", "Slope_O2.control", "Slope_O2_dif", 
                      "Slope_normalized", "Slope_O2_dif_normalized")

# Left join the summary_table with selected columns from slopes based on Specimen
SM2_table <- left_join(summary_table, slopes %>% select(all_of(selected_columns)), by = "Specimen") %>% as.data.frame()
SM2_table <- SM2_table[,c(1,4,2,3,5:ncol(SM2_table))]

write.table(SM2_table, "SMTable2.tsv", col.names = T, row.names = F, sep = "\t")
                        
# VIZ

pdf("Figures_respirometry/Kona2021-2023/slopes_raw.pdf", height=6, width=10)
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

pdf("Figures_respirometry/Kona2021-2023/slope_diffs.pdf", height=6, width=10)
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

pdf("Figures_respirometry/Kona2021-2023/slopes_volume.pdf", height=6, width=10)
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

pdf("Figures_respirometry/Kona2021-2023/slope_diffs_volume.pdf", height=6, width=10)
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

pdf("Figures_respirometry/Kona2021-2023/slopes_mgC.pdf", height=6, width=10)
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

pdf("Figures_respirometry/Kona2021-2023/slopes_mgC_dif.pdf", height=6, width=10)
slopes %>% filter(Treatment=="Intact", !is.na(Slopes_dif_mgC)) %>% 
  ggplot(aes(x=Species,y=-Slopes_dif_mgC))+
  geom_point(color="red",fill="red", alpha=0.7)+
  geom_point(data=slopes %>% filter(Treatment=="Anesthetized", !is.na(Slopes_dif_mgC)), color="blue",fill="red", alpha=0.7)+
  ylab("Net respiration rate (pgO2/min) / Specimen carbon (mg)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))
dev.off()

## SM Figure S5 ##

slopes %>% filter(-Slope_O2_dif_normalized/1000000 < 0.0025) %>% 
  ggplot(aes(x=Species,y=-Slope_O2_dif_normalized/1000000))+
  geom_boxplot(aes(fill=Treatment), size = 0.5)+
  ylab("Respiration rate (mgO2/min) per biovolume (ml)")+
  scale_fill_manual(values = c("lightblue","pink"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))

slopes %>%  ggplot(aes(x=Species,y=-1.025*60*31.25*Slope_O2_dif_normalized/1000000))+
  geom_boxplot(aes(fill=Treatment))+
  ylab("Respiration rate (µmolO2/g/h)")+
  ylim(NA,3)+
  scale_fill_manual(values = c("lightblue","pink"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))


## GROUP BY Specimen ##


### COST OF TRANSPORT ###
#Get swimming speeds from literature
#salplit <- read.csv("~/salp_ecomorphology/salplit.tsv", sep='\t', stringsAsFactors = F)
#swimraw <- salplit[,c(1,2,4)] %>% filter(Variable=="Mean swimming speed cms")
#swim <- data.frame(Species=swimraw$Species, Speed.cm.s=as.numeric(swimraw$Value))

#Get swimming speeds form preliminary pass in 2D
#swim <- read.csv("SalpPreliminaryPass.tsv",sep='\t', stringsAsFactors = F)[,c(1,15)]
#swim <- swim[!is.na(swim$Speed..cm.s.),]
#names(swim)[2] <- "Speed.cm.s"

#Get swimming speeds on Event Measure
swim <- read.csv("~/salp_ecomorphology/EMSpeeds_final_annotated.tsv",stringsAsFactors = F) #[,c(20,21,25)]
summarized_raw=swim
# rawswim <- read.csv("EMspeed_R_annotated.tsv",stringsAsFactors = F)
# summarized_raw <- rawswim[,c(1,20,21,22,24)] %>% group_by(Filename) %>%
#   summarize(
#     Mean_Speed_cms = mean(Speed_mms_abs)/1000,
#     Species = Species,
#     Architecture = Architecture,
#     Number.of.zooids = Zooid.number) %>% as.data.frame() %>% filter(!is.na(Species)) %>% unique()

ggplot(summarized_raw, aes(x=Species, y=Speed_mm_s))+
  geom_boxplot(aes(fill=Architecture))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90))

energetics <- mutate(slopes, Species=as.character(Species))
energetics$Species[energetics$Species=="Iasis (Weelia) cylindrica"] <- "Iasis cylindrica"
summarized_raw$Species[summarized_raw$Species=="Ritteriella sp."] <- "Ritteriella retracta"

#options: make both data sets species wise / 
  # keep specimen wise energetics and species wise speeds /
     # keep species wise energetics and specimen wise speeds /
        # keep specimen wise energetics and match by zooid number -- 
              ## build glm models of speed by zooid size and number for each species, then apply to the energetics 

summarized_raw[,c(3,5,4,9,10)]-> swimming_data
swimming_data %>% group_by(Species) %>% summarise(across(where(is.numeric), mean, na.rm = TRUE)) -> swimming_species

energetics_merged <- left_join(energetics, swimming_species, by="Species")

energetics_merged %>% mutate(Speed_mm_s = BLperSecond*Zooid.length..mm.) -> energetics_merged
energetics_merged %>% mutate(Speed.bl.p = BLperSecond/Pulses_per_second) -> energetics_merged

#energetics_merged %>% group_by(Species) %>% summarize(-Slope_O2_dif_normalized %>% mean(na.rm=T), Speed_mm_s %>% mean(na.rm=T)) %>% as.data.frame()


      #### CAUTION ##### WE ARE SWAPPING NAMES HERE
# swim$Species[swim$Species=="Pegea socia"] <- "Pegea sp."
# swim$Species[swim$Species=="Rittereilla amboinensis"] <- "Ritteriella amboinensis"
# swim$Species[swim$Species=="Thalia democratica"] <- "Thalia sp."

# energetics <- mutate(slopes, Species=as.character(Species))
# energetics$Species[energetics$Species=="Iasis (Weelia) cylindrica"] <- "Iasis cylindrica"
# energetics <- left_join(energetics, swim, by="Species") 
# energetics$Zooid.length..mm. <- as.numeric(energetics$Zooid.length..mm.)
# energetics$Number.of.zooids <- as.numeric(energetics$Number.of.zooids)
# energetics %>% mutate(Speed.body.s = 10*Speed.cm.s/Zooid.length..mm.) -> energetics #speed per body lengths
# energetics %>% group_by(Species) %>% summarize(-Slope_O2_dif_normalized %>% mean(na.rm=T))

#Define COST of LIVING: COT.abs as mgO2/(Zooid.vol*cm_moved)  v.v. COT.rel as mgO2/(Zooid.vol*bodylength)
# energetics_merged %>% filter(energetics_merged$Treatment=="Intact") %>% 
#   mutate(COT.abs.ml = -Slope_O2_dif_normalized/Speed.mm.s, 
#          COT.rel.ml = -Slope_O2_dif_normalized/BLperSecond, 
#          COT.blp.ml = -Slope_O2_dif_normalized/Speed.bl.p,
#          COT.abs.mgC = -Slopes_dif_mgC/Speed.mm.s, 
#          COT.rel.mgC = -Slopes_dif_mgC/BLperSecond) -> COT_intact
# COT_intact <- COT_intact[which(!is.na(COT_intact$COT.abs.ml)),]

# COT_intact$Species %>% factor(levels=c("Pegea sp.", "Helicosalpa virgula","Cyclosalpa affinis", 
#                                        "Cyclosalpa bakeri", "Cyclosalpa quadriluminis","Cyclosalpa polae",
#                                        "Cyclosalpa sewelli", "Ritteriella retracta", "Ritteriella amboinensis",
#                                        "Brooksia rostrata", "Thalia sp.", "Metcalfina hexagona", "Salpa fusiformis", "Salpa aspera", "Salpa maxima", 
#                                        "Iasis (Weelia) cylindrica", "Ihlea punctata", "Soestia zonaria")) -> COT_intact$Species

# COT_intact$Species %>% factor(levels=c("Pegea sp.","Cyclosalpa affinis", 
#                                        "Cyclosalpa bakeri", "Cyclosalpa polae",
#                                        "Cyclosalpa sewelli", "Ritteriella retracta", "Ritteriella amboinensis",
#                                        "Brooksia rostrata", "Salpa aspera", "Salpa maxima", 
#                                        "Iasis cylindrica")) -> COT_intact$Species
# 
# #raw cost of living by volume per cm
# pdf("Figures_respirometry/Kona2021-2023/rawCOT_intact_abs.pdf", height=6, width=10)
# COT_intact %>% filter(!is.na(COT.abs.ml) & !is.na(Species)) %>% 
#   ggplot(aes(x=Species,y=COT.abs.ml %>% log()))+
#   geom_boxplot(aes(fill=Speed.mm.s))+
#   theme_bw()+theme(axis.text.x = element_text(angle = 90))+ylab("Gross Cost of Transport (pgO2/Biovolume per mm moved)") #+ ylim(c(0, 7e-05))
# dev.off()
# 
# #raw cost of living by volume per body length
# pdf("Figures_respirometry/Kona2021-2023/rawCOT_intact_rel.pdf", height=6, width=10)
# COT_intact %>% filter(!is.na(COT.rel.ml) & !is.na(Species)) %>% 
#   ggplot(aes(x=Species,y=COT.rel.ml %>% log()))+
#   geom_boxplot(aes(fill=Speed.mm.s))+
#   #scale_fill_gradient(low="white",high="red")+
#   theme_bw()+theme(axis.text.x = element_text(angle = 90))+ylab("Gross Cost of Transport (pgO2/Biovolume per zooid length moved)")  #+ ylim(c(0, 6e-05))
# dev.off()
# 
# 
# #raw cost of living by volume per body length per pulse
# 
# COT_intact %>% filter(!is.na(COT.rel.ml) & !is.na(Species)) %>% 
#   ggplot(aes(x=Species,y=COT.blp.ml %>% log()))+
#   geom_boxplot()+
#   #scale_fill_gradient(low="white",high="red")+
#   theme_bw()+theme(axis.text.x = element_text(angle = 90))+ylab("Gross Cost of Transport (pgO2/Biovolume per zooid length moved)")
# 
# #raw cost of living by carbon per cm
# pdf("Figures_respirometry/Kona2021-2023/COT_mgC_abs.pdf", height=6, width=10)
# COT_intact %>% filter(!is.na(COT.abs.mgC) & !is.na(Species)) %>% 
#   ggplot(aes(x=Species,y=COT.abs.mgC %>% log()))+
#   geom_boxplot(aes(fill=Speed.mm.s))+
#   theme_bw()+theme(axis.text.x = element_text(angle = 90))+ylab("Gross Cost of Transport (pgO2/mgC per mm moved)") #+ ylim(c(0, 6e-05))
# dev.off()
# 
# #raw cost of living by carbon per bodylength
# pdf("Figures_respirometry/Kona2021-2023/COT_mgC_rel.pdf", height=6, width=10)
# COT_intact %>% filter(!is.na(COT.rel.mgC) & !is.na(Species)) %>% 
#   ggplot(aes(x=Species,y=COT.rel.mgC %>% log()))+
#   geom_boxplot(aes(fill=Speed.mm.s))+
#   theme_bw()+theme(axis.text.x = element_text(angle = 90))+ylab("Gross Cost of Transport (pgO2/mgC per zooid length moved)")
# dev.off()

#COT for paired species 
COT <- energetics_merged[,c("Species","Specimen","Colony.volume..ml.","Zooid.length..mm.","Number.of.zooids",
                            "Treatment","Paired","Slope_O2_dif_normalized","Slopes_dif_mgC","Speed_mm_s",
                            "BLperSecond","Pulses_per_second", "Speed.bl.p")]

COT %>% filter(Slope_O2_dif_normalized>-3000) -> COT
  
##### WARNING ##### SKETCHY STUFF TO REMOVE NEGATIVES !!!!
#COT$Slope_O2_dif_normalized[which(COT$Slope_O2_dif_normalized>0)] <- 0

COT %>% filter(Paired=="Yes") %>% 
  pivot_wider(names_from = Treatment, values_from = c(Slope_O2_dif_normalized, Slopes_dif_mgC)) %>% 
  mutate(COT.abs.ml = -(Slope_O2_dif_normalized_Intact-Slope_O2_dif_normalized_Anesthetized)/Speed_mm_s, 
         COT.rel.ml = -(Slope_O2_dif_normalized_Intact-Slope_O2_dif_normalized_Anesthetized)/BLperSecond) %>% 
  mutate(COT.abs.mgC = -(Slopes_dif_mgC_Intact-Slopes_dif_mgC_Anesthetized)/Speed_mm_s, 
         COT.rel.mgC = -(Slopes_dif_mgC_Intact-Slopes_dif_mgC_Anesthetized)/BLperSecond) -> COT_paired

COT_paired %>% filter(Slope_O2_dif_normalized_Anesthetized>Slope_O2_dif_normalized_Intact) -> COT_paired

#COT by biovolume per cm across species
COT_paired %>% filter(!is.na(COT.abs.ml)) %>% 
ggplot(aes(x=Species,y=(COT.abs.ml) ))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Net cost of Transport (pgO2/ml per mm moved)")

COT_paired %>% filter(!is.na(COT.abs.ml) & Species != "Thalia sp.") %>% 
  ggplot(aes(x=Species,y=(COT.abs.ml) ))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Net cost of Transport (pgO2/ml per mm moved)")

#COT by biovolume per zooid length across species
COT_paired %>% filter(!is.na(COT.rel.ml)) %>% 
  ggplot(aes(x=Species,y=(COT.rel.ml) ))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Net cost of Transport (pgO2/ml per zooid length moved)")

COT_paired %>% filter(!is.na(COT.abs.ml) & Species != "Thalia sp.") %>% 
  ggplot(aes(x=Species,y=(COT.rel.ml) ))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Net cost of Transport (pgO2/ml per zooid length moved)")

# #COT by mgC per cm across species
# COT_paired %>% filter(!is.na(COT.abs.mgC)) %>%
#   ggplot(aes(x=Species,y=(COT.abs.mgC) ))+
#   geom_boxplot()+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90))+
#   ylab("Net cost of Transport (pgO2/mgC per mm moved)")
# 
# #COT by mgC per zooid length across species
# COT_paired %>% filter(!is.na(COT.rel.mgC)) %>% 
#   ggplot(aes(x=Species,y=COT.rel.mgC ))+
#   geom_boxplot()+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90))+
#   ylab("Net cost of Transport (pgO2/mgC per zooid length moved)")

###  COT paired by species   ###
# COT %>% filter(!is.na(Species)) %>% 
#   pivot_wider(names_from = Treatment, values_from = c(Slope_O2_dif_normalized, Slopes_dif_mgC)) %>% as.data.frame() -> COTpivot

COT %>% filter(!is.na(Species)) %>% 
  pivot_wider(names_from = Treatment, values_from = c(Slope_O2_dif_normalized, Slopes_dif_mgC)) %>% 
  group_by(Species) %>% 
  summarise_at(vars("Speed_mm_s", "BLperSecond", "Slope_O2_dif_normalized_Intact", 
                    "Slope_O2_dif_normalized_Anesthetized", "Slopes_dif_mgC_Intact", 
                    "Slopes_dif_mgC_Anesthetized", "Zooid.length..mm.","Pulses_per_second"),
               function(x){mean(x, na.rm=T)}) %>% as.data.frame() -> COT_species

COT_species %>% filter(Slope_O2_dif_normalized_Anesthetized>Slope_O2_dif_normalized_Intact) -> COT_species

# Merge the mean values with the original data
COT_with_means <- COT %>%
  left_join(COT_species %>% select(Species, Slope_O2_dif_normalized_Anesthetized, Slopes_dif_mgC_Anesthetized), by = "Species")

COT_with_means %>% filter(Slope_O2_dif_normalized_Anesthetized>Slope_O2_dif_normalized) -> COT_with_means
COT_species %>% filter(Slope_O2_dif_normalized_Anesthetized>Slope_O2_dif_normalized_Intact) -> COT_species

##### WARNING ##### SKETCHY STUFF TO REMOVE NEGATIVES !!!!
#COT_species$Slope_O2_dif_normalized_Intact[which(COT_species$Slope_O2_dif_normalized_Intact>COT_species$Slope_O2_dif_normalized_Anesthetized)] <- COT_species$Slope_O2_dif_normalized_Anesthetized[which(COT_species$Slope_O2_dif_normalized_Intact>COT_species$Slope_O2_dif_normalized_Anesthetized)]

#Systemic corrections
## Calculate maximum anesthetized slope
## max(playDF$Slope_O2_dif_normalized_Anesthetized, na.rm = T) -> Anesthetized_systemic_correction
## playDF$Slope_O2_dif_normalized_Anesthetized <- playDF$Slope_O2_dif_normalized_Anesthetized - Anesthetized_systemic_correction
## playDF$Slope_O2_dif_normalized_Intact <- playDF$Slope_O2_dif_normalized_Intact - Anesthetized_systemic_correction

# COT_species[,c(1,4,5)] -> playDF
#   ## Calculate maximum negative deviation between intact and anesthetized across species
# max(playDF[which(playDF$Slope_O2_dif_normalized_Intact>
#                playDF$Slope_O2_dif_normalized_Anesthetized),"Slope_O2_dif_normalized_Intact"]-
#   playDF[which(playDF$Slope_O2_dif_normalized_Intact>
#                  playDF$Slope_O2_dif_normalized_Anesthetized),"Slope_O2_dif_normalized_Anesthetized"]) -> Intact_systemic_correction
# playDF$Slope_O2_dif_normalized_Intact <- playDF$Slope_O2_dif_normalized_Intact - Intact_systemic_correction
# playDF %>%  mutate(Percent=100*round((Slope_O2_dif_normalized_Intact-Slope_O2_dif_normalized_Anesthetized)/Slope_O2_dif_normalized_Intact,6))-> playDF
# playDF
# COT_species[,c("Slope_O2_dif_normalized_Intact","Slope_O2_dif_normalized_Anesthetized")] <- playDF[,c("Slope_O2_dif_normalized_Intact","Slope_O2_dif_normalized_Anesthetized")]

 COT_species %>% mutate(COT.abs.ml = -(Slope_O2_dif_normalized_Intact-Slope_O2_dif_normalized_Anesthetized)/Speed_mm_s, 
                        COT.rel.ml = -(Slope_O2_dif_normalized_Intact-Slope_O2_dif_normalized_Anesthetized)/BLperSecond) %>% 
   mutate(COT.abs.mgC = -(Slopes_dif_mgC_Intact-Slopes_dif_mgC_Anesthetized)/Speed_mm_s, 
          COT.rel.mgC = -(Slopes_dif_mgC_Intact-Slopes_dif_mgC_Anesthetized)/BLperSecond) %>% 
   mutate(COT.p.ml = 100 * (Slope_O2_dif_normalized_Anesthetized - Slope_O2_dif_normalized_Intact) / pmax(abs(Slope_O2_dif_normalized_Intact), 1e-6)) %>%
   mutate(COT.p.ml = pmax(-100, pmin(100, COT.p.ml))) %>%  
          mutate(COT.p.mgC = -100*(Slopes_dif_mgC_Anesthetized-Slopes_dif_mgC_Intact)/Slopes_dif_mgC_Intact) %>% 
   mutate(COL.abs.ml = -Slope_O2_dif_normalized_Intact/Speed_mm_s, 
          COL.rel.ml = -Slope_O2_dif_normalized_Intact/BLperSecond, 
          COL.abs.mgC = -Slope_O2_dif_normalized_Intact/Speed_mm_s, 
          COL.rel.mgC = -Slope_O2_dif_normalized_Intact/BLperSecond)-> COT_species 
 
 COT_with_means %>% mutate(COT.abs.ml = -(Slope_O2_dif_normalized-Slope_O2_dif_normalized_Anesthetized)/Speed_mm_s, 
                        COT.rel.ml = -(Slope_O2_dif_normalized-Slope_O2_dif_normalized_Anesthetized)/BLperSecond) %>% 
   mutate(COT.abs.mgC = -(Slopes_dif_mgC-Slopes_dif_mgC_Anesthetized)/Speed_mm_s, 
          COT.rel.mgC = -(Slopes_dif_mgC-Slopes_dif_mgC_Anesthetized)/BLperSecond) %>% 
   mutate(COT.p.ml = 100 * (Slope_O2_dif_normalized_Anesthetized - Slope_O2_dif_normalized) / pmax(abs(Slope_O2_dif_normalized), 1e-6)) %>%
   mutate(COT.p.ml = pmax(-100, pmin(100, COT.p.ml))) %>%  
   mutate(COT.p.mgC = -100*(Slopes_dif_mgC_Anesthetized-Slopes_dif_mgC)/Slopes_dif_mgC) %>% 
   mutate(COL.abs.ml = -Slope_O2_dif_normalized/Speed_mm_s, 
          COL.rel.ml = -Slope_O2_dif_normalized/BLperSecond, 
          COL.abs.mgC = -Slope_O2_dif_normalized/Speed_mm_s, 
          COL.rel.mgC = -Slope_O2_dif_normalized/BLperSecond)-> COT_with_means 
 
#order species as factors by architecture
COT_species$Species %>% factor(levels=c("Pegea sp.", "Helicosalpa virgula","Cyclosalpa affinis",
                                        "Cyclosalpa bakeri", "Cyclosalpa quadriluminis","Cyclosalpa polae",
                                        "Cyclosalpa sewelli", "Ritteriella retracta", "Ritteriella amboinensis",
                                        "Brooksia rostrata", "Thalia sp.", "Metcalfina hexagona", "Salpa fusiformis", "Salpa aspera", "Salpa maxima",
                                        "Iasis cylindrica", "Ihlea punctata", "Soestia zonaria")) -> COT_species$Species

COT_with_means$Species %>% factor(levels=c("Pegea sp.", "Helicosalpa virgula","Cyclosalpa affinis",
                                        "Cyclosalpa bakeri", "Cyclosalpa quadriluminis","Cyclosalpa polae",
                                        "Cyclosalpa sewelli", "Ritteriella retracta", "Ritteriella amboinensis",
                                        "Brooksia rostrata", "Thalia sp.", "Metcalfina hexagona", "Salpa fusiformis", "Salpa aspera", "Salpa maxima",
                                        "Iasis cylindrica", "Ihlea punctata", "Soestia zonaria")) -> COT_with_means$Species

COT_species <- unique(COT_species)


full_join(COT_with_means, swim[,3:4] %>% unique(), by="Species") -> COT_with_means
COT_with_means$Architecture[which(COT_with_means$Species == "Ritteriella retracta")] <- "Bipinnate"
COT_with_means$Architecture[which(COT_with_means$Species == "Thalia sp.")] <- "Oblique"
COT_with_means$Architecture[which(COT_with_means$Species == "Helicosalpa virgula")] <- "Helical"
COT_with_means$Architecture[which(COT_with_means$Species == "Ihlea punctata")] <- "Linear"
COT_with_means$Architecture[which(COT_with_means$Species == "Cyclosalpa quadriluminis")] <- "Whorl"
COT_with_means$Architecture[which(COT_with_means$Species == "Cyclosalpa affinis")] <- "Whorl"

full_join(COT_species, swim[,3:4], by="Species") -> COT_species
COT_species$Architecture[which(COT_species$Species == "Ritteriella retracta")] <- "Bipinnate"
COT_species$Architecture[which(COT_species$Species == "Thalia sp.")] <- "Oblique"
COT_species$Architecture[which(COT_species$Species == "Helicosalpa virgula")] <- "Helical"
COT_species$Architecture[which(COT_species$Species == "Ihlea punctata")] <- "Linear"
COT_species$Architecture[which(COT_species$Species == "Cyclosalpa quadriluminis")] <- "Whorl"
COT_species$Architecture[which(COT_species$Species == "Cyclosalpa affinis")] <- "Whorl"

COT_species %>% unique() -> COT_species

#Net respiration rate by biovolume in Swimmers and KO
# COT_species %>% filter(!is.na(COT.abs.ml) & Architecture != "Whorl chain") %>% 
#   ggplot(aes(x=Species))+
#   geom_point(aes(y=(-Slope_O2_dif_normalized_Intact)),color="red", alpha=0.7)+
#   geom_point(aes(y=(-Slope_O2_dif_normalized_Anesthetized)), color="blue",alpha=0.7)+
#   ylab("Net respiration rate (pgO2/min) per specimen biovolume (ml)")+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90))
# 
# COT_paired %>% filter(-Slope_O2_dif_normalized_Intact<2500 & -Slope_O2_dif_normalized_Anesthetized<2500) %>% 
#   ggplot(aes(x=Species))+
#   geom_point(aes(y=(-Slope_O2_dif_normalized_Intact)),color="red", alpha=0.7)+
#   geom_point(aes(y=(-Slope_O2_dif_normalized_Anesthetized)), color="blue",alpha=0.7)+
#   ylab("Net respiration rate (pgO2/min) per specimen biovolume (ml)")+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90))
# 
# COT_with_means %>% filter(-Slope_O2_dif_normalized<2500 & -Slope_O2_dif_normalized_Anesthetized<2500) %>% 
#   ggplot(aes(x=Species))+
#   geom_point(aes(y=(-Slope_O2_dif_normalized)),color="red", alpha=0.7)+
#   geom_point(aes(y=(-Slope_O2_dif_normalized_Anesthetized)), color="blue",alpha=0.7)+
#   ylab("Net respiration rate (pgO2/min) per specimen biovolume (ml)")+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90))

####################


## Difference Swim-KO barplot for presentations ##

architecture_order <- c("Transversal","Oblique","Linear","Bipinnate","Helical","Whorl","Cluster")
COT_species %>% filter(-(Slope_O2_dif_normalized_Intact - Slope_O2_dif_normalized_Anesthetized)>1) %>% 
  ggplot(aes(x = factor(Species %>% as.character(), levels = unique(Species[order(factor(Architecture, levels = architecture_order))])))) +
  geom_bar(
    aes(y = -(Slope_O2_dif_normalized_Intact - Slope_O2_dif_normalized_Anesthetized), fill = Architecture),
    stat = "identity"
  ) +
  ylab("Differential respiration rate (pgO2/min/biovolume_ml) Swimming-K.O.") +
  scale_fill_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                               "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                               "Cluster" = "magenta")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

COT_with_means %>% filter(-(Slope_O2_dif_normalized - Slope_O2_dif_normalized_Anesthetized)>1) %>% 
  ggplot(aes(x = factor(Species %>% as.character(), levels = unique(Species[order(factor(Architecture, levels = architecture_order))])), 
             y = -(Slope_O2_dif_normalized - Slope_O2_dif_normalized_Anesthetized)))+
  stat_summary(fun = mean, geom = "point", aes(col = Architecture), size = 3) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = "dodge", width = 0.25,aes(col=Architecture)) +  # Add error bars
  ylab("Differential respiration rate (pgO2/min/biovolume_ml) Swimming-K.O.") +
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                               "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                               "Cluster" = "magenta")) +
  theme_bw() +
  xlab("Species")+
  theme(axis.text.x = element_text(angle = 90, hjust=1))

#COT by biovolume per cm across Architecture
COT_with_means %>% filter(Slope_O2_dif_normalized>-900) %>% 
  filter(!is.na(COT.abs.ml) & Architecture != "Whorl chain") %>% 
  ggplot(aes(x = factor(Species %>% as.character(), levels = unique(Species[order(factor(Architecture, levels = architecture_order))])),
             y = COT.abs.ml/Colony.volume..ml.^0.75, col = Architecture, fill=Architecture)) +
  stat_summary(fun = "mean", geom = "bar", position = "dodge") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = "dodge", width = 0.25) +  # Add error bars
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                               "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                               "Cluster" = "magenta")) +
  scale_fill_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Cost of Transport (pgO2/ml per mm moved)") + xlab("Species") 

COT_with_means %>%
  filter(!is.na(COT.abs.ml)  & Architecture != "Whorl chain") %>% 
  ggplot(aes(x = factor(Species %>% as.character(), levels = unique(Species[order(factor(Architecture, levels = architecture_order))])),
             y = COT.abs.ml*4.75, fill=Architecture)) +
  stat_summary(fun = "mean", geom = "bar", position = "dodge") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = "dodge", width = 0.25) +  # Add error bars
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  scale_fill_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                               "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                               "Cluster" = "magenta")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(breaks = seq(0, max(COT_with_means$COT.abs.ml * 4.75, na.rm=T), by = 5)) +  # Add this line
  ylab("Cost of Transport (J/kg/m)") + xlab("Species")

# ggplot(COT_species %>% filter(!is.na(COT.abs.ml) & Species!="Brooksia rostrata" & Architecture != "Whorl chain"),
#        aes(x=factor(Species %>% as.character(), levels = unique(Species[order(factor(Architecture, levels = architecture_order))])),
#            y=COT.abs.ml, fill=Architecture))+ 
#   stat_summary(fun = "mean", geom = "bar", position = "dodge") +
#   scale_fill_manual(values=c("cyan4","magenta","darkorange1","green4","darkorchid4"))+
#   ylab("Net cost of Transport (pgO2/ml per zooid length moved)") + 
#   theme_dark() +
#   theme(
#     panel.grid.major = element_blank(),  # Remove major grid lines
#     panel.grid.minor = element_blank(),  # Remove minor grid lines
#     axis.line = element_line(color = "white"),  # Set axis line color to black
#     text = element_text(color = "white"),  # Set text color to white
#     axis.text = element_text(color = "white"),  # Set axis text color to white
#     panel.background = element_rect(fill = "black"),  # Set panel background color to black
#     axis.text.x = element_text(angle = 90, hjust=1)
#   )

########################

#### Figure 5A #####
architecture_order <- c("Transversal", "Linear", "Bipinnate", "Whorl", "Cluster","Helical")

F5A <- ggplot(COT_with_means %>% 
                filter(COT.abs.ml>-1000000 & Architecture != "Whorl chain"), 
       aes(factor(Species %>% as.character(), levels = unique(Species[order(factor(Architecture, levels = architecture_order))])), 
           y=COT.abs.ml, fill=Architecture)) +
  geom_jitter(aes(color=Architecture),alpha=0.7)+
  stat_summary(fun = "mean", geom = "bar", position = "dodge", alpha=0.5) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = "dodge", width = 0.25) +  # Add error bars
  scale_fill_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                               "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                               "Cluster" = "magenta")) +
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust=1))+
  ylab("Cost of Transport (pgO2/ml per mm moved)") + guides(fill="none",color="none")+
  xlab("Species")

#### Figure 5B #####

#COT by biovolume per bodylength across Species
F5B <- ggplot(COT_with_means %>% 
                filter(COT.rel.ml>-10000000 & Architecture != "Whorl chain"), 
       aes(factor(Species %>% as.character(), levels = unique(Species[order(factor(Architecture, levels = architecture_order))])), 
           y=COT.rel.ml, fill=Architecture)) +
  geom_jitter(aes(color=Architecture),alpha=0.7)+
  stat_summary(fun = "mean", geom = "bar", position = "dodge", alpha=0.5) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = "dodge", width = 0.25) +  # Add error bars
  scale_fill_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                               "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                               "Cluster" = "magenta")) +
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust=1))+
  ylab("Cost of Transport (pgO2/ml per zooid length moved)")+ guides(fill="none",color="none")+
  xlab("Species")

#### Figure 5C #####
architecture_order <- c("Transversal", "Linear", "Bipinnate", "Whorl", "Cluster","Helical")

F5C <- ggplot(COT_with_means %>% 
                filter(COT.abs.ml>-10000000 & Architecture != "Whorl chain"), 
              aes(Architecture, 
                  y=COT.abs.ml, fill=Architecture)) +
  geom_jitter(aes(color=Architecture),alpha=0.7)+
  stat_summary(fun = "mean", geom = "bar", position = "dodge", alpha=0.5) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = "dodge", width = 0.25) +  # Add error bars
  scale_fill_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                               "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                               "Cluster" = "magenta")) +
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust=1))+
  ylab("Cost of Transport (pgO2/ml per zooid length moved)")+
  xlab("Architecture")

#### Figure 5D #####

#COT by biovolume per bodylength across Species
F5D <- ggplot(COT_with_means %>% 
                filter(COT.rel.ml>-10000000 & Architecture != "Whorl chain"), 
              aes(Architecture, 
                  y=COT.rel.ml, fill=Architecture)) +
  geom_jitter(aes(color=Architecture),alpha=0.9)+
  stat_summary(fun = "mean", geom = "bar", position = "dodge", alpha=0.5) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = "dodge", width = 0.25) +  # Add error bars
  scale_fill_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                               "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                               "Cluster" = "magenta")) +
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                               "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                               "Cluster" = "magenta")) +
  theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust=1))+
  ylab("Cost of Transport (pgO2/ml per zooid length moved)")+
  xlab("Architecture")

wrap_plots(F5A, F5B, F5C, F5D)

anova_COTa_arch <- aov(COT.abs.ml ~ Architecture, data=COT_with_means %>% filter(!(Architecture %in% c("Whorl chain","Helical"))))
summary(anova_COTa_arch)
TukeyHSD(anova_COTa_arch, conf.level=.95) -> tukey_COTa_arch

anova_COTr_arch <- aov(COT.rel.ml ~ Architecture, data=COT_with_means %>% filter(!(Architecture %in% c("Whorl chain","Helical"))))
summary(anova_COTr_arch)
TukeyHSD(anova_COTr_arch, conf.level=.95) -> tukey_COTr_arch

####################

#### Figure 6A ####

#COT by biovolume per mm across absolute speeds
COT_aggregated_F6A <- COT_with_means %>%
  filter(COT.abs.ml > -10, Architecture != "Whorl chain") %>%
  group_by(Species, Architecture) %>%
  summarize(
    mean_COT = mean(COT.abs.ml, na.rm = TRUE),
    se_COT = sd(COT.abs.ml, na.rm = TRUE) / sqrt(n()),
    mean_Speed_mm_s = mean(Speed_mm_s, na.rm = TRUE),
    .groups = "drop"
  )

F6A <- ggplot() +
  # Raw data layer for individual points, ensuring we specify the full dataset
  geom_point(data = COT_with_means %>% 
               filter(COT.abs.ml > -10, Architecture != "Whorl chain"), 
             aes(x = Speed_mm_s, y = COT.abs.ml, col = Architecture), 
             alpha = 0.2) +
  
  # Smoothed line layer
  geom_smooth(data = COT_with_means %>% 
                filter(COT.abs.ml > -10, Architecture != "Whorl chain"),
              aes(x = Speed_mm_s, y = COT.abs.ml), 
              method = "lm", color = "black", formula = y ~ log(x)) +
  
  # Mean points and error bars from pre-aggregated data
  geom_point(data = COT_aggregated_F6A, 
             aes(x = mean_Speed_mm_s, y = mean_COT, col = Architecture), 
             size = 3) +
  geom_errorbar(data = COT_aggregated_F6A, 
                aes(x = mean_Speed_mm_s, y = mean_COT, ymin = mean_COT - se_COT, 
                    ymax = mean_COT + se_COT, col = Architecture), 
                width = 0.2) +
  
  # Color, labels, and theme
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", 
                                "Linear" = "darkorange1", "Bipinnate" = "cyan4", 
                                "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  ylab("Cost of Transport (pgO2/ml per mm moved)") +
  xlab("Speed (mm/s)") +
  theme_bw() +
  ylim(0, NA) +
  guides(color = "none")

##### Figure 6B #####

#COT by biovolume per cm across relative speeds

F6B <- ggplot(COT_with_means %>% 
                filter(COT.abs.ml > -10 & Architecture != "Whorl chain" & !is.na(Speed_mm_s) & !is.na(Species)) %>% group_by(Species), 
              aes(x = BLperSecond, y = COT.rel.ml)) +
  geom_smooth(method = "lm", color = "black", formula = y ~ log(x)) +
  stat_summary(fun = mean, geom = "point", aes(col = Architecture, group = Species), size = 2.7, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(col = Architecture, group = Species), width = 0.2, alpha = 0.7) +
  geom_point(data = COT_with_means %>% filter(COT.abs.ml > -10 & Architecture != "Whorl chain"), 
              aes(x = BLperSecond, y = COT.rel.ml, col = Architecture), 
              width = 0.2, height = 0, alpha = 0.2) +
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  ylab("Cost of Transport (pgO2/ml per zooid length moved)") + 
  xlab("Speed (zooids/s)") + 
  theme_bw() + 
  ylim(0, NA)

wrap_plots(F6A, F6B)

lm(log(COT.abs.ml) ~ Speed_mm_s, data=COT_with_means %>% filter(!is.na(COT.abs.ml) & Architecture != "Whorl chain")) %>% summary()
lm(log(COT.rel.ml) ~ BLperSecond, data=COT_with_means %>% filter(!is.na(COT.rel.ml) & Architecture != "Whorl chain")) %>% summary()

lm(COT.abs.ml ~ Speed_mm_s, data=COT_with_means %>% filter(!is.na(COT.abs.ml) & Architecture != "Whorl chain")) %>% summary()
lm(COT.rel.ml ~ BLperSecond, data=COT_with_means %>% filter(!is.na(COT.rel.ml) & Architecture != "Whorl chain")) %>% summary()

#########################

# COT_species <- mutate(COT_species, percentSwim = -5000*(Slope_O2_dif_normalized_Intact - Slope_O2_dif_normalized_Anesthetized) )
# 
# ggplot(COT_species %>% filter(!is.na(percentSwim)), aes(x=Speed.cm.s,y=percentSwim))+
#   geom_point()+
#   geom_text(label=COT_species %>% filter(!is.na(percentSwim)) %>% .$Species)+
#   theme_bw()

#COT by biovolume per zooid length across relative speeds
ggplot(COT_species %>% filter(!is.na(COT.rel.ml)), aes(x=BLperSecond,y=log(COT.rel.ml)))+
  geom_point(aes(color=BLperSecond))+
  theme_bw()+
  geom_text(label=COT_species %>% filter(!is.na(COT.rel.ml)) %>% .$Species, hjust=0.4, vjust=-0.7)+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Cost of Transport (pgO2/ml per body length moved)")

#### SM Figure 6 ####

# % cost (by biovolume) invested in swimming across species order: COT.p, color: Speed_cm
architecture_order <- c("Transversal", "Oblique", "Linear", "Bipinnate", "Helical", "Whorl", "Cluster")

ggplot(COT_with_means %>% filter(!is.na(COT.p.ml) & Species != "Brooksia rostrata" & Architecture != "Whorl chain"),
       aes(x = factor(Species %>% as.character(), levels = unique(Species[order(factor(Architecture, levels = architecture_order))])), 
           y = COT.p.ml, fill = Architecture)) +
  stat_summary(fun = "mean", geom = "bar", position = "dodge") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = "dodge", width = 0.25) +  # Add error bars
  scale_fill_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                               "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                               "Cluster" = "magenta")) +
  theme_bw() +
  xlab("Species") +
  ylab("Proportion of Metabolic Cost Spent on Swimming (%)") +theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

###########################

###### SM FIGURE 7 ####### # % cost (by biovolume) invested in swimming across pulstion rates , color: Species

COT_means <- COT_with_means %>%
  filter(COT.p.ml>0 & !is.na(Species) & Architecture != "Whorl chain") %>%
  group_by(Species, Architecture) %>%
  summarise(mean_COT = mean(COT.p.ml),
            mean_Pulses_per_second = mean(Pulses_per_second))

ggplot(COT_with_means %>% 
         filter(COT.p.ml>0 & !is.na(Species) & Architecture != "Whorl chain" & Species != "Brooksia rostrata"), 
       aes(x = Pulses_per_second, y = COT.p.ml)) +
  stat_summary(fun = mean, geom = "point", aes(col = Architecture), size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(col = Architecture), width = 0.2) +
  geom_smooth(method = "lm", color = "black") +
  geom_text(data = COT_means %>% filter(Architecture %in% unique(COT_with_means$Architecture)), 
            aes(label = Species, x = mean_Pulses_per_second, y = mean_COT), 
            hjust = 0.4, vjust = -0.7, color = "black", cex = 3) +
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  xlab("Pulsation rate (pulses/s)") +
  ylab("Proportion of Metabolic Cost Spent on Swimming (%)") + 
  theme_bw()

lm(COT.p.ml ~ Pulses_per_second, 
    data = COT_with_means %>% 
      filter(COT.p.ml>0 & !is.na(Species) & Architecture != "Whorl chain"), 
    family = gaussian(link = "identity")) %>% summary()

COT_with_means %>% filter(!is.na(COT.p.ml) & Species != "Brooksia rostrata" & Architecture != "Whorl chain") %>% 
  .[,c("Species","COT.p.ml")] %>% group_by(Species) %>% 
  summarise(mean(COT.p.ml))


#############

glm(COT.p.ml ~ Speed_mm_s, data = COT_with_means %>%  
      filter(!is.na(COT.p.ml) & !is.na(Species) & Architecture != "Whorl chain"), 
    family = gaussian(link = "identity")) %>% summary()

glm(COT.p.ml ~ BLperSecond, data = COT_with_means %>%  
      filter(!is.na(COT.p.ml) & !is.na(Species) & Architecture != "Whorl chain"), 
    family = gaussian(link = "identity")) %>% summary()

##### SM FIGURE 8 #####

# % cost (by biovolume) invested in swimming across COT_vol_mm , color: Species

Sm8A <- ggplot(COT_with_means %>%  filter(!is.na(COT.p.ml) & !is.na(Species)), 
               aes(y = COT.abs.ml,x=COT.p.ml))+
  stat_summary(fun = mean, geom = "point", aes(col = Architecture), size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(col = Architecture), width = 0.2) +
  geom_smooth(method="lm", color="black")+
 # geom_text(label=COT_species %>% filter(!is.na(COT.p.ml) & !is.na(Species) & Architecture != "Whorl chain" & Species != "Thalia sp.") %>% 
  #            .$Species, hjust=0.4, vjust=-0.7, color="black")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  xlim(0,100)+
  ylab("Cost of Transport (pgO2/mgC per mm moved)")+
  xlab("Proportion of Metabolic Cost Spent on Swimming (%)") + guides(color="none")

# % cost (by biovolume) invested in swimming across COT_vol_mm , color: Species

Sm8B <- ggplot(COT_with_means %>%  filter(!is.na(COT.p.ml) & !is.na(Species)), aes(y = COT.rel.ml,x=COT.p.ml))+
  stat_summary(fun = mean, geom = "point", aes(col = Architecture), size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(col = Architecture), width = 0.2) +
  geom_smooth(method="lm", color="black")+
  theme_bw()+
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  xlim(0,100)+
  ylab("Cost of Transport (pgO2/mgC per zooid length moved)")+
  xlab("Proportion of Metabolic Cost Spent on Swimming (%)")

wrap_plots(Sm8A, Sm8B)

lm(COT.p.ml ~ COT.abs.ml, data = COT_with_means %>%  
      filter(!is.na(COT.p.ml) & !is.na(Species)), 
    family = gaussian(link = "identity")) %>% summary()

lm(COT.p.ml ~ COT.rel.ml, data = COT_with_means %>%  
      filter(!is.na(COT.p.ml) & !is.na(Species)), 
    family = gaussian(link = "identity")) %>% summary()

##########################

##### SM FIGURE 9 #####

# % cost (by biovolume) invested in swimming across COT_vol_mm , color: Species

Sm9A <- ggplot(COT_with_means %>%  filter(!is.na(COT.p.ml) & !is.na(Species)), 
               aes(y = Speed_mm_s,x=COT.p.ml))+
  stat_summary(fun = mean, geom = "point", aes(col = Architecture), size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(col = Architecture), width = 0.2) +
  geom_smooth(method="lm", color="black")+
  # geom_text(label=COT_species %>% filter(!is.na(COT.p.ml) & !is.na(Species) & Architecture != "Whorl chain" & Species != "Thalia sp.") %>% 
  #            .$Species, hjust=0.4, vjust=-0.7, color="black")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  xlim(0,100)+
  ylab("Swimming speed (mm/s)")+
  xlab("Proportion of Metabolic Cost Spent on Swimming (%)") + guides(color="none")

# % cost (by biovolume) invested in swimming across COT_vol_mm , color: Species

Sm9B <- ggplot(COT_with_means %>%  filter(!is.na(COT.p.ml) & !is.na(Species)), aes(y = BLperSecond,x=COT.p.ml))+
  stat_summary(fun = mean, geom = "point", aes(col = Architecture), size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(col = Architecture), width = 0.2) +
  geom_smooth(method="lm", color="black")+
  theme_bw()+
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  xlim(0,100)+
  ylab("Swimming speed (zooids/s)")+
  xlab("Proportion of Metabolic Cost Spent on Swimming (%)")

wrap_plots(Sm9A, Sm9B)

lm(COT.p.ml ~ Speed_mm_s, data = COT_with_means %>%  
      filter(!is.na(COT.p.ml) & !is.na(Species)), 
    family = gaussian(link = "identity")) %>% summary()

lm(COT.p.ml ~ BLperSecond, data = COT_with_means %>%  
      filter(!is.na(COT.p.ml) & !is.na(Species)), 
    family = gaussian(link = "identity")) %>% summary()

##################

# % cost (by biovolume) invested in swimming across basal rate , color: Species
ggplot(COT_species %>%  filter(!is.na(COT.p.ml) & Species != "Thalia sp."), aes(x = -Slope_O2_dif_normalized_Anesthetized,y=COT.p.ml))+
  geom_point(aes(size=Zooid.length..mm.,color=Species))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  #ylim(c(0,100))+
  ylab("% Cost of Transport Volume-normalized")+xlab("Basal rate")

## SUMMARY TABLES ##

COT_with_means %>% filter(!is.na(Species)) %>% 
  group_by(Species) %>%
  summarise(
    Nspecimens = n_distinct(Specimen),
    Nmeasurements = n()
  ) 

slopes %>% filter(-Slope_O2_dif_normalized/1000000 < 0.0025) %>% 
  group_by(Species,Treatment) %>% 
  summarise(SLN = mean(-Slope_O2_dif_normalized, na.rm=T)) %>% as.data.frame()

slopes %>% filter(-Slope_O2_dif_normalized/1000000 < 0.0025) %>% 
  group_by(Species) %>% 
  summarise(CVol = mean(Zooid.length..mm., na.rm=T)) %>% as.data.frame()

COT_with_means %>% filter(-Slope_O2_dif_normalized/1000000 < 0.0025) %>% 
  group_by(Architecture) %>% 
  summarise(meanCOTabs = mean(COT.abs.ml, na.rm=T)) %>% as.data.frame()

COT_with_means %>% filter(-Slope_O2_dif_normalized/1000000 < 0.0025) %>% 
  group_by(Species) %>% 
  summarise(meanCOTp = mean(COT.p.ml, na.rm=T)) %>% as.data.frame()


