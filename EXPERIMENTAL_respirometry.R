
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

slopes %>% filter(Slope_O2_dif<0) -> slopes
slopes$Slope_O2_dif = slopes$Slope_O2_dif-0.000967
slopes$Slope_O2_dif_normalized = slopes$Slope_O2_dif_normalized-802.86

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


## SM Figure 4 ##

slopes %>%  ggplot(aes(x=Species,y=-Slope_O2_dif_normalized/1000000))+
  geom_boxplot(aes(fill=Treatment))+
  ylab("Respiration rate (mgO2/min) per specimen biovolume (ml)")+
  ylim(NA,0.002)+
  scale_fill_manual(values = c("lightblue","pink"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))

swim <- read.csv("~/salp_ecomorphology/EMSpeeds_final_annotated.tsv",stringsAsFactors = F) #[,c(20,21,25)]
summarized_raw=swim


energetics <- mutate(slopes, Species=as.character(Species))
energetics$Species[energetics$Species=="Iasis (Weelia) cylindrica"] <- "Iasis cylindrica"
summarized_raw$Species[summarized_raw$Species=="Ritteriella sp."] <- "Ritteriella retracta"

summarized_raw[,c(3,5,4,9,10)]-> swimming_data
swimming_data %>% group_by(Species) %>% summarise(across(where(is.numeric), mean, na.rm = TRUE)) -> swimming_species

energetics_merged <- left_join(energetics, swimming_species, by="Species")

energetics_merged %>% mutate(Speed_mm_s = BLperSecond*Zooid.length..mm.) -> energetics_merged
energetics_merged %>% mutate(Speed.bl.p = BLperSecond/Pulses_per_second) -> energetics_merged

COT <- energetics_merged[,c("Species","Specimen","Colony.volume..ml.","Zooid.length..mm.","Number.of.zooids",
                            "Treatment","Paired","Slope_O2_dif_normalized","Slopes_dif_mgC","Speed_mm_s",
                            "BLperSecond","Pulses_per_second", "Speed.bl.p")]

##### WARNING ##### SKETCHY STUFF TO REMOVE NEGATIVES !!!!
#COT$Slope_O2_dif_normalized[which(COT$Slope_O2_dif_normalized>0)] <- 0

COT %>% filter(Paired=="Yes") %>% 
  pivot_wider(names_from = Treatment, values_from = c(Slope_O2_dif_normalized, Slopes_dif_mgC)) %>% 
  mutate(COT.abs.ml = -(Slope_O2_dif_normalized_Intact-Slope_O2_dif_normalized_Anesthetized)/Speed_mm_s, 
         COT.rel.ml = -(Slope_O2_dif_normalized_Intact-Slope_O2_dif_normalized_Anesthetized)/BLperSecond) %>% 
  mutate(COT.abs.mgC = -(Slopes_dif_mgC_Intact-Slopes_dif_mgC_Anesthetized)/Speed_mm_s, 
         COT.rel.mgC = -(Slopes_dif_mgC_Intact-Slopes_dif_mgC_Anesthetized)/BLperSecond) -> COT_paired


COT %>% filter(!is.na(Species)) %>% 
  pivot_wider(names_from = Treatment, values_from = c(Slope_O2_dif_normalized, Slopes_dif_mgC)) %>% 
  group_by(Species) %>% 
  summarise_at(vars("Speed_mm_s", "BLperSecond", "Slope_O2_dif_normalized_Intact", 
                    "Slope_O2_dif_normalized_Anesthetized", "Slopes_dif_mgC_Intact", 
                    "Slopes_dif_mgC_Anesthetized", "Zooid.length..mm.","Pulses_per_second"),
               function(x){mean(x, na.rm=T)}) %>% as.data.frame() -> COT_species

# Merge the mean values with the original data
COT_with_means <- COT %>%
  left_join(COT_species %>% select(Species, Slope_O2_dif_normalized_Anesthetized, Slopes_dif_mgC_Anesthetized), by = "Species")

COT_with_means %>% mutate(DIF = Slope_O2_dif_normalized - Slope_O2_dif_normalized_Anesthetized) %>% .$DIF %>%  max(na.rm=T)
COT_with_means[which(COT_with_means$Slope_O2_dif_normalized > COT_with_means$Slope_O2_dif_normalized_Anesthetized),"Slope_O2_dif_normalized"]<-COT_with_means[which(COT_with_means$Slope_O2_dif_normalized > COT_with_means$Slope_O2_dif_normalized_Anesthetized),"Slope_O2_dif_normalized"]-2016.389


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


#COT_species <- filter(COT_species, Species != "Cyclosalpa bakeri")  ### REMOVED OUTLIER WONKY SPP !!!
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

dim(COT_with_means)
dim(COT_with_means %>% filter(Slope_O2_dif_normalized_Anesthetized>Slope_O2_dif_normalized))

COT_with_means %>% filter(Slope_O2_dif_normalized_Anesthetized>Slope_O2_dif_normalized) -> COT_with_means



