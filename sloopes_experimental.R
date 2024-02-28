
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
slopes %>% mutate(Slope_normalized = 1000000*Slope_O2/(Colony.volume..ml.^0.75), Slope_O2_dif_normalized = 1000000*Slope_O2_dif/(Colony.volume..ml.^0.75)) -> slopes

slopes$Species[which(slopes$Species == "Thalia cicar")] <- "Thalia sp."
slopes$Species[which(slopes$Species == "Ritteriella sp.")] <- "Ritteriella retracta"
slopes$Species %>% factor(levels=c("Pegea sp.", "Helicosalpa virgula","Cyclosalpa affinis", 
                                   "Cyclosalpa bakeri", "Cyclosalpa quadriluminis","Cyclosalpa polae",
                                   "Cyclosalpa sewelli", "Ritteriella retracta", "Ritteriella amboinensis",
                                   "Brooksia rostrata", "Thalia sp.", "Metcalfina hexagona", "Salpa fusiformis", "Salpa aspera", "Salpa maxima", 
                                   "Iasis (Weelia) cylindrica", "Ihlea punctata", "Soestia zonaria")) -> slopes$Species

slopes %>%  ggplot(aes(x=Species,y=-Slope_O2_dif_normalized/1000000))+
  geom_boxplot(aes(fill=Treatment))+
  ylab("Respiration rate (mgO2/min) per specimen biovolume (ml)")+
  scale_fill_manual(values = c("lightblue","pink"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))

slopes %>%  ggplot(aes(x=Species,y=-Slope_O2_dif))+
  geom_boxplot(aes(fill=Treatment))+
  ylab("Respiration rate (pgO2/min)")+
  scale_fill_manual(values = c("lightblue","pink"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))

slopes %>%  ggplot(aes(x=Species,y=-60*31.25*Slope_O2_dif/1000000))+
  geom_boxplot(aes(fill=Treatment))+
  ylab("Respiration rate (µmolO2/h)")+
  scale_fill_manual(values = c("lightblue","pink"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))

#Get swimming speeds on Event Measure
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
                            "Treatment","Paired","Slope_O2_dif","Slopes_dif_mgC","Speed_mm_s",
                            "BLperSecond","Pulses_per_second", "Speed.bl.p")]

COT %>% filter(!is.na(Species)) %>% 
  pivot_wider(names_from = Treatment, values_from = c(Slope_O2_dif, Slopes_dif_mgC)) %>% 
  group_by(Species) %>% 
  summarise_at(vars("Speed_mm_s", "BLperSecond", "Slope_O2_dif_Intact", 
                    "Slope_O2_dif_Anesthetized", "Slopes_dif_mgC_Intact", 
                    "Slopes_dif_mgC_Anesthetized", "Zooid.length..mm.","Pulses_per_second"),
               function(x){mean(x, na.rm=T)}) %>% as.data.frame() -> COT_species

COT_species %>% filter(Slope_O2_dif_Anesthetized>Slope_O2_dif_Intact) -> COT_species

# Merge the mean values with the original data
COT_with_means <- COT %>%
  left_join(COT_species %>% select(Species, Slope_O2_dif_Anesthetized, Slopes_dif_mgC_Anesthetized), by = "Species")

COT_with_means %>% filter(Slope_O2_dif_Anesthetized>Slope_O2_dif) -> COT_with_means

COT_with_means %>% mutate(COT.abs.ml = -(Slope_O2_dif-Slope_O2_dif_Anesthetized)/Speed_mm_s, 
                          COT.rel.ml = -(Slope_O2_dif-Slope_O2_dif_Anesthetized)/BLperSecond) %>% 
  mutate(COT.p.ml = 100 * (Slope_O2_dif_Anesthetized - Slope_O2_dif) / pmax(abs(Slope_O2_dif), 1e-6)) %>%
  mutate(COT.p.ml = pmax(-100, pmin(100, COT.p.ml))) %>%  
  mutate(COT.p.mgC = -100*(Slopes_dif_mgC_Anesthetized-Slopes_dif_mgC)/Slopes_dif_mgC)-> COT_with_means 

COT_with_means$Species %>% factor(levels=c("Pegea sp.", "Helicosalpa virgula","Cyclosalpa affinis",
                                           "Cyclosalpa bakeri", "Cyclosalpa quadriluminis","Cyclosalpa polae",
                                           "Cyclosalpa sewelli", "Ritteriella retracta", "Ritteriella amboinensis",
                                           "Brooksia rostrata", "Thalia sp.", "Metcalfina hexagona", "Salpa fusiformis", "Salpa aspera", "Salpa maxima",
                                           "Iasis cylindrica", "Ihlea punctata", "Soestia zonaria")) -> COT_with_means$Species

full_join(COT_with_means, swim[,3:4] %>% unique(), by="Species") -> COT_with_means
COT_with_means$Architecture[which(COT_with_means$Species == "Ritteriella retracta")] <- "Bipinnate"
COT_with_means$Architecture[which(COT_with_means$Species == "Thalia sp.")] <- "Oblique"
COT_with_means$Architecture[which(COT_with_means$Species == "Helicosalpa virgula")] <- "Helical"
COT_with_means$Architecture[which(COT_with_means$Species == "Ihlea punctata")] <- "Linear"
COT_with_means$Architecture[which(COT_with_means$Species == "Cyclosalpa quadriluminis")] <- "Whorl"
COT_with_means$Architecture[which(COT_with_means$Species == "Cyclosalpa affinis")] <- "Whorl"

#COT per cm across Architecture
COT_with_means %>%
  filter(!is.na(COT.abs.ml)  & Architecture != "Whorl chain") %>% 
  ggplot(aes(x = factor(Species %>% as.character(), levels = unique(Species[order(factor(Architecture, levels = architecture_order))])),
             y = COT.abs.ml, col = Architecture, fill=Architecture)) +
  stat_summary(fun = "mean", geom = "bar", position = "dodge") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = "dodge", width = 0.25) +  # Add error bars
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  scale_fill_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                               "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                               "Cluster" = "magenta")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Cost of Transport (pgO2 per mm moved)") + xlab("Species")

COT_with_means %>%
  filter(!is.na(COT.rel.ml)  & Architecture != "Whorl chain") %>% 
  ggplot(aes(x = factor(Species %>% as.character(), levels = unique(Species[order(factor(Architecture, levels = architecture_order))])),
             y = COT.rel.ml, col = Architecture, fill=Architecture)) +
  stat_summary(fun = "mean", geom = "bar", position = "dodge") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = "dodge", width = 0.25) +  # Add error bars
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  scale_fill_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                               "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                               "Cluster" = "magenta")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Cost of Transport (pgO2 per zooid length moved)") + xlab("Species")

##

ggplot(COT_with_means %>% 
         filter(COT.abs.ml>-10 & Architecture != "Whorl chain"), 
       aes(x=Colony.volume..ml.,y=COT.abs.ml)) +
  stat_summary(fun = mean, geom = "point", aes(col = Architecture), size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(col = Architecture), width = 0.2) +
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  ylab("Cost of Transport (pgO2 per mm moved)") +
  xlab("Colony volume (ml)")+ 
  theme_bw()

ggplot(COT_with_means %>% 
         filter(COT.rel.ml>-10 & Architecture != "Whorl chain"), 
       aes(x=Colony.volume..ml.,y=COT.rel.ml)) +
  stat_summary(fun = mean, geom = "point", aes(col = Architecture), size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(col = Architecture), width = 0.2) +
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  ylab("Cost of Transport (pgO2 per body length moved)") +
  xlab("Colony volume (ml)")+ 
  theme_bw()

COT_with_means %>%
  filter(!is.na(COT.abs.ml)  & Architecture != "Whorl chain") %>% filter(Species!="Thalia sp.") %>% 
  ggplot(aes(x = factor(Species %>% as.character(), levels = unique(Species[order(factor(Architecture, levels = architecture_order))])),
             y = 5*COT.abs.ml/Colony.volume..ml.*1025, fill=Architecture)) +
  stat_summary(fun = "mean", geom = "bar", position = "dodge") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = "dodge", width = 0.25) +  # Add error bars
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  scale_fill_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                               "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                               "Cluster" = "magenta")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Cost of Transport (J/kg/m)") + xlab("Species")
