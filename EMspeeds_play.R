# BACKBONE CODE FOR EM SPEEDS #

## Install packages ##

# install.packages("tidyverse")
# install.packages("magrittr")
# install.packages("patchwork")
# install.packages("data.table")

## Load packages ##

library(tidyverse)
library(magrittr)
library(patchwork)
require(data.table)

## Acquire list of filenames with speed data in your folder ##

list_tsv_files <- list.files(path = "~/Downloads/EMspeedsKai_July24_2023/") #change path to the folder where you store them

## Set working directory to that folder ##

setwd("~/Downloads/EMspeedsKai_July24_2023/") #change path to your folder

## Read the tab-separated-data data skipping the first 4 rows for every file in the folder ##

df = do.call(rbind, lapply(list_tsv_files, function(x) read.csv(x, sep='\t', header = T, stringsAsFactors = F, skip=4)))

## Identify the junk particle rows by the zooid number zero AND chop off the columns that are not useful ##

df[,c(1:3,6:8,31)] %>% filter(Number %in% c("0","1")) %>% mutate(is_junk = (Number=="0") == FALSE) -> keyvars #is_junk identifies the "junk" reference particle

## Make a tidy data frame by splitting the data into zooid and junk rows and then full-joining them side by side ##

tidy <- full_join(group_split(keyvars, is_junk)[[1]][,1:7], 
                  group_split(keyvars, is_junk)[[2]][,1:7], 
                  by=c("Filename","Frame","Time..mins."), 
                  suffix= c("Zooid","Junk")) %>% as.data.frame() %>% .[,c(-11)]

## Define a function to calculate time lags between rows to get the delta of time between measures, 
## get the vectors between the zooid and the junk particle, 
## calculate XYZ components of the distance, 
## use pythagoras to estimate the hypothenuse and divide by the delta of time

wrangle <- function(df){
  df %>% arrange(Frame) %>%
    mutate(Delta_time_s = (Time..mins. - lag(Time..mins., default = first(Time..mins.)))*60) %>% 
    mutate(Xd_mm = X..mm.Zooid - X..mm.Junk, Yd_mm = Y..mm.Zooid - Y..mm.Junk, Zd_mm =  Z..mm.Zooid - Z..mm.Junk) %>% 
    mutate(Delta_X_mm = Xd_mm - lag(Xd_mm, default = first(Xd_mm)), 
           Delta_Y_mm = Yd_mm - lag(Yd_mm, default = first(Yd_mm)),  
           Delta_Z_mm = Zd_mm - lag(Zd_mm, default = first(Zd_mm))) %>% 
    mutate(Distance_mm = sqrt(Delta_X_mm^2 + Delta_Y_mm^2 + Delta_Z_mm^2)) %>% 
    mutate(Speed_mm_s = Distance_mm/Delta_time_s) %>% filter(Speed_mm_s>0) %>% return()
}

## Run the wrangle function for the first (frontal) zooid ##

speed_data <- group_split(tidy, Filename) %>% lapply(wrangle) %>% do.call(rbind, .) %>% as.data.frame()

## Get video inventory metadata and match it by file name to the speeds to obtain: Species, Architecture, & Number of zooids ##

inventory <- read.csv("../Stereo_inventory_24July23.tsv", sep='\t', header = T, stringsAsFactors = F) #change path to where you saved it
names(inventory)[3] <- "Filename"
inventory$Species %>% str_replace_all("Iasis cylindrica \\(yellowtail\\)","Iasis cylindrica") -> inventory$Species
inventory$Species %>% str_replace_all("Metcalfina hexagona\\?","Metcalfina hexagona") -> inventory$Species
inventory$Species %>% str_replace_all("Soestia zonaria\\?","Soestia zonaria") -> inventory$Species
inventory$Species %>% str_replace_all("cf. ","") -> inventory$Species
inventory_clean <- inventory[,c(3,5,6,7,15)] #Keeping columns: filename, species, architecture, number of zooids

## Full join speeds with inventory metadata ##

speed_data$Filename %>% str_replace_all(".MOV", "") %>% str_replace_all(" - Copy", "") -> speed_data$Filename
speed_annotated <- full_join(speed_data, inventory_clean, by = "Filename" ) %>% 
  filter(Speed_mm_s != Inf) %>% filter(Species != "NA") %>%
  mutate(Speed_mms_abs = abs(Speed_mm_s)) #transform negative vectorial speeds into absolute scalars

## Load morphometrics data ##

mm <- read.csv("../ABPmorphometrics_18Sept23.tsv", sep='\t', header = T, stringsAsFactors = F) #change path to your file
mm <- mm[,c(1,11:14,17:19)]
names(mm)[1] <- "Filename"
coalesce_by_column <- function(df) {
  return(dplyr::coalesce(!!! as.list(df)))
}
mm %>% group_by(Filename) %>% summarise_all(coalesce_by_column) -> mm
mm <- filter(mm, Filename != "")

## Full-join morphometrics annotations to the speeds ##

speed_annotated <- full_join(speed_annotated, mm, by = "Filename")

## Load pulsation rates ##

pulsation <- read.csv("../PulsationRates_Oct16_23.tsv", sep='\t', header = T, stringsAsFactors = F) #change path to your file
pulsation <- pulsation[,c(1,16)]
names(pulsation) <- c("Filename", "Pulses_per_second")

## Full join pulsation rates to the speeds ##

speed_annotated <- full_join(speed_annotated, pulsation, by = "Filename" )

## Estimate normalized speed in body_lengths per pulse ##

speed_annotated <- mutate(speed_annotated, BLperSecond = Speed_mms_abs/Zooid_length_mm )
speed_annotated <- mutate(speed_annotated, BLperPulse = BLperSecond/Pulses_per_second )

## Collapse into the mean per specimen

speed_annotated[,c(1,19:22,24:34)] %>% unique() -> speed_collapsed
speed_collapsed %>% group_by(Filename, Species, Architecture) %>% 
  summarise(across(.cols = c(Speed_mm_s, Zooid.number, Speed_mms_abs, Zooid_length_mm, Pulses_per_second, BLperSecond, BLperPulse), 
                   ~ mean(., na.rm = TRUE))) %>% filter(!is.na(Species)) %>% as.data.frame() -> speed_collapsed

## Plot 


  #KAW

speed_collapsed %>% filter(!is.na(Species)) %>%  ggplot(aes(x=Species, y=Speed_mms_abs)) + geom_boxplot(aes(fill=Architecture))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90)) + ylab("Speed (mm/s)")

speed_collapsed %>% filter(!is.na(Species)) %>%  ggplot(aes(x=Species, y=BLperPulse)) + geom_boxplot(aes(fill=Architecture))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90)) + ylab("Speed (body lengths/pulsation)")

speed_collapsed %>% filter(!is.na(Architecture)) %>%  ggplot(aes(x=Architecture, y=Speed_mms_abs)) + geom_boxplot(aes(fill=Architecture))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90)) + ylab("Speed (mm/s)")

speed_collapsed %>% filter(!is.na(Architecture)) %>%  ggplot(aes(x=Architecture, y=BLperPulse)) + geom_boxplot(aes(fill=Architecture))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90)) + ylab("Speed (body lengths/pulsation)")

speed_collapsed %>% filter(!is.na(Species)) %>%  ggplot(aes(x=Pulses_per_second, y=Speed_mms_abs)) + geom_point(aes(col=Architecture))+
  theme_bw() + ylab("Speed (mm/s)") + xlab("Pulsation rate (pulses/s)") + geom_abline()

speed_collapsed %>% filter(!is.na(Species)) %>%  ggplot(aes(x=Pulses_per_second, y=BLperSecond)) + geom_point(aes(col=Architecture))+
  theme_bw() + ylab("Speed (body lengths/s)") + xlab("Pulsation rate (pulses/s)") + geom_abline()

glm(Speed_mms_abs ~ Pulses_per_second, data = speed_collapsed, family = gaussian(link = "identity")) %>% summary()
glm(BLperSecond ~ Pulses_per_second, data = speed_collapsed, family = gaussian(link = "identity")) %>% summary()

speed_collapsed %>% filter(!is.na(Species)) %>%  ggplot(aes(x=Zooid_length_mm, y=Speed_mms_abs)) + geom_point(aes(col=Architecture))+
  theme_bw() + ylab("Speed (mm/s)") + xlab("Zooid length (mm)") + geom_abline() +ylim(0,150)

speed_collapsed %>% filter(!is.na(Species)) %>%  ggplot(aes(x=Zooid_length_mm, y=Speed_mms_abs/Pulses_per_second)) + geom_point(aes(col=Architecture))+
  theme_bw() + ylab("Speed (mm/pulsation)") + xlab("Zooid length (mm)") + geom_abline() +ylim(0,75)

write.csv(speed_collapsed, "EMSpeeds_final_annotated.tsv")

# Create the ggplot with different ablines for each Architecture
ggplot(data = speed_collapsed , aes(x = Zooid_length_mm, y = Speed_mms_abs)) +
  geom_point(data=speed_collapsed %>% filter(Architecture=="Linear"),col="green", na.rm = TRUE) +  # Add na.rm = TRUE to remove missing values
  geom_line(data=speed_collapsed %>% filter(Architecture=="Linear"), col="green")  +
  geom_point(data=speed_collapsed %>% filter(Architecture=="Transversal"),col="blue", na.rm = TRUE) +
  geom_line(data=speed_collapsed %>% filter(Architecture=="Transversal"), col="blue")  +
  geom_point(data=speed_collapsed %>% filter(Architecture=="Whorl"),col="purple", na.rm = TRUE) +
  geom_line(data=speed_collapsed %>% filter(Architecture=="Whorl"), col="purple")  +
  geom_point(data=speed_collapsed %>% filter(Architecture=="Cluster"),col="brown", na.rm = TRUE) +
  geom_line(data=speed_collapsed %>% filter(Architecture=="Cluster"), col="brown")  +
  geom_point(data=speed_collapsed %>% filter(Architecture=="Bipinnate"),col="pink", na.rm = TRUE) +
  geom_line(data=speed_collapsed %>% filter(Architecture=="Bipinnate"), col="pink")  +
  ylim(0,150)+
  theme_bw()

# Create an empty list to store model fits
model_list <- list()

# List of Architectures
architectures <- c("Linear", "Transversal", "Cluster", "Whorl", "Bipinnate")

# Fit GLM models for each Architecture
for (arch in architectures) {
  model <- glm(Speed_mms_abs ~ Pulses_per_second, data = speed_collapsed %>% filter(Architecture == arch), family = gaussian(link = "identity"))
  model_list[[arch]] <- model
}

# Create the ggplot
ggplot(data = speed_collapsed %>% filter(!is.na(Species)), aes(x = Pulses_per_second, y = Speed_mms_abs)) +
  geom_point(aes(col = Architecture), na.rm = TRUE) +  # Add na.rm = TRUE to remove missing values
  lapply(architectures, function(arch) {
    geom_smooth(
      method = "glm",
      formula = y ~ x,
      data = speed_collapsed %>% filter(Architecture == arch),
      col=ifelse(arch=="Linear", "green", ifelse(arch=="Transversal", "blue", ifelse(arch=="Whorl", "purple", ifelse(arch=="Cluster", "brown", ifelse(arch=="Bipinnate", "magenta", "grey")))))
    )
  }) +
  ylim(0,150)+
  theme_bw()

# Create an empty list to store model fits
model_list <- list()

# List of Architectures
species <- unique(speed_collapsed$Species)

# Fit GLM models for each Species
for (sp in species) {
  model <- glm(Speed_mms_abs ~ Pulses_per_second, data = speed_collapsed %>% filter(Species == sp), family = gaussian(link = "identity"))
  model_list[[sp]] <- model
}

# Create the ggplot
ggplot(data = speed_collapsed %>% filter(!is.na(Species)), aes(x = Pulses_per_second, y = Speed_mms_abs)) +
  geom_point(aes(col = Species), na.rm = TRUE) +  # Add na.rm = TRUE to remove missing values
  lapply(species, function(sp) {
    geom_smooth(
      method = "glm",
      formula = y ~ x,
      data = speed_collapsed %>% filter(Species == sp))
  }) +
  ylim(0,150)+
  theme_bw()


speed_collapsed %>% filter(!is.na(Species)) %>%  ggplot(aes(x=Zooid_length_mm, y=Speed_mms_abs)) + geom_point(aes(col=Architecture))+
  theme_bw() + geom_abline()

glm(Speed_mms_abs ~ Pulses_per_second+Zooid_length_mm+Architecture, data = speed_collapsed, family = gaussian(link = "identity")) %>% summary()

  #ABP

speed_collapsed %>% filter(!is.na(Species)) %>%  ggplot(aes(Species, BLperPulse)) + geom_boxplot(aes(fill=Architecture))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90)) + ylab("Speed (zooid_lengths/Pulsation)")

speed_collapsed %>% filter(!is.na(Species)) %>% group_by(Architecture) %>% ggplot(aes(Zooid.number, BLperPulse)) + geom_point(aes(shape=Species, col=Species))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90))  + facet_wrap(~Architecture, scales="free") +  scale_color_manual(values = rep(c("red", "blue"), each = 14 / 2)) +
  scale_shape_manual(values = 1:14) 

speed_collapsed %>% filter(!is.na(Species) & Architecture != "Whorl chain") %>% ggplot(aes(Zooid.number, BLperPulse)) + geom_point(aes(col=Architecture))+
  theme_bw()

speed_collapsed %>% filter(!is.na(Species) & Architecture != "Whorl chain") %>% group_by(Architecture) %>% 
  ggplot(aes(Zooid.number, BLperPulse)) + geom_point()+
  theme_bw()+facet_wrap(~Architecture, scales="free")

speed_collapsed %>%
  filter(!is.na(Species) & Architecture != "Whorl chain") %>% mutate(CSAmode = case_when(
    Architecture %in% c("Transversal", "Cluster", "Whorl") ~ "CSA_Scaling",
    Architecture %in% c("Linear", "Bipinnate") ~ "CSA_Constant",
    TRUE ~ as.character(Architecture)
  )) %>% 
  group_by(CSAmode) %>%
  ggplot(aes(Zooid.number, BLperPulse)) +
  geom_point(aes(col=Architecture)) +
  scale_color_manual(values=c("red","blue","green","brown","purple"))+
  geom_smooth(method = 'lm', se = TRUE) +  # Add regression line
  theme_bw() +
  ylim(c(0,NA))+
  ylab("Speed (Body lengths per pulse)")+
  xlab("Number of zooids")+
  facet_wrap(~CSAmode, scales = "free_x")

speed_collapsed %>%
  filter(!is.na(Species) & Architecture != "Whorl chain") %>% mutate(CSAmode = case_when(
    Architecture %in% c("Transversal", "Cluster", "Whorl") ~ "CSA_Scaling",
    Architecture %in% c("Linear", "Bipinnate") ~ "CSA_Constant",
    TRUE ~ as.character(Architecture)
  )) -> speed_groupedArch


### GLM mega-analysis

glm(Speed_mms_abs ~ Zooid_length_mm + Pulses_per_second + Zooid.number + CSAmode, data=speed_groupedArch) %>% summary() #significant
glm(Speed_mms_abs ~ Zooid_length_mm + Pulses_per_second + Zooid.number * CSAmode, data=speed_groupedArch) %>% summary()
glm(Speed_mms_abs ~ Zooid_length_mm + Zooid.number:CSAmode, data=speed_groupedArch) %>% summary()
glm(Speed_mms_abs ~ Pulses_per_second + Zooid.number:CSAmode, data=speed_groupedArch) %>% summary()

library(caret)

speed_collapsed %>%
  filter(!is.na(Species) & Architecture != "Whorl chain") -> speed_filtered

for(a in unique(speed_filtered$Architecture)){
  speed_a <- filter(speed_filtered, Architecture == a)
  print(a)
  glm(BLperPulse ~ Zooid.number, data = speed_a) %>% summary() %>% .$coefficients %>% .[1,4] %>% print()
}

mean_BL <- speed_collapsed %>% group_by(Species) %>% mutate(mean(BLperPulse, na.rm=T))
names(mean_BL)[11] <- "meanBLperPulse"
mean_BL %>% filter(!is.na(Species) & Architecture != "Whorl chain" & !is.nan(BLperPulse)) %>% group_by(Architecture) %>% 
    ggplot(aes(Species, BLperPulse/Zooid.number)) + 
    geom_boxplot(aes(fill=meanBLperPulse), na.rm=T)+
    theme_bw()+
    ylab("Speed (Body lengths per pulse) / Zooid number")+
    theme(axis.text.x = element_text(angle = 90))+
    facet_wrap(~Architecture, scales="free_x", ncol=5)

