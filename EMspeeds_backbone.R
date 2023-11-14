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

speed_collapsed %>% filter(!is.na(Species)) %>%  ggplot(aes(Species, Speed_mms_abs)) + geom_boxplot(aes(fill=Architecture))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90))

speed_collapsed %>% filter(!is.na(Species)) %>%  ggplot(aes(Species, BLperPulse)) + geom_boxplot(aes(fill=Architecture))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90)) + ylab("Speed (zooid_lengths/Pulsation)")

speed_collapsed %>% filter(!is.na(Species)) %>% group_by(Species) %>% ggplot(aes(Zooid.number, BLperPulse)) + geom_point()+
  theme_bw()+theme(axis.text.x = element_text(angle = 90))  + facet_wrap(~Species, scales="free")
