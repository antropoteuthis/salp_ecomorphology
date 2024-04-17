# BACKBONE CODE FOR ALL SPEEDS #

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
library(purrr)
library(ggsignif)
library(ggpubr)

## Acquire list of filenames with speed data in your folder ##

folder_path <- "~/salp_ecomorphology/EventMeasure_TXT_FINAL/"
list_csv_files <- list.files(path = folder_path, pattern = "\\.TXT$", full.names = T)

#change path to the folder where you store them
df_list <- map_dfr(list_csv_files, function(file_path) {
  # Read the CSV file skipping the first 4 rows
  df <- read.csv(file_path, sep = '\t', header = TRUE, stringsAsFactors = FALSE, skip = 4,fill = TRUE, check.names = FALSE)[,c(1:3,6:8,12:13,31)]
  df <- filter(df, `RMS (mm)`< 2)
  df <- df[,c(-7,-8)]
  df$Number <- as.character(df$Number)
  # Extract the file name without extension
  file_name <- tools::file_path_sans_ext(basename(file_path))

  # Add a new column "Filename" with the file name
  df <- mutate(df, Filename = file_name)

  return(df)
})

# added an RMS cutoff of 750

## Identify the junk particle rows by the zooid number zero AND chop off the columns that are not useful ##

df_list %>% filter(Number %in% c("0","1")) %>% 
  mutate(is_junk = (Number=="0") == FALSE) -> keyvars #is_junk identifies the "junk" reference particle

## Make a tidy data frame by splitting the data into zooid and junk rows and then full-joining them side by side ##

tidy <- full_join(group_split(keyvars, is_junk)[[1]][,1:7], 
                  group_split(keyvars, is_junk)[[2]][,1:7], 
                  by=c("Filename","Frame","Time (mins)"), 
                  suffix= c("Zooid","Junk")) %>% as.data.frame() %>% .[,c(-7,-11)]

## Define a function to calculate time lags between rows to get the delta of time between measures, 
## get the vectors between the zooid and the junk particle, 
## calculate XYZ components of the distance, 
## use pythagoras to estimate the hypothenuse and divide by the delta of time

wrangle <- function(df){
  df %>% arrange(Frame) %>%
    mutate(Delta_time_s = (`Time (mins)` - lag(`Time (mins)`, default = first(`Time (mins)`)))*60) %>% 
    mutate(Xd_mm = `X (mm)Zooid` - `X (mm)Junk`, Yd_mm = `Y (mm)Zooid` - `Y (mm)Junk`, Zd_mm =  `Z (mm)Zooid` - `Z (mm)Junk`) %>% 
    mutate(Delta_X_mm = Xd_mm - lag(Xd_mm, default = first(Xd_mm)), 
           Delta_Y_mm = Yd_mm - lag(Yd_mm, default = first(Yd_mm)),  
           Delta_Z_mm = Zd_mm - lag(Zd_mm, default = first(Zd_mm))) %>% 
    mutate(Distance_mm = sqrt(Delta_X_mm^2 + Delta_Y_mm^2 + Delta_Z_mm^2)) %>% 
    mutate(Speed_mm_s = Distance_mm/Delta_time_s) %>% filter(Speed_mm_s>0) %>% return()
}

## Run the 3D wrangle function for the first (frontal) zooid ##

speed_data <- group_split(tidy, Filename) %>% lapply(wrangle) %>% do.call(rbind, .) %>% as.data.frame()

## 2D data ##

list_2D_files <- list.files(path = "~/salp_ecomorphology/Speeds2D_csv_FINAL/")
setwd("~/salp_ecomorphology/Speeds2D_csv_FINAL/")
df2d = do.call(dplyr::bind_rows, lapply(list_2D_files, function(x){
  tb <- read.csv(x, sep=',', header = T, stringsAsFactors = F)
  tb$Filename <- basename(x)
  tb <- dplyr::mutate(tb, is_junk = (X.1 %% 2))
  return(tb)
          }))
keyvars2D <-df2d[order(df2d$Slice),c(6,7,11:13)]
 #is_junk identifies the "junk" reference particle

tidy2D <- full_join(dplyr::group_split(keyvars2D, is_junk)[[1]][,1:4], 
                    dplyr::group_split(keyvars2D, is_junk)[[2]][,1:4], 
                  by=c("Filename","Slice"), 
                  suffix= c("Zooid","Junk")) %>% 
  as.data.frame()

names(tidy2D)[3]<-"Frame"
tidy2D <- tidy2D[order(tidy2D$Filename),c(4,3,1,2,5,6)]
#fps = 30

inventory2D <- read.csv("~/salp_ecomorphology/Inventory_FINAL.csv") %>% 
  dplyr::filter(!is.na(FPS)) %>% 
  .[,c("Video.name..or.master..or.left.","Species","Architecture","FPS","Zooid.number")]

names(inventory2D)[1] <- "Filename"
tidy2D$Filename <- str_replace_all(tidy2D$Filename, ".csv", "")
inventory2D$Filename <- str_replace_all(inventory2D$Filename, ".mp4", "")
dplyr::left_join(tidy2D, inventory2D[,c(1,4)], by = "Filename") -> tidy2D

wrangle2D <- function(df){
    df %>% mutate(Delta_time_s = (Frame - lag(Frame, default = first(Frame)))/FPS) %>% 
    mutate(Xd_mm = XZooid - XJunk, Yd_mm = YZooid - YJunk) %>% 
    mutate(Delta_X_mm = Xd_mm - lag(Xd_mm, default = first(Xd_mm)), 
           Delta_Y_mm = Yd_mm - lag(Yd_mm, default = first(Yd_mm))) %>% 
    mutate(Distance_mm = sqrt(Delta_X_mm^2 + Delta_Y_mm^2)) %>% 
    mutate(Speed_mm_s = Distance_mm/Delta_time_s) %>% return()
}

speed_data2d <- group_split(tidy2D, Filename) %>% lapply(wrangle2D) %>% do.call(rbind, .) %>% as.data.frame() 
left_join(speed_data2d, inventory2D[,c(1:3,5)], by = "Filename") %>% filter(!is.na(Speed_mm_s)) -> speed_data2d


## Get video inventory metadata and match it by file name to the speeds to obtain: Species, Architecture, & Number of zooids ##

inventory <- read.csv("~/salp_ecomorphology/Inventory_FINAL.csv", header = T, stringsAsFactors = F) #change path to where you saved it
names(inventory)[3] <- "Filename"
inventory$Species %>% str_replace_all("Iasis cylindrica \\(yellowtail\\)","Iasis cylindrica") -> inventory$Species
inventory$Species %>% str_replace_all("Metcalfina hexagona\\?","Metcalfina hexagona") -> inventory$Species
inventory$Species %>% str_replace_all("Soestia zonaria\\?","Soestia zonaria") -> inventory$Species
inventory$Species %>% str_replace_all("cf. ","") -> inventory$Species
inventory_clean <- inventory[,c(3,5,6,7,15)] #Keeping columns: filename, species, architecture, number of zooids,FPS

## Full join speeds with inventory metadata ##

speed_annotated <- full_join(speed_data, inventory_clean, by = "Filename" ) %>% 
  filter(Speed_mm_s != Inf) %>% filter(Species != "NA") %>%
  mutate(Speed_mms_abs = abs(Speed_mm_s)) #transform negative vectorial speeds into absolute scalars

## Load morphometrics data ##

mm <- read.csv("~/salp_ecomorphology/Morphometrics_FINAL.csv", header = T, stringsAsFactors = F) #change path to your file
mm <- mm[,c(1,11:14,17:19)]
names(mm)[1] <- "Filename"
coalesce_by_column <- function(df) {
  return(dplyr::coalesce(!!! as.list(df)))
}
mm %>% group_by(Filename) %>% summarise_all(coalesce_by_column) -> mm
mm <- filter(mm, Filename != "")

##Merge 2D and 3D data

names(speed_data2d)[3:6] <- c("X..mm.Zooid", "Y..mm.Zooid", "X..mm.Junk", "Y..mm.Junk")
mutate(speed_data2d, Speed_mms_abs = abs(Speed_mm_s))->speed_data2d
bind_rows(speed_annotated, speed_data2d) -> speed_annotated

## Full-join morphometrics annotations to the speeds ##

speed_annotated <- full_join(speed_annotated, mm, by = "Filename")

## Load pulsation rates ##

pulsation <- read.csv("~/salp_ecomorphology/Pulsation_FINAL.csv", header = T, stringsAsFactors = F) #change path to your file
pulsation <- pulsation[,c(1,16)]
names(pulsation) <- c("Filename", "Pulses_per_second")
pulsation$Filename %>% str_remove_all(".mp4") -> pulsation$Filename

## Full join pulsation rates to the speeds ##

speed_annotated <- full_join(speed_annotated, pulsation, by = "Filename" )

## Estimate normalized speed in body_lengths per pulse ##

speed_annotated <- mutate(speed_annotated, BLperSecond = Speed_mms_abs/Zooid_length_mm )
speed_annotated <- mutate(speed_annotated, BLperPulse = BLperSecond/Pulses_per_second )

## Collapse into the mean per specimen

speed_annotated[,c(1,18:21,23,29:37)] %>% unique() -> speed_collapsed
speed_collapsed %>% group_by(Filename, Species, Architecture) %>% 
  summarise(across(.cols = c(Speed_mm_s, Zooid.number, Speed_mms_abs, Zooid_length_mm, Pulses_per_second, BLperSecond, BLperPulse), 
                   ~ mean(., na.rm = TRUE))) %>% filter(!is.na(Species)) %>% as.data.frame() -> speed_collapsed

## Plot 


  #KAW

### FIG 2A ###

architecture_order <- c("Transversal", "Oblique", "Linear", "Bipinnate", "Helical", "Whorl", "Cluster")

F2A <- speed_annotated %>% filter(!is.na(Species) & Architecture != "Whorl chain") %>%  
  ggplot(aes(x = factor(Species %>% as.character(), levels = unique(Species[order(factor(Architecture, levels = architecture_order))])), 
             y=Speed_mms_abs)) + geom_boxplot(aes(fill=Architecture))+guides(fill = "none") +
  scale_fill_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                               "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                               "Cluster" = "magenta")) +
  ylab("Speed (mm/s)") + xlab("Species")+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=1))

### FIG 2B ###

F2B <- speed_annotated %>% filter(!is.na(Species) & !is.na(BLperPulse) & Architecture != "Whorl chain") %>%  
  ggplot(aes(x = factor(Species %>% as.character(), levels = unique(Species[order(factor(Architecture, levels = architecture_order))])), 
             y=BLperPulse)) + geom_boxplot(aes(fill=Architecture))+guides(fill = "none") +
  scale_fill_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                               "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                               "Cluster" = "magenta")) +
  theme(axis.text.x = element_text(angle = 90)) + ylab("Speed (zooids/pulse)") +
  xlab("Species") +theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust=1))


###############

### FIG 2C ###

F2C <- speed_annotated %>% filter(!is.na(Architecture) & Architecture != "Whorl chain") %>% 
  ggplot(aes(x=Architecture, y=Speed_mms_abs)) + 
  geom_boxplot(aes(fill=Architecture))+
  scale_fill_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                               "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                               "Cluster" = "magenta")) + theme_bw()+
 theme(axis.text.x = element_text(angle = 90, hjust=1)) + ylab("Speed (mm/s)")+ guides(fill = "none") +
  geom_signif(comparisons = list(c("Bipinnate","Cluster"), c("Linear","Bipinnate"), c("Linear","Helical")), map_signif_level=TRUE, test="t.test")

### FIG 3B ###

F2D <- speed_annotated %>% filter(!is.na(Architecture) & BLperPulse<5 & Architecture != "Whorl chain") %>%  ggplot(aes(x=Architecture, y=BLperPulse)) + 
  geom_boxplot(aes(fill=Architecture))+
  scale_fill_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                               "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                               "Cluster" = "magenta")) +theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1)) + ylab("Speed (zooids/pulse)") +
  geom_signif(comparisons = list(c("Bipinnate","Cluster"), c("Linear","Bipinnate"), c("Linear","Cluster"), c("Linear","Helical")), map_signif_level=TRUE, test="t.test")

wrap_plots(F2A, F2B, F2C, F2D)

speed_annotated %>%
  filter(Architecture %in% c("Linear", "Cluster")) %>%
  group_by(Architecture) %>%
  t.test(BLperPulse ~ Architecture, data = .) %>% list() %>% .[[1]] %>% print()

speed_annotated %>%
  filter(Architecture %in% c("Linear", "Bipinnate")) %>%
  group_by(Architecture) %>%
  t.test(BLperPulse ~ Architecture, data = .) %>% list() %>% .[[1]] %>% print()

speed_annotated %>%
  filter(Architecture %in% c("Helical", "Bipinnate")) %>%
  group_by(Architecture) %>%
  t.test(BLperPulse ~ Architecture, data = .) %>% list() %>% .[[1]] %>% print()

speed_annotated %>%
  filter(Architecture %in% c("Transversal", "Whorl")) %>%
  group_by(Architecture) %>%
  t.test(BLperPulse ~ Architecture, data = .) %>% list() %>% .[[1]] %>% print()

speed_annotated %>%
  filter(Architecture %in% c("Helical", "Bipinnate")) %>%
  group_by(Architecture) %>%
  t.test(Speed_mm_s ~ Architecture, data = .) %>% list() %>% .[[1]] %>% print()

##############

#### SM Figure 1A #####

Sm1A <- speed_annotated %>%
  filter(!is.na(Species) & Architecture != "Whorl chain") %>%
  ggplot(aes(x = Pulses_per_second, y = Speed_mms_abs)) +
  stat_summary(fun = mean, geom = "point", aes(col = Architecture), size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(col = Architecture), width = 0.2) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  theme_bw() +
  ylab("Speed (mm/s)") +
  xlab("Pulsation rate (pulses/s)") +
  guides(color = "none")


#### SM Figure 1B #####

Sm1B <- speed_annotated %>%
  filter(!is.na(Species) & Architecture != "Whorl chain") %>%
  ggplot(aes(x = Pulses_per_second, y = BLperSecond)) +
  stat_summary(fun = mean, geom = "point", aes(col = Architecture), size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(col = Architecture), width = 0.2) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  theme_bw() +
  ylab("Speed (zooids/s)") +
  xlab("Pulsation rate (pulses/s)")


wrap_plots(Sm1A, Sm1B)

#######################

# GLMs of speed vs pulsation rate

lm(Speed_mms_abs ~ Pulses_per_second, data = speed_annotated, family = gaussian(link = "identity")) %>% summary()
lm(BLperSecond ~ Pulses_per_second, data = speed_annotated, family = gaussian(link = "identity")) %>% summary()

#### SM Figure 2A #####

speed_annotated %>% filter(!is.na(Species) & Architecture != "Whorl chain") %>%  
  ggplot(aes(x=Zooid_length_mm, y=Speed_mms_abs)) + geom_point(aes(col=Architecture)) + geom_smooth(method="lm",color="white") +
  scale_color_manual(values = c("cyan4", "red1", "darkorange1", "magenta", "gold1", "darkorchid4", "green4")[c(1, 4, 3, 5, 2, 7, 6)]) +
  theme_bw() + ylab("Speed (mm/s)") + xlab("Zooid length (mm)") + guides(color = "none") + ylim(0,200)+theme_dark()+
  theme( text = element_text(color = "white"), axis.text = element_text(color = "white"),  # Set text color to white
         panel.background = element_rect(fill = "black"))

Sm2A <- speed_annotated %>%
  filter(!is.na(Species) & Architecture != "Whorl chain") %>%
  ggplot(aes(x = Zooid_length_mm, y = Speed_mms_abs)) +
  stat_summary(fun = mean, geom = "point", aes(col = Architecture), size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(col = Architecture), width = 0.2) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  theme_bw() + ylab("Speed (mm/s)") + xlab("Zooid length (mm)") + guides(color = "none")

#### SM Figure 2B #####

speed_annotated %>% filter(!is.na(Species) & Architecture != "Whorl chain") %>%  
  ggplot(aes(x=Zooid_length_mm, y=Speed_mms_abs/Pulses_per_second)) + geom_point(aes(col=Architecture)) + geom_smooth(method="lm",color="white") +
  scale_color_manual(values = c("cyan4", "red1", "darkorange1", "magenta", "gold1", "darkorchid4", "green4")[c(1, 4, 3, 5, 2, 7, 6)]) +
  theme_bw() + ylab("Speed (mm/pulse)") + xlab("Zooid length (mm)") +theme_dark()+
  theme( text = element_text(color = "white"), axis.text = element_text(color = "white"),  # Set text color to white
         panel.background = element_rect(fill = "black"))

Sm2B <- speed_annotated %>%
  filter(!is.na(Species) & Architecture != "Whorl chain") %>%
  ggplot(aes(x = Zooid_length_mm, y=Speed_mms_abs/Pulses_per_second)) +
  stat_summary(fun = mean, geom = "point", aes(col = Architecture), size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(col = Architecture), width = 0.2) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  theme_bw() + ylab("Speed (mm/pulse)") + xlab("Zooid length (mm)")
  
wrap_plots(Sm2A, Sm2B)

###########################

# GLMs of speed vs zooid length

glm(Speed_mms_abs ~ Zooid_length_mm, data = speed_annotated, family = gaussian(link = "identity")) %>% summary()
glm(Speed_mms_abs/Pulses_per_second ~ Zooid_length_mm, data = speed_annotated, family = gaussian(link = "identity")) %>% summary()

### EXPORT ALL SPEEDS ####

write.csv(speed_collapsed, "~/salp_ecomorphology/EMSpeeds_final_annotated.tsv")

###########################

glm(Speed_mms_abs ~ Pulses_per_second+Zooid_length_mm+Architecture, data = speed_annotated, family = gaussian(link = "identity")) %>% summary()

#### FIGURE 4 ######

speed_annotated %>%
  filter(!is.na(Species) & Architecture != "Whorl chain" & !is.na(BLperPulse)) %>% mutate(CSAmode = case_when(
    Architecture %in% c("Transversal", "Cluster", "Whorl", "Oblique") ~ "Scaling Frontal Area",
    Architecture %in% c("Linear", "Bipinnate", "Helical") ~ "Constant Frontal Area",
    TRUE ~ as.character(Architecture)
  )) %>% 
  group_by(CSAmode) %>% 
  ggplot(aes(Zooid.number, BLperPulse)) +
  stat_summary(fun = mean, geom = "point", aes(col = Architecture), size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(col = Architecture), width = 0.2) +
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  geom_smooth(col = "black", method = 'lm', se = TRUE) +  # Add regression line
  theme_bw() +
  ylab("Speed (zooids/pulse)")+
  xlab("Number of zooids")+
  facet_wrap(~CSAmode, scales = "free_x")

glm(BLperPulse ~ Zooid.number, data=speed_annotated %>% filter(Architecture %in% c("Linear","Bipinnate", "Helical"))) %>% summary()
glm(BLperPulse ~ Zooid.number, data=speed_annotated %>% filter(!(Architecture %in% c("Linear","Bipinnate", "Helical")))) %>% summary()

#########

speed_annotated %>%
  filter(!is.na(Species) & Architecture != "Whorl chain") %>% mutate(CSAmode = case_when(
    Architecture %in% c("Transversal", "Cluster", "Whorl", "Oblique") ~ "CSA_Scaling",
    Architecture %in% c("Linear", "Bipinnate","Helical") ~ "CSA_Constant",
    TRUE ~ as.character(Architecture)
  )) -> speed_groupedArch


### GLM mega-analysis

# glm(Speed_mms_abs ~ Zooid_length_mm + Pulses_per_second + Zooid.number + CSAmode, data=speed_groupedArch) %>% summary() #significant
# glm(Speed_mms_abs ~ Zooid_length_mm + Pulses_per_second + Zooid.number * CSAmode, data=speed_groupedArch) %>% summary()
# glm(Speed_mms_abs ~ Zooid_length_mm + Zooid.number:CSAmode, data=speed_groupedArch) %>% summary()
# glm(Speed_mms_abs ~ Pulses_per_second + Zooid.number:CSAmode, data=speed_groupedArch) %>% summary()

# glm(Speed_mms_abs ~ Zooid_length_mm + Pulses_per_second + Zooid.number:CSAmode, data=speed_groupedArch) %>% summary()
# glm(Speed_mms_abs ~ Zooid_length_mm + Pulses_per_second + Zooid.number + CSAmode + Zooid.number:CSAmode, data=speed_groupedArch) %>% summary()
glm(Speed_mms_abs ~ Zooid_length_mm + Pulses_per_second + Zooid.number + CSAmode, data=speed_groupedArch) %>% summary()
glm(Speed_mms_abs ~ Zooid_length_mm + Pulses_per_second + Zooid.number + CSAmode, data=speed_groupedArch) %>% anova()

# Assuming 'model' is your GLM model
null_deviance <- 1213636

# Deviance reduction by each factor
deviance_reduction <- c(70229, 42748, 9859, 323369) #edit with anova values

# Calculate the percentage of deviance explained for each factor
percentage_deviance_explained <- deviance_reduction / null_deviance * 100

# Display the results
for (i in seq_along(percentage_deviance_explained)) {
  cat(paste("Percentage of Deviance Explained by factor", i, ":", percentage_deviance_explained[i], "%\n"))
}

#glm(BLperPulse ~ Zooid.number:CSAmode, data=speed_groupedArch) %>% summary()

### GLM by architecture

glm(Speed_mms_abs ~ Zooid_length_mm + Pulses_per_second + Zooid.number, data=speed_annotated %>% filter(Architecture=="Linear")) %>% summary()
glm(Speed_mms_abs ~ Zooid_length_mm + Pulses_per_second + Zooid.number, data=speed_annotated %>% filter(Architecture=="Cluster")) %>% summary()
glm(Speed_mms_abs ~ Zooid_length_mm + Pulses_per_second + Zooid.number, data=speed_annotated %>% filter(Architecture=="Bipinnate")) %>% summary()
glm(Speed_mms_abs ~ Pulses_per_second + Zooid.number, data=speed_annotated %>% filter(Architecture=="Whorl")) %>% summary()
glm(Speed_mms_abs ~ Pulses_per_second + Zooid.number, data=speed_annotated %>% filter(Architecture=="Transversal")) %>% summary()

#### FIGURE S4 ####

speed_annotated %>%
  filter(!is.na(Species) & Architecture != "Whorl chain" & Architecture != "Oblique" &  Architecture != "Helical" & !is.na(BLperPulse)) %>%
  group_by(Architecture) %>%
  ggplot(aes(Zooid.number, BLperPulse)) +
  stat_summary(fun = mean, geom = "point", aes(col = Architecture), size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(col = Architecture), width = 0.2) +
  geom_smooth(aes(color = Architecture), method = 'lm', se = TRUE) +  # Match color to Architecture
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  theme_bw() +
  ylim(c(0, NA)) +
  facet_wrap(~Architecture, scales = "free") +
  ylab("Speed (zooids/pulse)") +xlab("Number of zooids") + guides(color = "none") 

lm(BLperPulse ~ Zooid.number, data=speed_annotated %>% filter(Architecture=="Whorl")) %>% summary()
lm(BLperPulse ~ Zooid.number, data=speed_annotated %>% filter(Architecture=="Cluster")) %>% summary()
lm(BLperPulse ~ Zooid.number, data=speed_annotated %>% filter(Architecture=="Transversal")) %>% summary()
lm(BLperPulse ~ Zooid.number, data=speed_annotated %>% filter(Architecture=="Linear")) %>% summary()
lm(BLperPulse ~ Zooid.number, data=speed_annotated %>% filter(Architecture=="Bipinnate")) %>% summary()

#####
### SM FIGURE 3 ###
# BLperPulse vs Zooid # by Species
speed_annotated %>%
  filter(!is.na(Species) & Architecture != "Whorl chain" & !is.na(BLperPulse)) %>%
  group_by(Species) %>%
  ggplot(aes(Zooid.number, BLperPulse))+
  stat_summary(fun = mean, geom = "point", aes(col = Architecture), size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(col = Architecture), width = 0.2) +
  geom_smooth(aes(color = Architecture), method = 'lm', se = TRUE) +  # Match color to Architecture
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  theme_bw() +
  ylim(c(0, NA)) +
  facet_wrap(~Species, scales = "free") +
  ylab("Speed (zooids/pulse)") + xlab("Number of zooids") + guides(color = "none") 

#DVZSangle

dvsz <- read.csv("~/salp_ecomorphology/DV_Zooid.stolon.angle.tsv",sep='\t', stringsAsFactors = F)[,c(5,8,7,10)]
names(dvsz)[4] <- "Dorsoventral_Zooid_Stolon_Angle"
dvsz$Species[which(dvsz$Species == "Thalia longicauda" | dvsz$Species == "Thalia cicar")] <- "Thalia sp."

meanDVZS <- dvsz[,-3] %>%
  group_by(Species, Architecture) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>% as.data.frame() #Get the mean DVZS angle for each species, since the specimens dont match

speed_dvzs <- full_join(speed_annotated, meanDVZS, by = c("Species","Architecture")) %>% filter(Architecture!= "Whorl chain") #join with the speeds

speed_dvzs %>% ggplot(aes(x=Dorsoventral_Zooid_Stolon_Angle, y=BLperPulse)) +
  stat_summary(fun = mean, geom = "point", aes(col = Architecture), size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(col = Architecture), width = 0.2) +
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +  # Set color palette for each Architecture
  geom_smooth(method = 'lm', se = TRUE, color="black")+theme_bw()+xlab("Dorsoventral Zooid Rotation Angle")+ylab("Speed (body lengths/pulsation)")

#### FIGURE 3 #####

   ## means with crossbars relative speed #
DVR <- speed_dvzs %>%
  filter(!is.na(Speed_mms_abs)) %>%
  group_by(Species, Architecture) %>%
  summarise(
    mean_Dorsoventral_Zooid_Stolon_Angle = mean(Dorsoventral_Zooid_Stolon_Angle, na.rm = TRUE),
    mean_BLperPulse = mean(BLperPulse, na.rm = TRUE),
    sd_BLperPulse = sd(BLperPulse, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = mean_Dorsoventral_Zooid_Stolon_Angle, y = mean_BLperPulse, col = Architecture)) +
  geom_smooth(method = 'lm', se = TRUE, color = "black") +
  geom_point(cex=3) +
  geom_errorbar(aes(ymin = mean_BLperPulse - sd_BLperPulse, ymax = mean_BLperPulse + sd_BLperPulse),
    width = 2.5,
    position = position_dodge(0.2)
  ) +
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  theme_bw() +
  xlab("Dorsoventral Zooid Rotation Angle") +
  ylab("Speed (body lengths/pulsation)") 

## means with crossbars absolute speed #

DVA <- speed_dvzs %>%
  filter(!is.na(Speed_mms_abs)) %>%
  group_by(Species, Architecture) %>%
  summarise(
    mean_Dorsoventral_Zooid_Stolon_Angle = mean(Dorsoventral_Zooid_Stolon_Angle, na.rm = TRUE),
    mean_Speed_mm_s = mean(Speed_mms_abs, na.rm = TRUE),
    sd_Speed_mm_s = sd(Speed_mms_abs, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = mean_Dorsoventral_Zooid_Stolon_Angle, y = mean_Speed_mm_s, col = Architecture)) +
  geom_smooth(method = 'lm', se = TRUE, color = "black") +
  geom_point(cex=3) +
  geom_errorbar(aes(ymin = mean_Speed_mm_s - sd_Speed_mm_s, ymax = mean_Speed_mm_s + sd_Speed_mm_s),width = 2.5, position = position_dodge(0.2)) +
  scale_color_manual(values = c("Transversal" = "green4", "Oblique" = "red1", "Linear" = "darkorange1", 
                                "Bipinnate" = "cyan4", "Helical" = "gold1", "Whorl" = "darkorchid4", 
                                "Cluster" = "magenta")) +
  theme_bw() +
  xlab("Dorsoventral Zooid Rotation Angle") +
  ylab("Speed (mm/s)")+ guides(color="none")

wrap_plots(DVA,DVR)

lm(Speed_mms_abs ~ Dorsoventral_Zooid_Stolon_Angle , data = speed_dvzs) %>% summary()
lm(BLperPulse ~ Dorsoventral_Zooid_Stolon_Angle , data = speed_dvzs) %>% summary()

############## ## ###

### JET motion

speed_annotated %>% filter(Architecture != "Whorl chain" & !is.na(Architecture) ) %>% 
  glm(Speed_mms_abs ~ Jet_motion_angle, data=.) %>% summary()

speed_annotated %>% filter(Architecture != "Whorl chain" & !is.na(Architecture) ) %>% 
  ggplot(aes(x=Jet_motion_angle, y=Speed_mm_s)) +
  geom_point(aes(col=Architecture)) +
  geom_smooth(method = 'lm', se = TRUE)+ylim(0,250)

speed_annotated %>% filter(Architecture != "Whorl chain" & !is.na(Architecture) ) %>% 
  glm(BLperPulse ~ Jet_motion_angle, data=.) %>% summary()

speed_annotated %>% filter(Architecture != "Whorl chain" & !is.na(Architecture) ) %>% 
  ggplot(aes(x=Jet_motion_angle, y=BLperSecond)) +
  geom_point(aes(col=Architecture)) +
  geom_smooth(method = 'lm', se = TRUE)+ylim(0,25)

speed_annotated %>% filter(Architecture != "Whorl chain" & !is.na(Architecture) ) %>% 
  ggplot(aes(x=Jet_motion_angle, y=BLperPulse)) +
  geom_point(aes(col=Architecture)) +
  geom_smooth(method = 'lm', se = TRUE)+ylim(0,7)

## SUMMARY TABLES ##

speed_annotated %>% filter(!is.na(Species)) %>% 
  group_by(Species) %>%
  summarise(
    Nspecimens = n_distinct(Filename),
    Nmeasurements = n()
  ) 

speed_annotated %>% filter(Speed_mm_s > -90000) %>% 
  group_by(Filename, Species) %>%
  summarise(Nmeasurements = n()) %>%
  arrange(Species) %>% as.data.frame()

speed_annotated %>% filter(!is.na(Species)) %>% 
  group_by(Species) %>%
  summarise(
    meanSpeed = mean(Speed_mms_abs, na.rm=T)
  ) 

speed_annotated %>% filter(!is.na(Species)) %>% 
  group_by(Architecture) %>%
  summarise(
    meanSpeed = mean(Speed_mms_abs, na.rm=T)
  ) 

speed_annotated %>% filter(!is.na(Species)) %>% 
  group_by(Species) %>%
  summarise(
    meanSpeed = mean(BLperPulse, na.rm=T)
  ) 

speed_annotated %>% filter(!is.na(Species)) %>% 
  group_by(Architecture) %>%
  summarise(
    meanSpeed = mean(BLperPulse, na.rm=T)
  ) 


## SM Table 1

speed_annotated %>% filter(!is.na(Species)) %>% 
  group_by(Filename) %>%
  summarise(
    N = n(),
    Timespan_s = sum(Delta_time_s, na.rm = T)
  ) -> NandT

inventory[which(inventory_clean$Filename %in% speed_annotated$Filename),] -> sub_inv
sub_inv[,c(3,2,8,15,5,6,7)] -> sub_inv
full_join(sub_inv, mm[,c(1,4,5)], by = "Filename") -> morph_inv
full_join(morph_inv, pulsation, by = "Filename") -> mp_inv
full_join(mp_inv, speed_collapsed[,c(1,4)], by = "Filename") %>% filter(!is.na(Speed_mm_s) & Architecture != "Whorl chain") -> sp_inv
sp_inv %>% mutate(STEREO. = ifelse(substr(Filename, 1, 3) == "A00", "3D", "2D")) %>% .[,-3] -> sp_inv
full_join(sp_inv,NandT, by = "Filename") %>% .[,c(1:5,11,12,6:10)] -> sp_inv
names(sp_inv) <- c("Video file", "Type", "FPS", "Species", "Architecture", "Number of measurements", "Timespan of measurements (s)", "Number of zooids", "Mean zooid length (mm)", 
                   "Mean zooid width (mm)", "Pulsation rate (pulses/s)", 
                   "Mean swimming speed (mm/s)")
sp_inv %>% arrange(Architecture, Species) -> sp_inv
write_csv(sp_inv, "SMTable1_Specimens.csv")
