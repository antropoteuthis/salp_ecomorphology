library(tidyverse)
library(magrittr)
library(patchwork)
require(data.table)

##Import arguments
#myargs = commandArgs(trailingOnly=TRUE)
#file = myargs[1]
#raw <- read.csv(paste("~/salp_ecomorphology/", file, collapse=""), sep='\t', header = T, stringsAsFactors = F)

list_csv_files <- list.files(path = "~/Pictures/Brad_random/speed_data")
setwd("~/Pictures/Brad_random/speed_data")
df = do.call(rbind, lapply(list_csv_files, function(x) read.csv(x, sep=',', header = T, stringsAsFactors = F)))
df <-cbind(df[order(df$Slice),6:8], 1:(row_number(df)-1))
names(df)[4]<-"oddeven"
df %>% mutate(is_junk = (oddeven %% 2)) -> keyvars #is_junk identifies the "junk" reference particle

tidy <- full_join(group_split(keyvars, is_junk)[[1]][,1:3], 
                  group_split(keyvars, is_junk)[[2]][,1:3], 
                  by=c("Slice"), 
                  suffix= c("Zooid","Junk")) %>% 
  as.data.frame()

names(tidy)[3]<-"Frame"
fps = 30

wrangle <- function(df){
  df %>% 
    mutate(Delta_time_s = (Frame - lag(Frame, default = first(Frame)))/fps) %>% 
    mutate(Xd_mm = XZooid - XJunk, Yd_mm = YZooid - YJunk) %>% 
    mutate(Delta_X_mm = Xd_mm - lag(Xd_mm, default = first(Xd_mm)), 
           Delta_Y_mm = Yd_mm - lag(Yd_mm, default = first(Yd_mm))) %>% 
    mutate(Distance_mm = sqrt(Delta_X_mm^2 + Delta_Y_mm^2)) %>% 
    mutate(Speed_mm_s = 0.1*Distance_mm/Delta_time_s) %>% return()
}

speed_data <- wrangle(tidy)
