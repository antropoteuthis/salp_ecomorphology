library(tidyverse)

##Import arguments
#myargs = commandArgs(trailingOnly=TRUE)
#file = myargs[1]
#raw <- read.csv(paste("~/salp_ecomorphology/", file, collapse=""), sep='\t', header = T, stringsAsFactors = F)

raw <- read.csv("~/salp_ecomorphology/A001C0167_20210709184129_0001.csv", sep=',', header = T, stringsAsFactors = F, skip=4)
raw[,c(1:3,6:8,31)] %>% mutate(is_junk = (1/Number<Inf) == FALSE) -> keyvars

tidy <- full_join(group_split(keyvars, is_junk)[[1]][,1:7], 
          group_split(keyvars, is_junk)[[2]][,1:7], 
          by=c("Filename","Frame","Time..mins."), 
          suffix= c("Zooid","Junk")) %>% 
          as.data.frame() %>% .[,c(-11)]

zooid_list <- group_split(tidy, NumberZooid)

wrangle <- function(df){
  df %>% arrange(Frame) %>%
    mutate(Delta_time_s = Time..mins. - lag(Time..mins., default = first(Time..mins.))) %>% 
    mutate(Xd_mm = X..mm.Zooid - X..mm.Junk, Yd_mm = Y..mm.Zooid - Y..mm.Junk, Zd_mm =  Z..mm.Zooid - Z..mm.Junk) %>% 
    mutate(Delta_X_mm = Xd_mm - lag(Xd_mm, default = first(Xd_mm)), 
           Delta_Y_mm = Yd_mm - lag(Yd_mm, default = first(Yd_mm)),  
           Delta_Z_mm = Zd_mm - lag(Zd_mm, default = first(Zd_mm))) %>% 
    mutate(Distance_mm = sqrt(Delta_X_mm^2 + Delta_Y_mm^2 + Delta_Z_mm^2)) %>% 
    mutate(Speed_cm_s = 0.1*Distance_mm/Delta_time_s) %>% return()
}

wrangled_list <- lapply(zooid_list,wrangle)
speed_data <- do.call(rbind, wrangled_list)

ggplot(speed_data, aes(Time..mins., Speed_cm_s))+geom_line(aes(color=NumberZooid %>% as.factor()))


