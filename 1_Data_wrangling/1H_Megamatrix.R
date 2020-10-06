#1H_Megamatrix
#Author: Mary Lofton
#Date: 30SEP20

#joining all drivers and response variables into a single matrix for analysis

#load packages
#install.packages('pacman')
pacman::p_load(tidyverse, lubridate)
rm(list=ls())

cs <- read_csv("./00_Data_files/Community_structure.csv") %>%
  mutate(Year = year(Date))

dist <- read_csv("./00_Data_files/FP_DistributionMetrics.csv") %>%
  mutate(Year = year(Date))

bd <- read_csv("./00_Data_files/Biodiversity.csv") %>%
  mutate(Year = year(Date))

mega1 <- left_join(cs, dist, by = c("Year","Date")) 
mega <- left_join(mega1, bd, by = c("Year","Date")) 

write.csv(mega, "./2_Data_analysis/megamatrix.csv",row.names = FALSE)
