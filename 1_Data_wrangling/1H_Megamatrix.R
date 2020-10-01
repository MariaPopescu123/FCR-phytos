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


ggplot(data = mega, aes(x = Peak_depth_m, y = rel_abund_cyano, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_depth_m, y = rel_abund_crypto, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_depth_m, y = rel_abund_green, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_depth_m, y = rel_abund_brown, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_depth_m, y = richness, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_depth_m, y = shannon, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_magnitude_ugL, y = richness, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)+
  theme_classic()
ggplot(data = mega, aes(x = Peak_magnitude_ugL, y = richness))+
  geom_point()+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)+
  theme_classic()
median(mega$Peak_magnitude_ugL, na.rm = TRUE)
mean(mega$Peak_magnitude_ugL, na.rm = TRUE)
sd(mega$Peak_magnitude_ugL, na.rm = TRUE)



ggplot(data = mega, aes(x = Peak_magnitude_ugL, y = shannon, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_magnitude_ugL, y = rel_abund_cyano, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_magnitude_ugL, y = rel_abund_crypto, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_magnitude_ugL, y = rel_abund_green, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_magnitude_ugL, y = rel_abund_brown, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_width_m, y = rel_abund_cyano, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_width_m, y = rel_abund_crypto, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_width_m, y = rel_abund_brown, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_width_m, y = rel_abund_green, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_width_m, y = richness, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_width_m, y = shannon, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()
