#1K_WRT
#Author: Mary Lofton
#Date: 22FEB21

#load packages
#install.packages('pacman')
pacman::p_load(tidyverse, lubridate, data.table)
rm(list=ls())

#read in discharge data
inf <- read_csv("./0_Data_files/discharge.csv") %>%
  select(DateTime, WVWA_Flow_cms) %>%
  mutate(Date = date(DateTime)) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") %>%
  filter(year(Date) %in% c(2016:2019)) %>%
  group_by(Date) %>%
  summarize(mean_inf_cms = mean(WVWA_Flow_cms, na.rm = TRUE))

#calculate WRT
FCR_vol_cm <- 310000

inf$WRT_sec <- FCR_vol_cm/inf$mean_inf_cms
inf$WRT_day <- inf$WRT_sec/(60*60*24)
hist(inf$WRT_day)

ggplot(data = inf, aes(x = Date, y = WRT_day))+
  geom_point()+
  geom_line()+
  theme_classic()

write.csv(inf, "./0_Data_files/WRT.csv",row.names = FALSE)

