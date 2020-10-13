#1A_HOX_EM_operation
#Author: Mary Lofton
#Date: 22JUN20

#creates binary timeseries reflecting operation of hypolimnetic oxygenation system
#and epilimnetic mixing system in FCR from 2016-2019

####NOTE starting dates of HOx are sometimes before phyto sampling starts;
#ending dates of HOx are sometimes after phyto sampling ends
#dates here reflect when first phyto samples collected if HOx was already on
#or last samples collected if HOx stayed on after last phyto sample

#2016
#HOx operation
#2016-04-18 through 2016-10-09

#EM operation
#2016-05-30 (6 hrs)
#2016-06-27 (24 hrs)
#2016-07-25 (intermittant)
#2016-07-26 (intermittant)
#2016-07-27 (intermittant)

#2017
#HOx operation
#2017-04-18 through 2017-10-25

#EM operation
#2017-05-29 (24 hrs)
#2017-07-10 (8-10 hrs)
#2017-07-11 (8-10 hrs)
#2017-07-12 (8-10 hrs)

#2018
#HOx operation
#2018-04-24 through 2018-07-19

#2019 
#HOx operation
#2019-06-03 through 2019-06-17
#2019-07-08 through 2019-07-22
#2019-08-05 through 2019-08-19
#2019-09-02 through 2019-11-20

#load packages
#install.packages('pacman')
pacman::p_load(tidyverse, lubridate)
rm(list=ls())


#vector of EM operation dates
EM_dates <- as.Date(c("2016-05-30", "2016-06-06","2016-06-13","2016-06-27",
              "2016-07-05","2016-07-11","2016-07-25","2016-08-01","2016-08-08",
              "2017-05-29","2017-06-05","2017-06-12", "2017-07-10",
              "2017-07-17","2017-07-24"))
#these are including dates of EM + 2 wks after
              
#get sampling dates and add columns for HOx/EM operation
hoxem <- read_csv("./00_Data_files/EDI_phytos/phytoplankton.csv") %>%
  select(Date) %>%
  distinct() %>%
  add_column(HOx = 0) %>%
  mutate(HOx = ifelse((Date > "2016-04-18" & Date < "2016-10-09") | 
                            (Date > "2017-04-18" & Date < "2017-10-25") | 
                            (Date > "2018-04-24" & Date < "2018-07-19") |
                            (Date > "2019-06-03" & Date <= "2019-06-17") | 
                            (Date > "2019-07-08" & Date <= "2019-07-22") |
                            (Date > "2019-08-05" & Date <= "2019-08-19") |
                            (Date > "2019-09-02" & Date <= "2019-11-20"),1,0),
         EM = ifelse(Date %in% EM_dates,1,0))

#check to make sure 1s populated properly
plot(hoxem$Date, hoxem$HOx)
plot(hoxem$Date, hoxem$EM)

#write to file
write.csv(hoxem, file = "./00_Data_files/HOx_EM_operation.csv", row.names = FALSE)

#preliminary visualization
hoxem <- read_csv("./00_Data_files/HOx_EM_operation.csv") %>%
  mutate(Year = as.factor(year(Date)))

ggplot(data = hoxem, aes(x = Date, y = HOx))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = hoxem, aes(x = Date, y = EM))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()
