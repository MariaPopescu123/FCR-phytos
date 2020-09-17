#1D_PAR
#Author: Mary Lofton
#Date: 16SEP20

#pull in photosynthetically active radiation (PAR) data and calculate:
#1. Kd, or light attenuation coefficient
#2. % surface light at the depth of Cmax - NOT CALCULATING BECAUSE DON'T HAVE INCIDENT LIGHT
#3. % surface light at the depth of nutrients - NOT CALCULATING BECAUSE DON'T HAVE INCIDENT LIGHT

#load packages
#install.packages('pacman')
pacman::p_load(tidyverse, lubridate)
rm(list=ls())

#read in sample dates and depths of phyto samples
sample_info <- read_csv("./00_Data_files/EDI_phytos/phytoplankton.csv") %>%
  select(Date, Depth_m) %>%
  distinct()
sample_info$number <- 1:100

#read in FP data so can match temp profiles to hour FP profiles were taken
replacement_dates <- as.Date(c("2016-07-12","2018-05-24","2019-07-03","2019-07-11","2019-07-18","2019-10-22"))

fp_sample <- read_csv("./00_Data_files/FP.csv")%>%
  mutate(Date = date(DateTime),
         Hour = hour(DateTime)) %>%
  filter(Reservoir == "FCR" & Site == 50) %>%
  filter(Date %in% sample_info$Date | Date %in% replacement_dates) %>%
  arrange(Date)%>%
  select(Date,Hour,Depth_m,TotalConc_ugL) %>%
  filter(complete.cases(.)) %>%
  filter(!(Date == "2016-05-30" & Hour == 10),
         !(Date == "2016-05-30" & Hour == 11),
         !(Date == "2016-05-30" & Hour == 12),
         !(Date == "2016-05-30" & Hour == 13),
         !(Date == "2016-05-30" & Hour == 14),
         !(Date == "2016-05-30" & Hour == 18),
         !(Date == "2016-06-27" & Hour == 14),
         !(Date == "2016-07-25" & Hour == 1),
         !(Date == "2016-07-25" & Hour == 2),
         !(Date == "2016-07-25" & Hour == 12),
         !(Date == "2016-07-25" & Hour == 13),
         !(Date == "2016-07-25" & Hour == 14),
         !(Date == "2016-07-25" & Hour == 15),
         !(Date == "2016-07-25" & Hour == 16),
         !(Date == "2017-07-10" & Hour == 9),
         !(Date == "2017-07-10" & Hour == 13),
         !(Date == "2017-07-24" & Hour == 10),
         !(Date == "2017-08-21" & Hour == 14),
         !(Date == "2017-08-21" & Hour == 15),
         !(Date == "2018-07-16" & Hour == 9),
         !(Date == "2019-05-13" & Hour == 15),
         !(Date == "2019-10-04" & Hour == 17),
         !(Date == "2019-10-11" & Hour == 15)) %>%
  mutate(Year = year(Date)) %>%
  distinct() %>%
  filter(Depth_m <= 9.8) #revisit this later but a lot of 2019 casts have junk at the bottom
#check 2017-07-05
#check 2017-05-29
fp_sample <- fp_sample[-c(2612:2733,3240:3385),]
fp_sample <- fp_sample %>%
  distinct(Date,Hour)
fp_sample <- fp_sample[-81,]

#read in CTD data and limit to PAR
ctd <- read_csv("./00_Data_files/CTD.csv",
                col_types = cols(PAR_umolm2s = col_double())) %>%
  select(Reservoir, Site, Date, Depth_m, PAR_umolm2s) %>%
  filter(Reservoir == "FCR" & Site == 50 & date(Date) %in% sample_info$Date) %>%
  mutate(Hour = hour(Date),
         Date = date(Date))

#select casts that are closest in time to FP casts
#Write a function that returns the closest value
#xv is vector, sv is specific value
closest<-function(xv, sv){
  xv[which.min(abs(xv-sv))]}

final <- ctd[0,]

ctd_dates <- unique(ctd$Date)
for (i in 1:length(ctd_dates)){
  ctd_sample <- subset(ctd, ctd$Date == ctd_dates[i])
  fp_hour <- subset(fp_sample, fp_sample$Date == ctd_dates[i])
  ctd_profile <- ctd_sample[ctd_sample[, "Hour"] == closest(ctd_sample$Hour,fp_hour$Hour[1]),]
  final <- bind_rows(final, ctd_profile)
}

final <- final %>%
  select(Date, Depth_m, PAR_umolm2s)

final <- final[complete.cases(final),]

check <- final %>%
  filter(Depth_m < 0)

ctd_par <- final

#visualize PAR data, retrieve Kd and % light at Cmax and nutrient max
#only calculating Kd for now because don't have incident light for most profiles

##HUGE ISSUE: NO INCIDENT LIGHT DATA FOR CTD CASTS :-( :-(
par_dates <- unique(ctd_par$Date)

final <- matrix(NA, nrow = length(par_dates), ncol = 2)

for (i in 1:length(par_dates)){
  par_profile <- subset(ctd_par, ctd_par$Date == par_dates[i])
  final[i,1] <- as.character(par_dates[i])
  plot(par_profile$Depth_m, log(par_profile$PAR_umolm2s/max(par_profile$PAR_umolm2s, na.rm = TRUE)*100))
  mod <- lm(log(par_profile$PAR_umolm2s/max(par_profile$PAR_umolm2s, na.rm = TRUE)*100)~par_profile$Depth_m, )
  kd <- unlist(mod$coefficients[2])*-1
  final[i,2] <- kd
  
}

final <- data.frame(final) %>%
  mutate(Date = as.Date(X1),CTD_Kd = as.double(X2)) %>%
  select(Date, CTD_Kd)
ctd_kd <- final


#read in YSI data and limit to PAR
ysi <- read_csv("./00_Data_files/YSI.csv") %>%
  select(Reservoir, Site, DateTime, Depth_m, PAR_umolm2s) %>%
  mutate(Date = as.POSIXct(DateTime, format = "%m/%d/%y %H:%M")) %>%
  select(-DateTime)%>%
  filter(Reservoir == "FCR" & Site == 50 & date(Date) %in% sample_info$Date) %>%
  mutate(Hour = hour(Date),
         Date = date(Date))

# #check for sampling day alignment btwn. phytos and YSI profiles
# ysi_dates <- data.frame(unique(ysi$Date))
# colnames(ysi_dates) <- c("DateTime")
# ysi_dates$date <- date(ysi_dates$DateTime)
# ysi_dates1 <- data.frame(unique(ysi_dates$date))
# colnames(ysi_dates1) <- c("Date")
# ysi_dates1$YSI_sample <- TRUE
# check1 <- left_join(check, ysi_dates1, by = "Date")

#select casts that are closest in time to FP casts
#Write a function that returns the closest value
#xv is vector, sv is specific value
closest<-function(xv, sv){
  xv[which.min(abs(xv-sv))]}

final <- ysi[0,]

ysi_dates <- unique(ysi$Date)
for (i in 1:length(ysi_dates)){
  ysi_sample <- subset(ysi, ysi$Date == ysi_dates[i])
  fp_hour <- subset(fp_sample, fp_sample$Date == ysi_dates[i])
  ysi_profile <- ysi_sample[ysi_sample[, "Hour"] == closest(ysi_sample$Hour,fp_hour$Hour[1]),]
  final <- bind_rows(final, ysi_profile)
}

final <- final %>%
  select(Date, Depth_m, PAR_umolm2s)

ysi_par <- data.frame(final)%>%
  mutate(PAR_umolm2s = as.double(PAR_umolm2s))%>%
  filter(!is.na(PAR_umolm2s))

check <- ysi_par %>%
  filter(Depth_m == -0.1)
#oofta - only 32 of 72 casts have incident light and actually it's only
#recorded for about two thirds of those - yikes...

#visualize PAR data, retrieve Kd and % light at Cmax and nutrient max
#only calculating Kd for now because don't have incident light for most profiles

par_dates <- unique(ysi_par$Date)

final <- matrix(NA, nrow = length(par_dates), ncol = 2)

for (i in 1:length(par_dates)){
  par_profile <- subset(ysi_par, ysi_par$Date == par_dates[i])
  final[i,1] <- as.character(par_dates[i])
  plot(par_profile$Depth_m, log((par_profile$PAR_umolm2s+0.001)/max(par_profile$PAR_umolm2s, na.rm = TRUE)*100))
  mod <- lm(log((par_profile$PAR_umolm2s+0.001)/max(par_profile$PAR_umolm2s, na.rm = TRUE)*100)~par_profile$Depth_m, )
  kd <- unlist(mod$coefficients[2])*-1
  final[i,2] <- kd
  
}

final <- data.frame(final) %>%
  mutate(Date = as.Date(X1),YSI_Kd = as.double(X2)) %>%
  select(Date, YSI_Kd)
ysi_kd <- final

Kd <- data.frame(sample_info$Date)
colnames(Kd) <- "Date"
Kd1 <- left_join(Kd, ctd_kd, by = "Date")
Kd2 <- left_join(Kd1, ysi_kd, by = "Date")

check <- Kd2 %>%
  filter(is.na(CTD_Kd) & is.na(YSI_Kd))

check <- Kd2 %>%
  filter(YSI_Kd > 1.25)

#visualization
Kd <- read_csv("./00_Data_files/Kd.csv")

Kd_plot <- Kd %>%
  gather(CTD_Kd:YSI_Kd, key = "sensor",value = "Kd") %>%
  mutate(Year = year(Date))

ggplot(data = Kd_plot, aes(x = Kd, group = sensor, color = sensor, fill = sensor))+
  geom_density(alpha = 0.3)+
  theme_classic()

ggplot(data = Kd_plot, aes(x = Date, y = Kd, color = sensor, group = sensor))+
  facet_wrap(vars(Year), scales = "free")+
  geom_point(size = 2)+
  geom_line(size = 1)+
  theme_classic()

write.csv(Kd2, "./00_Data_files/Kd.csv",row.names = FALSE)
