#1E_PAR
#Author: Mary Lofton
#Date: 16SEP20

#pull in photosynthetically active radiation (PAR) data and calculate:
#1. Kd, or light attenuation coefficient
#2. % surface light at the depth of Cmax - NOT CALCULATED YET 
#3. % surface light at the thermocline - FOR NOW USING LIGHT AT 0.1 M AS INCIDENT LIGHT

#load packages
#install.packages('pacman')
pacman::p_load(tidyverse, lubridate, data.table)
rm(list=ls())


#read in sample dates and depths of phyto samples
sample_info <- read_csv("./0_Data_files/EDI_phytos/phytoplankton.csv") %>%
  select(Date, Depth_m) %>%
  distinct()
sample_info$number <- 1:100

#read in FP data so can match temp profiles to hour FP profiles were taken
replacement_dates <- as.Date(c("2016-07-12","2018-05-24","2019-07-03","2019-07-11","2019-07-18","2019-10-22"))

fp_sample <- read_csv("./0_Data_files/FP.csv")%>%
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
# data  <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.200.11&entityid=d771f5e9956304424c3bc0a39298a5ce"
# 
# destination <- "./0_Data_files"
# 
# download.file(data, destfile = "./0_Data_files/CTD.csv", method='libcurl')
ctd <- fread("./0_Data_files/CTD.csv")
ctd <- data.frame(ctd) %>%
  select(Reservoir, Site, Date, Depth_m, PAR_umolm2s) %>%
  filter(Reservoir == "FCR" & Site == 50 & date(Date) %in% sample_info$Date) %>%
  mutate(Hour = hour(Date),
         Date = date(Date),
         PAR_umolm2s = as.double(PAR_umolm2s))

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

ctd_par <- final %>%
  arrange(Date)

#visualize PAR data, retrieve Kd and % light at Cmax and nutrient max
#only calculating Kd for now because don't have incident light for most profiles

#read in data to pull thermocline depth
wts <- read_csv("./0_Data_files/WtrTemp_Stability.csv") %>%
  mutate(Year = year(Date)) %>%
  mutate(Temp_Depth_m = ifelse((is.na(CTD_Depth_m) & is.na(YSI_Depth_m)),SCC_Depth_m,
                               ifelse((is.na(CTD_Depth_m)),YSI_Depth_m,CTD_Depth_m)),
         Temp_C = ifelse((is.na(CTD_Temp_C) & is.na(YSI_Temp_C)),SCC_Temp_C,
                         ifelse((is.na(CTD_Temp_C)),YSI_Temp_C,CTD_Temp_C)),
         schmidt.stability = ifelse((is.na(CTD_schmidt.stability) & is.na(YSI_schmidt.stability)),SCC_schmidt.stability,
                                    ifelse((is.na(CTD_schmidt.stability)),YSI_schmidt.stability,CTD_schmidt.stability)),
         thermo.depth = ifelse((is.na(CTD_thermo.depth) & is.na(YSI_thermo.depth)),SCC_thermo.depth,
                               ifelse((is.na(CTD_thermo.depth)),YSI_thermo.depth,CTD_thermo.depth))) %>%
  select(Year,Date, Temp_Depth_m,Temp_C,schmidt.stability,thermo.depth)


par_dates <- unique(ctd_par$Date)
thermo_CTD_dates <- wts %>%
  filter(Date %in% par_dates)

#read in cmax depths
cmax <- read_csv("./0_Data_files/FP_DistributionMetrics.csv") %>%
  select(Date, Peak_depth_m) %>%
  filter(Date %in% par_dates)

par_grabs <- sample_info %>%
  filter(Date %in% par_dates)

final <- matrix(NA, nrow = length(par_dates), ncol = 6)

for (i in 1:length(par_dates)){
  par_profile <- subset(ctd_par, ctd_par$Date == par_dates[i])
  final[i,1] <- as.character(par_dates[i])
  plot(par_profile$Depth_m, log(par_profile$PAR_umolm2s/max(par_profile$PAR_umolm2s, na.rm = TRUE)*100))
  mod <- lm(log(par_profile$PAR_umolm2s/max(par_profile$PAR_umolm2s, na.rm = TRUE)*100)~par_profile$Depth_m, )
  kd <- unlist(mod$coefficients[2])*-1
  final[i,2] <- kd
  
  if(!is.na(thermo_CTD_dates$thermo.depth[i])){
    thermo.light.df <- par_profile[par_profile[, "Depth_m"] == closest(par_profile$Depth_m,thermo_CTD_dates$thermo.depth[i]),]
    thermo.light <- thermo.light.df$PAR_umolm2s/max(par_profile$PAR_umolm2s, na.rm = TRUE)*100
    final[i,3] <- thermo.light
  } else {
    final[i,3] <- NA
  }
  
  if(!is.na(cmax$Peak_depth_m[i])){
    cmax.light.df <- par_profile[par_profile[, "Depth_m"] == closest(par_profile$Depth_m,cmax$Peak_depth_m[i]),]
    cmax.light <- cmax.light.df$PAR_umolm2s/max(par_profile$PAR_umolm2s, na.rm = TRUE)*100
    final[i,4] <- cmax.light
  } else {
    final[i,4] <- NA
  }
  
  if(!is.na(par_grabs$Depth_m[i])){ 
    cmax.light.df <- par_profile[par_profile[, "Depth_m"] == closest(par_profile$Depth_m,par_grabs$Depth_m[i]),]
    #final[i,5] <- cmax.light.df$Depth_m
    cmax.light <- cmax.light.df$PAR_umolm2s/max(par_profile$PAR_umolm2s, na.rm = TRUE)*100
    final[i,5] <- cmax.light[1]
  } else {
    final[i,5] <- NA
  }
  
  par_profile$perc_par <- par_profile$PAR_umolm2s/max(par_profile$PAR_umolm2s, na.rm = TRUE)*100
  pz_depth.df <- par_profile[par_profile[, "perc_par"] == closest(par_profile$perc_par,1),]
  pz_depth <- pz_depth.df$Depth_m
  final[i,6] <- pz_depth
  
  
}

final <- data.frame(final) %>%
  mutate(Date = as.Date(X1),CTD_Kd = as.double(X2),
         CTD_perc_light_thermocline = as.double(X3),
         CTD_perc_light_Cmax = as.double(X4),
         CTD_pz_depth_m = as.double(X6),
         CTD_perc_light_grab = as.double(X5)) %>%
  select(Date, CTD_Kd, CTD_perc_light_thermocline, CTD_perc_light_Cmax,
         CTD_pz_depth_m, CTD_perc_light_grab)
ctd_kd <- final

#check <- left_join(ctd_kd,sample_info,by = "Date")
#read in YSI data and limit to PAR
ysi <- read_csv("./0_Data_files/YSI.csv") %>%
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
thermo_YSI_dates <- wts %>%
  filter(Date %in% par_dates)

#read in cmax depths
cmax <- read_csv("./0_Data_files/FP_DistributionMetrics.csv") %>%
  select(Date, Peak_depth_m) %>%
  filter(Date %in% par_dates)

par_grabs <- sample_info %>%
  filter(Date %in% par_dates)

final <- matrix(NA, nrow = length(par_dates), ncol = 6)

for (i in 1:length(par_dates)){
  par_profile <- subset(ysi_par, ysi_par$Date == par_dates[i])
  final[i,1] <- as.character(par_dates[i])
  plot(par_profile$Depth_m, log((par_profile$PAR_umolm2s+0.001)/max(par_profile$PAR_umolm2s, na.rm = TRUE)*100))
  mod <- lm(log((par_profile$PAR_umolm2s+0.001)/max(par_profile$PAR_umolm2s, na.rm = TRUE)*100)~par_profile$Depth_m, )
  kd <- unlist(mod$coefficients[2])*-1
  final[i,2] <- kd
  
  if(!is.na(thermo_YSI_dates$thermo.depth[i])){
    thermo.light.df <- par_profile[par_profile[, "Depth_m"] == closest(par_profile$Depth_m,thermo_YSI_dates$thermo.depth[i]),]
    thermo.light <- thermo.light.df$PAR_umolm2s/max(par_profile$PAR_umolm2s, na.rm = TRUE)*100
    final[i,3] <- thermo.light
  } else {
    final[i,3] <- NA
  }
  
  if(!is.na(cmax$Peak_depth_m[i])){
    cmax.light.df <- par_profile[par_profile[, "Depth_m"] == closest(par_profile$Depth_m,cmax$Peak_depth_m[i]),]
    cmax.light <- cmax.light.df$PAR_umolm2s/max(par_profile$PAR_umolm2s, na.rm = TRUE)*100
    final[i,4] <- cmax.light
  } else {
    final[i,4] <- NA
  }
  
  if(!is.na(par_grabs$Depth_m[i])){
    cmax.light.df <- par_profile[par_profile[, "Depth_m"] == closest(par_profile$Depth_m,par_grabs$Depth_m[i]),]
    cmax.light <- cmax.light.df$PAR_umolm2s/max(par_profile$PAR_umolm2s, na.rm = TRUE)*100
    final[i,5] <- cmax.light
  } else {
    final[i,5] <- NA
  }
  
  par_profile$perc_par <- par_profile$PAR_umolm2s/max(par_profile$PAR_umolm2s, na.rm = TRUE)*100
  pz_depth.df <- par_profile[par_profile[, "perc_par"] == closest(par_profile$perc_par,1),]
  pz_depth <- pz_depth.df$Depth_m[1]
  final[i,6] <- pz_depth
  
}

final <- data.frame(final) %>%
  mutate(Date = as.Date(X1),YSI_Kd = as.double(X2),
         YSI_perc_light_thermocline = as.double(X3),
         YSI_perc_light_Cmax = as.double(X4),
         YSI_pz_depth_m = as.double(X6),
         YSI_perc_light_grab = as.double(X5)) %>%
  select(Date, YSI_Kd, YSI_perc_light_thermocline,
         YSI_perc_light_Cmax,YSI_pz_depth_m,YSI_perc_light_grab)
ysi_kd <- final[-c(27,34),]

#add Secchi as last resort
secchi <- read_csv("./0_Data_files/Secchi.csv") %>%
  filter(Reservoir == "FCR" & Site == 50 & date(DateTime) %in% sample_info$Date)%>%
  mutate(Date = date(DateTime))

secchi$Secchi_Kd = 1.7/secchi$Secchi_m
secchi$Secchi_pz_depth_m = 2.8*secchi$Secchi_m

secchi <- secchi %>%
  select(Date, Secchi_Kd, Secchi_pz_depth_m)

hist(secchi$Secchi_Kd)
hist(secchi$Secchi_pz_depth_m)
hist(ctd_kd$CTD_Kd)
hist(ctd_kd$CTD_pz_depth_m)
hist(ysi_kd$YSI_Kd)
hist(ysi_kd$YSI_pz_depth_m)

#combine Kd from all sources into a single data frame

Kd <- data.frame(sample_info$Date)
colnames(Kd) <- "Date"
Kd1 <- left_join(Kd, ctd_kd, by = "Date")
Kd2 <- left_join(Kd1, ysi_kd, by = "Date")
Kd3 <- left_join(Kd2, secchi, by = "Date") 
Kd3 <- Kd3[-c(30,31),]

check <- Kd3 %>%
  filter(is.na(CTD_Kd) & is.na(YSI_Kd) & is.na(Secchi_Kd))

write.csv(Kd3, "./0_Data_files/Kd.csv",row.names = FALSE)

#visualization
Kd <- read_csv("./0_Data_files/Kd.csv")

Kd_plot <- Kd %>%
  select(Date, CTD_Kd, YSI_Kd)%>%
  gather(CTD_Kd,YSI_Kd, key = "sensor",value = "Kd") %>%
  mutate(Year = year(Date))

ggplot(data = Kd_plot, aes(x = Kd, group = sensor, color = sensor, fill = sensor))+
  geom_density(alpha = 0.3)+
  theme_classic()

ggplot(data = Kd_plot, aes(x = Date, y = Kd, color = sensor, group = sensor))+
  facet_wrap(vars(Year), scales = "free")+
  geom_point(size = 2)+
  geom_line(size = 1)+
  theme_classic()

perc_light_plot <- Kd %>%
  select(Date, CTD_perc_light_thermocline, YSI_perc_light_thermocline)%>%
  gather(CTD_perc_light_thermocline,YSI_perc_light_thermocline, key = "sensor",value = "perc_light_thermocline") %>%
  mutate(Year = year(Date))

ggplot(data = perc_light_plot, aes(x = perc_light_thermocline, group = sensor, color = sensor, fill = sensor))+
  geom_density(alpha = 0.3)+
  theme_classic()

ggplot(data = perc_light_plot, aes(x = Date, y = perc_light_thermocline, color = sensor, group = sensor))+
  facet_wrap(vars(Year), scales = "free")+
  geom_point(size = 2)+
  geom_line(size = 1)+
  theme_classic()

