#1F_Photic_zone_WtrTemp_DO_pH
#Author: Mary Lofton
#Date: 22FEB21

##load packages
#install.packages('pacman')
pacman::p_load(tidyverse, lubridate, data.table,zoo)
rm(list=ls())

#read in sample dates and depths of phyto samples
sample_info <- read_csv("./0_Data_files/phytoplankton.csv") %>%
  select(Date, Depth_m) %>%
  distinct()
sample_info$number <- 1:100

#read in FP data so can match temp profiles to hour FP profiles were taken
replacement_dates <- as.Date(c("2016-07-12","2018-05-24","2019-07-03","2019-07-11","2019-07-18","2019-10-22"))

fp_sample <- read.csv("./0_Data_files/FP.csv")%>%
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

#read in CTD data and limit to temperature, DO, pH
ctd <- fread("./0_Data_files/CTD.csv")
ctd <- tibble(ctd) %>%
  select(Reservoir, Site, Date, Depth_m, Temp_C, DO_mgL, pH) %>%
  filter(Reservoir == "FCR" & Site == 50 & date(Date) %in% sample_info$Date) %>%
  mutate(Hour = hour(Date),
         Date = date(Date))

# #check for sampling day alignment btwn. phytos and CTD profiles
# ctd_dates <- data.frame(unique(ctd$Date))
# colnames(ctd_dates) <- c("DateTime")
# ctd_dates$date <- date(ctd_dates$DateTime)
# ctd_dates1 <- data.frame(unique(ctd_dates$date))
# colnames(ctd_dates1) <- c("Date")
# ctd_dates1$CTD_sample <- TRUE
# check <- left_join(sample_info, ctd_dates1, by = "Date")

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
  select(Date, Depth_m, Temp_C, DO_mgL, pH)

#trim depths to depth of photic zone
pz <- read_csv("./0_Data_files/Kd.csv") %>%
  mutate(pz_depth_m = ifelse(is.na(CTD_pz_depth_m) & is.na(YSI_pz_depth_m) & is.na(Secchi_pz_depth_m),NA,
                             ifelse(is.na(CTD_pz_depth_m) & is.na(YSI_pz_depth_m),Secchi_pz_depth_m,
                                    ifelse(is.na(CTD_pz_depth_m),YSI_pz_depth_m,CTD_pz_depth_m)))) %>%
  select(Date, pz_depth_m) %>%
  filter(Date %in% sample_info$Date) %>%
  mutate(pz_depth_m_interp = na.approx(pz_depth_m, method = "linear"))

final <- ctd[0,]

ctd_dates <- unique(ctd$Date)
for (i in 1:length(ctd_dates)){
  ctd_sample <- subset(ctd, ctd$Date == ctd_dates[i])
  pz_sample <- subset(pz, pz$Date == ctd_dates[i])
  ctd_profile <- ctd_sample %>%
    filter(Depth_m <= pz_sample$pz_depth_m_interp)%>%
    mutate(Depth_m = pz_sample$pz_depth_m_interp)%>%
    group_by(Reservoir, Site, Date) %>%
    summarize(Temp_C = mean(Temp_C, na.rm = TRUE),
              DO_mgL = mean(DO_mgL, na.rm = TRUE),
              pH = mean(pH, na.rm = TRUE))
  final <- bind_rows(final, ctd_profile)
}
#final[NaN]<-NA

final <- bind_cols(final, pz_sample$pz_depth_m_interp)

CTD_pz_vars <- final %>%
  select(-Depth_m,-Hour,-Reservoir,-Site)
colnames(CTD_pz_vars)<- c("Date","pz_Temp_C","pz_DO_mgL","pz_pH","Interp_pz_depth_m")

#read in YSI data and limit to temperature, DO, pH
ysi <- read.csv("./0_Data_files/YSI.csv") %>%
  select(Reservoir, Site, DateTime, Depth_m, Temp_C, DO_mgL, pH) %>%
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
  select(Date, Depth_m, Temp_C, DO_mgL,pH)

final <- ysi[0,]

ysi_dates <- unique(ysi$Date)
for (i in 1:length(ysi_dates)){
  ysi_sample <- subset(ysi, ysi$Date == ysi_dates[i])
  pz_sample <- subset(pz, pz$Date == ysi_dates[i])
  ysi_profile <- ysi_sample %>%
    filter(Depth_m <= pz_sample$pz_depth_m_interp)%>%
    mutate(Depth_m = pz_sample$pz_depth_m_interp)%>%
    group_by(Reservoir, Site, Date) %>%
    summarize(Temp_C = mean(Temp_C, na.rm = TRUE),
              DO_mgL = mean(DO_mgL, na.rm = TRUE),
              pH = mean(pH, na.rm = TRUE))
  final <- bind_rows(final, ysi_profile)
}

final <- bind_cols(final, pz_sample$pz_depth_m_interp)

ysi_pz_vars <- final %>%
  select(-Depth_m,-Hour,-Reservoir,-Site)
colnames(ysi_pz_vars)<- c("pz_Temp_C","pz_DO_mgL","pz_pH","Date","Interp_pz_depth_m")
ysi_pz_vars <- ysi_pz_vars[,-5]

#read in SCC data
scc <- read_csv("./0_Data_files/SCC.csv") %>%
  select(DateTime:ThermistorTemp_C_9) %>%
  filter(date(DateTime) %in% sample_info$Date)

fun1 <- function(lst, n){
  sapply(lst, `[`, n)
}

colnames(scc)[2:11] <- fun1(strsplit(colnames(scc)[2:11], split = "_"),3)
colnames(scc)[1:2] <- c("datetime",0.1)
colnames(scc)[2:11] <- paste("wtr", colnames(scc)[2:11], sep = "_")

scc1 <- scc %>%
  mutate(Date = date(datetime),Hour = hour(datetime))

#Write a function that returns the closest value
#xv is vector, sv is specific value
closest<-function(xv, sv){
  xv[which.min(abs(xv-sv))]}

final <- scc1[0,]

scc_dates <- unique(scc1$Date)
for (i in 1:length(scc_dates)){
  scc_sample <- subset(scc1, scc1$Date == scc_dates[i])
  fp_hour <- subset(fp_sample, fp_sample$Date == scc_dates[i])
  scc_profile <- scc_sample[scc_sample[, "Hour"] == closest(scc_sample$Hour,fp_hour$Hour[1]),]
  final <- bind_rows(final, scc_profile[1,])
}

final <- final[complete.cases(final), ]

final <- final %>%
  mutate(datetime = Date) %>%
  select(-Date, -Hour)

scc_new <- final

wtr_scc <- scc_new

scc_long <- wtr_scc %>%
  gather(wtr_0.1:wtr_9, key = "Depth",value = "Temp_C")%>%
  mutate(Depth_m = fun1(strsplit(Depth, split = "_"),2)) %>%
  select(-Depth)

scc_long <- arrange(scc_long, datetime)

final <- scc_long[0,]

scc_dates <- unique(scc_long$datetime)
for (i in 1:length(scc_dates)){
  scc_sample <- subset(scc_long, scc_long$datetime == scc_dates[i])
  pz_sample <- subset(pz, pz$Date == scc_dates[i])
  scc_profile <- scc_sample %>%
    filter(Depth_m <= pz_sample$pz_depth_m_interp)%>%
    mutate(Depth_m = pz_sample$pz_depth_m_interp)%>%
    group_by(datetime) %>%
    summarize(Temp_C = mean(Temp_C, na.rm = TRUE))
  final <- bind_rows(final, scc_profile)
}

scc_pz_vars <- final %>%
  select(-Depth_m)
colnames(scc_pz_vars)<- c("Date","pz_Temp_C")

#combine data from all three sources (CTD, YSI, SCC thermistor string) into a single data frame
Tmetrics <- data.frame(sample_info$Date)
colnames(Tmetrics) <- "Date"

Tmetrics$Date <- as.Date(Tmetrics$Date)
CTD_pz_vars$Date <- as.Date(CTD_pz_vars$Date)

Tmetrics1 <- left_join(Tmetrics, CTD_pz_vars, by = "Date")  
colnames(Tmetrics1)[2:5] <- paste("CTD", colnames(Tmetrics1)[2:5], sep = "_")
Tmetrics2 <- left_join(Tmetrics1, ysi_pz_vars, by = "Date")
colnames(Tmetrics2)[6:8] <- paste("YSI", colnames(Tmetrics2)[6:8], sep = "_")
Tmetrics3 <- left_join(Tmetrics2, scc_pz_vars, by = "Date")
colnames(Tmetrics3)[9] <- paste("SCC", colnames(Tmetrics3)[9], sep = "_")



write.csv(Tmetrics3,"./0_Data_files/pz_WtrTemp_DO_pH.csv",row.names = FALSE)

