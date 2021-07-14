##1G_FP_Distribution_metrics
##Author: Mary Lofton
##Date: 27SEP18

#load packages
pacman::p_load(tidyverse, lubridate)

#clear environment
rm(list=ls())

#read in sample dates and depths of phyto samples
sample_info <- read_csv("./0_Data_files/EDI_phytos/phytoplankton.csv") %>%
  select(Date, Depth_m) %>%
  distinct()
sample_info$number <- 1:100

#load fp data to get biomass value
fp <- read_csv("./0_Data_files/FP.csv")%>%
  mutate(Date = date(DateTime),
         Hour = hour(DateTime)) %>%
  filter(Reservoir == "FCR" & Site == 50 & date(DateTime) %in% sample_info$Date) 

#check for sampling day alignment btwn. phytos and FP profiles
fp_dates <- data.frame(unique(fp$Date))
colnames(fp_dates) <- c("Date")
fp_dates$FP_sample <- TRUE
check <- left_join(sample_info, fp_dates, by = "Date")

#missing days: 2016-07-11, 2017-01-19, 2018-05-28, 2019-01-21, 2019-07-01,
#2019-07-08, 2019-07-15, 2019-10-23

#possible replacement days: 2016-07-12, 2018-05-24, 2019-07-03, 2019-07-11, 2019-07-18,
#2019-10-22
replacement_dates <- as.Date(c("2016-07-12","2018-05-24","2019-07-03","2019-07-11","2019-07-18","2019-10-22"))


#pull individual depth sample information and clean profiles
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


#check to make sure don't have duplicate casts on same day
check <- fp_sample %>%
  distinct(Date,Hour)
#the "duplicate" on 07-03 is ok because the cast was taken over 2 hours

#calculate depth distribution characteristics for each profile
dates <- unique(fp_sample$Date)

final <- matrix(NA,nrow = 98, ncol = 5)

#use i = 63 for supplement figure example of peak width

for (i in 1:length(dates)){
  
  profile <- fp_sample %>%
    filter(Date == dates[i])
  
  final[i,1]<- as.character(dates[i])
  final[i,2]<- unlist(profile[which.max(profile$TotalConc_ugL),"TotalConc_ugL"])
  final[i,3]<- unlist(profile[which.max(profile$TotalConc_ugL),"Depth_m"])
  final[i,4]<- max(profile$TotalConc_ugL, na.rm = TRUE) - mean(profile$TotalConc_ugL, na.rm = TRUE)
  
  #peak width calculation
  max_depth <- unlist(profile[which.max(profile$TotalConc_ugL),"Depth_m"])
  conc_med <- mean(profile$TotalConc_ugL, na.rm = TRUE)
  
  peak.top.temp <- subset(profile, Depth_m <= max_depth & TotalConc_ugL <= conc_med)
  if(nrow(peak.top.temp) == 0){
    peak.top = min(profile$Depth_m, na.rm = TRUE)
  } else {
    peak.top <- unlist(peak.top.temp[which.min(abs(peak.top.temp$Depth_m - max_depth)),"Depth_m"])
  }

  peak.bottom.temp <- subset(profile, Depth_m >= max_depth & TotalConc_ugL <= conc_med)
  if(nrow(peak.bottom.temp) == 0){
    peak.bottom = max(profile$Depth_m, na.rm = TRUE)
  } else {
    peak.bottom <- unlist(peak.bottom.temp[which.min(abs(peak.bottom.temp$Depth_m - max_depth)),"Depth_m"])
    }
  

  peak.width = abs(peak.top - peak.bottom)
  final[i,5] <- peak.width
  
  max_conc <- max(profile$TotalConc_ugL, na.rm = TRUE)
  max_conc_depth <- as.numeric(final[i,3])
  
}

final <- data.frame(final)

colnames(final) <- c("Date","Max_biomass_ugL","Peak_depth_m","Peak_magnitude_ugL","Peak_width_m")

#tweaks to match up dates with phytoplankton dates
#in cases where profiles were not available from
#the same day we took profiles taken +/- 3
#days from day of phytoplankton sample (this occurred for 6 samples)
final <- final %>%
  mutate(Date = ifelse(Date == "2016-07-12","2016-07-11",
                       ifelse(Date == "2018-05-24","2018-05-28",
                              ifelse(Date == "2019-07-03","2019-07-01",
                                     ifelse(Date == "2019-07-11","2019-07-08",
                                            ifelse(Date == "2019-07-18","2019-07-15",
                                                   ifelse(Date == "2019-10-22","2019-10-23",Date)))))))

final <- final %>%
  mutate(Date = as.Date(Date),
         Max_biomass_ugL = as.numeric(Max_biomass_ugL),
         Peak_depth_m = as.numeric(Peak_depth_m),
         Peak_magnitude_ugL = as.numeric(Peak_magnitude_ugL),
         Peak_width_m = as.numeric(Peak_width_m))
final$Peak_width_m[97:98] <- NA
#note there are only 98 rows here but should be ok when left-join to phyto samples
write.csv(final, "./0_Data_files/FP_DistributionMetrics.csv", row.names = FALSE)

