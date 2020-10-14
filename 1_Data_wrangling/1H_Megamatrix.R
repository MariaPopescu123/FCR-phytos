#1H_Megamatrix
#Author: Mary Lofton
#Date: 30SEP20

#joining all drivers and response variables into a single matrix for analysis

#load packages
#install.packages('pacman')
pacman::p_load(tidyverse, lubridate)
rm(list=ls())

#Hox/EM operation
hoxem <- read_csv("./00_Data_files/HOX_EM_operation.csv") %>%
  mutate(Year = year(Date))

#chemistry
chem <- read_csv("./00_Data_files/chem_vars.csv") %>%
  mutate(Year = year(Date))

#wtr temp and stability
wts <- read_csv("./00_Data_files/WtrTemp_Stability.csv") %>%
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



#PAR
par <- read_csv("./00_Data_files/Kd.csv") %>%
  mutate(Year = year(Date)) %>%
  mutate(Kd = ifelse(is.na(CTD_Kd),YSI_Kd,CTD_Kd)) %>%
  select(Year,Date,Kd)

#FP distribution metrics
dist <- read_csv("./00_Data_files/FP_DistributionMetrics.csv") %>%
  mutate(Year = year(Date))

#biodiversity metrics
bd <- read_csv("./00_Data_files/Biodiversity.csv") %>%
  mutate(Year = year(Date))

#community structure
cs <- read_csv("./00_Data_files/Community_structure.csv") %>%
  mutate(Year = year(Date))

#join all vars
mega <- left_join(hoxem, chem, by = c("Year","Date")) 
mega1 <- left_join(mega, wts, by = c("Year","Date")) 
mega2 <- left_join(mega1, par, by = c("Year","Date")) 
mega3 <- left_join(mega2, dist, by = c("Year","Date"))
mega4 <- left_join(mega3, bd, by = c("Year","Date")) 
mega5 <- left_join(mega4, cs, by = c("Year","Date")) 
mega6 <- mega5[,c(4,1:3,5:47)]

#write megamatrix to file
#write.csv(mega6, "./2_Data_analysis/megamatrix.csv",row.names = FALSE)


#DATA WRANGLING FOR FP DISTRIBUTION METRIC ARIMAS

#check depths: samples that are > 1.1 m from Peak_depth_m are replaced w/ NAs for
#that variable

#1.1 m is selected based on the resolution of the chem dataset, where the two
#depth samples that are farthest apart are 2.2 m apart

check <- mega6 %>%
  filter(abs(Peak_depth_m - Chem_Depth_m)>=1.1) %>%
  select(Date, Peak_depth_m, Chem_Depth_m)
#2017-05-01
mega7 <- mega6 %>%
  mutate(Cmax_SRP_ugL = ifelse(Date == "2017-05-01",NA,Cmax_SRP_ugL),
         Cmax_DIN_ugL = ifelse(Date == "2017-05-01",NA,Cmax_DIN_ugL),
         Cmax_DOC_mgL = ifelse(Date == "2017-05-01",NA,Cmax_DOC_mgL))

check1 <- mega6 %>%
  filter(abs(Peak_depth_m - Temp_Depth_m)>=1.1) %>%
  select(Date, Peak_depth_m, Temp_Depth_m)
#2017-06-05, 2017-07-24, 2017-08-14
bad_dates <- c("2017-06-05","2017-07-24","2017-08-14")
mega8 <- mega7 %>%
  mutate(Temp_C = ifelse(Date %in% bad_dates,NA,Temp_C))

#eliminate phyto vars that won't be included as drivers and limit to date of FP
#samples
colnames(mega8)
mega9 <- mega8[,1:23] %>%
  filter(!is.na(Peak_depth_m))

#add rows of NA in between observations that are 2 wks or more apart
#after row 24 (2016-10-07)
#after row 45 (2017-09-25)
#after row 65 (2018-09-17)
#after row 66 (2019-02-27)
#after row 67 (2019-03-18)
#after row 68 (2019-04-01)

myNAs <- mega9 %>%
  filter(Peak_depth_m == 12)
myNAs[1:6,] <- NA
myNAs$Year <- c(2017,2018,2019,2019,2019,2019)
myNAs$Date <- as.Date(c("2017-01-01","2018-01-01","2019-01-01","2019-03-07","2019-03-25","2019-04-08"))

mega10 <- rbind(mega9,myNAs) %>%
  arrange(Year,Date)

#write megamatrix for FP ARIMA to file
write.csv(mega10, "./2_Data_analysis/FP_megamatrix.csv",row.names = FALSE)


##NEED TO GO BACK AND FIND DEPTH-SPECIFIC THINGS BASED ON PHYTO SAMPLE
#SEPARATELY AND THEN CHECK DEPTH CORRELATION OF DRIVERS WITH PHYTO SAMPLES
#AND MAKE SEPARATE MATRIX FOR PHYTO ARIMA