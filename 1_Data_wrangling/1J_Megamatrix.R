#1J_Megamatrix
#Author: Mary Lofton
#Date: 30SEP20

#joining all drivers and response variables into a single matrix for analysis

#load packages
#install.packages('pacman')
pacman::p_load(tidyverse, lubridate)
rm(list=ls())

#Hox/EM operation
hoxem <- read_csv("./0_Data_files/HOX_EM_operation.csv") %>%
  mutate(Year = year(Date))

#chemistry

chem <- read_csv("./0_Data_files/chem_vars.csv") %>%
  mutate(Year = year(Date))

#wtr temp and stability and DO
wts <- read_csv("./0_Data_files/WtrTemp_Stability_DO_pH.csv") %>%
  mutate(Year = year(Date)) %>%
  mutate(Temp_DO_Depth_m = ifelse((is.na(CTD_Depth_m) & is.na(YSI_Depth_m)),SCC_Depth_m,
                               ifelse(is.na(CTD_Depth_m),YSI_Depth_m,CTD_Depth_m)),
         Grab_Temp_Depth_m = ifelse((is.na(CTD_Grab_T_Depth_m) & is.na(YSI_Grab_T_Depth_m)),SCC_Grab_T_Depth_m,
                                  ifelse(is.na(CTD_Grab_T_Depth_m),YSI_Grab_T_Depth_m,CTD_Grab_T_Depth_m)),
         Temp_C_grab = ifelse((is.na(CTD_Grab_Depth_Temp_C) & is.na(YSI_Grab_Depth_Temp_C)),SCC_Grab_Depth_Temp_C,
                         ifelse(is.na(CTD_Grab_Depth_Temp_C),YSI_Grab_Depth_Temp_C,CTD_Grab_Depth_Temp_C)),
         Temp_C_Cmax = ifelse((is.na(CTD_Temp_C) & is.na(YSI_Temp_C)),SCC_Temp_C,
                         ifelse(is.na(CTD_Temp_C),YSI_Temp_C,CTD_Temp_C)),
         schmidt.stability = ifelse((is.na(CTD_schmidt.stability) & is.na(YSI_schmidt.stability)),SCC_schmidt.stability,
                                    ifelse(is.na(CTD_schmidt.stability),YSI_schmidt.stability,CTD_schmidt.stability)),
         n2 = ifelse((is.na(CTD_n2) & is.na(YSI_n2)),SCC_n2,
                                    ifelse(is.na(CTD_n2),YSI_n2,CTD_n2)),
         lake.number = ifelse((is.na(CTD_lake.number) & is.na(YSI_lake.number)),SCC_lake.number,
                     ifelse(is.na(CTD_lake.number),YSI_lake.number,CTD_lake.number)),
         wedderburn.number = ifelse((is.na(CTD_wedderburn.number) & is.na(YSI_wedderburn.number)),SCC_wedderburn.number,
                     ifelse(is.na(CTD_wedderburn.number),YSI_wedderburn.number,CTD_wedderburn.number)),
         thermo.depth = ifelse((is.na(CTD_thermo.depth) & is.na(YSI_thermo.depth)),SCC_thermo.depth,
                               ifelse(is.na(CTD_thermo.depth),YSI_thermo.depth,CTD_thermo.depth)),
         DO_mgL = ifelse((is.na(CTD_DO_mgL) & is.na(YSI_DO_mgL)),NA,
                         ifelse(is.na(CTD_DO_mgL),YSI_DO_mgL,CTD_DO_mgL))) %>%
  select(Year,Date, Temp_DO_Depth_m,Grab_Temp_Depth_m,Temp_C_grab,Temp_C_Cmax,schmidt.stability,n2, lake.number, wedderburn.number, thermo.depth, DO_mgL)



#PAR
par <- read_csv("./0_Data_files/Kd.csv") %>%
  mutate(Year = year(Date)) %>%
  mutate(Kd = ifelse(is.na(CTD_Kd) & is.na(YSI_Kd) & is.na(Secchi_Kd),NA,
                             ifelse(is.na(CTD_Kd) & is.na(YSI_Kd),Secchi_Kd,
                                    ifelse(is.na(CTD_Kd),YSI_Kd,CTD_Kd))),
         pz_depth_m = ifelse(is.na(CTD_pz_depth_m) & is.na(YSI_pz_depth_m) & is.na(Secchi_pz_depth_m),NA,
                             ifelse(is.na(CTD_pz_depth_m) & is.na(YSI_pz_depth_m),Secchi_pz_depth_m,
                                    ifelse(is.na(CTD_pz_depth_m),YSI_pz_depth_m,CTD_pz_depth_m))),
         perc_light_thermocline = ifelse(is.na(CTD_perc_light_thermocline),YSI_perc_light_thermocline,CTD_perc_light_thermocline),
         perc_light_Cmax = ifelse(is.na(CTD_perc_light_Cmax),YSI_perc_light_Cmax,CTD_perc_light_Cmax),
         perc_light_grab = ifelse(is.na(CTD_perc_light_grab),YSI_perc_light_grab,CTD_perc_light_grab)) %>%
  select(Year,Date,Kd,perc_light_thermocline, perc_light_Cmax,pz_depth_m,perc_light_grab)

#photic zone temp, DO, pH
pz.vars <- read_csv("./0_Data_files/pz_WtrTemp_DO_pH.csv") %>%
  mutate(Year = year(Date),
         pz_Temp_C = ifelse(is.na(CTD_pz_Temp_C)&is.na(YSI_pz_Temp_C),SCC_pz_Temp_C,
                            ifelse(is.na(CTD_pz_Temp_C),YSI_pz_Temp_C,CTD_pz_Temp_C)),
         pz_DO_mgL = ifelse(is.na(CTD_pz_DO_mgL),YSI_pz_DO_mgL,CTD_pz_DO_mgL),
         interp_pz_depth_m = CTD_Interp_pz_depth_m) %>%
  select(Year, Date, pz_Temp_C, pz_DO_mgL, interp_pz_depth_m)

#WRT
wrt <- read_csv("./0_Data_files/WRT.csv")%>%
  mutate(Year = year(Date)) %>%
  select(Year, Date, WRT_day)

#FP distribution metrics
dist <- read_csv("./0_Data_files/FP_DistributionMetrics.csv") %>%
  mutate(Year = year(Date))

#biodiversity metrics
bd <- read_csv("./0_Data_files/Biodiversity.csv") %>%
  mutate(Year = year(Date))

#community structure
cs <- read_csv("./0_Data_files/Community_structure.csv") %>%
  mutate(Year = year(Date))

#join all vars
mega <- left_join(hoxem, chem, by = c("Year","Date")) 
mega1 <- left_join(mega, wts, by = c("Year","Date")) 
mega2 <- left_join(mega1, par, by = c("Year","Date")) 
mega3 <- left_join(mega2, dist, by = c("Year","Date"))
mega3.1 <- left_join(mega3, pz.vars, by = c("Year","Date"))
mega3.2 <- left_join(mega3.1, wrt, by = c("Year","Date"))
mega4 <- left_join(mega3.2, bd, by = c("Year","Date")) 
mega5 <- left_join(mega4, cs, by = c("Year","Date")) 
mega6 <- mega5[,c(6,1:5,7:73)]

#write megamatrix to file
write.csv(mega6, "./2_Data_analysis/megamatrix.csv",row.names = FALSE)


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
  filter(abs(Peak_depth_m - Temp_DO_Depth_m)>=1.1) %>%
  select(Date, Peak_depth_m, Temp_DO_Depth_m)
#2017-06-05, 2017-07-24, 2017-08-14
bad_dates <- c("2017-06-05","2017-07-24","2017-08-14")
mega8 <- mega7 %>%
  mutate(Temp_C_Cmax = ifelse(Date %in% bad_dates,NA,Temp_C_Cmax),
         DO_mgL = ifelse(Date %in% bad_dates,NA,DO_mgL))

#eliminate phyto vars that won't be included as drivers and limit to date of FP
#samples
colnames(mega8)
mega9 <- mega8[,c(1:49)] %>%
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

#subset phyto samples that are within 1.1 m of Cmax
#using 1.1 to be consistent w/ chem sample cutoff
mega11 <- mega8 %>%
  filter(abs(Peak_depth_m - Phyto_Depth_m) <= 1.1)
drivers <- mega8 %>%
  filter(abs(Peak_depth_m - Phyto_Depth_m) > 1.1) 
colnames(drivers)
drivers[,c(50:72)] <- NA

mega12 <- bind_rows(mega11, drivers) %>%
  arrange(Date)

#then left join to dates from mega10 in such way that end up with NAs 
#in all appropriate places to run an AR1
final_dates <- tibble(mega10$Date)
colnames(final_dates) <- "Date"
mega13 <- left_join(final_dates, mega12, by = "Date") %>%
  mutate(Temp_C_grab = ifelse(Date == "2017-08-14",NA,Temp_C_grab))

check <- mega13 %>%
  filter(abs(Grab_Temp_Depth_m - Phyto_Depth_m) >= 1.1 | abs(Grab_Chem_Depth_m - Phyto_Depth_m) >= 1.1) %>%
  select(Date,Phyto_Depth_m,Grab_Temp_Depth_m,Grab_Chem_Depth_m,Temp_C_grab)
        
#limit to columns needed for phyto comm structure ARIMAs
colnames(mega13)


#write megamatrix for phyto comm structure ARIMAs to file
write.csv(mega13, "./2_Data_analysis/CS_megamatrix.csv",row.names = FALSE)

