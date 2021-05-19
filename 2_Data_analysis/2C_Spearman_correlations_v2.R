#2C_Spearman_correlations
#Date: 18APR21
#Author: Mary Lofton

#chemistry
chem <- read_csv("./00_Data_files/chem_vars.csv") %>%
  mutate(Year = year(Date))

#wtr temp and stability and DO
wts <- read_csv("./00_Data_files/WtrTemp_Stability_DO_pH.csv") %>%
  mutate(Year = year(Date)) %>%
  mutate(Temp_DO_Depth_m = ifelse((is.na(CTD_Depth_m) & is.na(YSI_Depth_m)),SCC_Depth_m,
                                  ifelse(is.na(CTD_Depth_m),YSI_Depth_m,CTD_Depth_m)),
         Temp_C = ifelse((is.na(CTD_Temp_C) & is.na(YSI_Temp_C)),SCC_Temp_C,
                         ifelse(is.na(CTD_Temp_C),YSI_Temp_C,CTD_Temp_C)),
         schmidt.stability = ifelse((is.na(CTD_schmidt.stability) & is.na(YSI_schmidt.stability)),SCC_schmidt.stability,
                                    ifelse(is.na(CTD_schmidt.stability),YSI_schmidt.stability,CTD_schmidt.stability)),
         n2 = ifelse((is.na(CTD_n2) & is.na(YSI_n2)),SCC_n2,
                     ifelse(is.na(CTD_n2),YSI_n2,CTD_n2)),
         thermo.depth = ifelse((is.na(CTD_thermo.depth) & is.na(YSI_thermo.depth)),SCC_thermo.depth,
                               ifelse(is.na(CTD_thermo.depth),YSI_thermo.depth,CTD_thermo.depth)),
         DO_mgL = ifelse((is.na(CTD_DO_mgL) & is.na(YSI_DO_mgL)),NA,
                         ifelse(is.na(CTD_DO_mgL),YSI_DO_mgL,CTD_DO_mgL)),
         pH = ifelse((is.na(CTD_pH) & is.na(YSI_pH)),NA,
                     ifelse(is.na(CTD_pH),YSI_pH,CTD_pH))) %>%
  select(Year,Date, Temp_DO_Depth_m,Temp_C,schmidt.stability,n2, thermo.depth, DO_mgL,pH)



#PAR
par <- read_csv("./00_Data_files/Kd.csv") %>%
  mutate(Year = year(Date)) %>%
  mutate(Kd = ifelse(is.na(CTD_Kd) & is.na(YSI_Kd) & is.na(Secchi_Kd),NA,
                     ifelse(is.na(CTD_Kd) & is.na(YSI_Kd),Secchi_Kd,
                            ifelse(is.na(CTD_Kd),YSI_Kd,CTD_Kd))),
         pz_depth_m = ifelse(is.na(CTD_pz_depth_m) & is.na(YSI_pz_depth_m) & is.na(Secchi_pz_depth_m),NA,
                             ifelse(is.na(CTD_pz_depth_m) & is.na(YSI_pz_depth_m),Secchi_pz_depth_m,
                                    ifelse(is.na(CTD_pz_depth_m),YSI_pz_depth_m,CTD_pz_depth_m))),
         perc_light_thermocline = ifelse(is.na(CTD_perc_light_thermocline),YSI_perc_light_thermocline,CTD_perc_light_thermocline),
         perc_light_Cmax = ifelse(is.na(CTD_perc_light_Cmax),YSI_perc_light_Cmax,CTD_perc_light_Cmax)) %>%
  select(Year,Date,Kd,perc_light_thermocline, perc_light_Cmax,pz_depth_m)

#photic zone temp, DO, pH
pz.vars <- read_csv("./00_Data_files/pz_WtrTemp_DO_pH.csv") %>%
  mutate(Year = year(Date),
         pz_Temp_C = ifelse(is.na(CTD_pz_Temp_C)&is.na(YSI_pz_Temp_C),SCC_pz_Temp_C,
                            ifelse(is.na(CTD_pz_Temp_C),YSI_pz_Temp_C,CTD_pz_Temp_C)),
         pz_DO_mgL = ifelse(is.na(CTD_pz_DO_mgL),YSI_pz_DO_mgL,CTD_pz_DO_mgL),
         pz_pH = ifelse(is.na(CTD_pz_pH),YSI_pz_pH,CTD_pz_pH),
         interp_pz_depth_m = CTD_Interp_pz_depth_m) %>%
  select(Year, Date, pz_Temp_C)

#FP distribution metrics
dist <- read_csv("./00_Data_files/FP_DistributionMetrics.csv") %>%
  mutate(Year = year(Date))

#join all vars
mega <- left_join(dist, chem, by = c("Year","Date")) 
mega1 <- left_join(mega, wts, by = c("Year","Date")) 
mega2 <- left_join(mega1, par, by = c("Year","Date")) 
mega3 <- left_join(mega2, pz.vars, by = c("Year","Date")) 

colnames(mega3)
#peak depth cols
#peak depth = 3, drivers = c(11,13,15,17,18,19,27,33,34)

#peak width cols
#peak width = 5, drivers = c(17,18,19,22,23,27)

#max biomass cols
#max biomass = 2, drivers = c(8,9,10,21,29)

#peak depth correlations
driver_cols <- c(11,13,15,17,18,19,27,33,34)
pd_final <- matrix(NA, nrow = length(driver_cols), ncol = 7)

for (i in 1:length(driver_cols)){
  
  pd_final[i,1] <- "Peak_depth_m"
  pd_final[i,2] <- colnames(mega3)[driver_cols[i]]
  pd_final[i,3] <- round(cor(mega3[,3],mega3[,driver_cols[i]],method = "spearman",use = "complete.obs"),2)
  
  d <- subset(mega3, mega3$Year == 2016)
  pd_final[i,4] <- round(cor(d[,3],d[,driver_cols[i]],method = "spearman",use = "complete.obs"),2)
  
  d <- subset(mega3, mega3$Year == 2017)
  pd_final[i,5] <- round(cor(d[,3],d[,driver_cols[i]],method = "spearman",use = "complete.obs"),2)
  
  d <- subset(mega3, mega3$Year == 2018)
  pd_final[i,6] <- round(cor(d[,3],d[,driver_cols[i]],method = "spearman",use = "complete.obs"),2)
  
  d <- subset(mega3, mega3$Year == 2019)
  pd_final[i,7] <- round(cor(d[,3],d[,driver_cols[i]],method = "spearman",use = "complete.obs"),2)
  

}

#peak width correlations
driver_cols <- c(17,18,19,22,23,27)
pw_final <- matrix(NA, nrow = length(driver_cols), ncol = 7)

for (i in 1:length(driver_cols)){
  
  pw_final[i,1] <- "Peak_width_m"
  pw_final[i,2] <- colnames(mega2)[driver_cols[i]]
  pw_final[i,3] <- round(cor(mega2[,5],mega2[,driver_cols[i]],method = "spearman",use = "complete.obs"),2)
  
  d <- subset(mega2, mega2$Year == 2016)
  pw_final[i,4] <- round(cor(d[,5],d[,driver_cols[i]],method = "spearman",use = "complete.obs"),2)
  
  d <- subset(mega2, mega2$Year == 2017)
  pw_final[i,5] <- round(cor(d[,5],d[,driver_cols[i]],method = "spearman",use = "complete.obs"),2)
  
  d <- subset(mega2, mega2$Year == 2018)
  pw_final[i,6] <- round(cor(d[,5],d[,driver_cols[i]],method = "spearman",use = "complete.obs"),2)
  
  d <- subset(mega2, mega2$Year == 2019)
  pw_final[i,7] <- round(cor(d[,5],d[,driver_cols[i]],method = "spearman",use = "complete.obs"),2)
  
}

#max biomass correlations
driver_cols <- c(8,9,10,21,29)
mb_final <- matrix(NA, nrow = length(driver_cols), ncol = 7)

for (i in 1:length(driver_cols)){
  
  mb_final[i,1] <- "Max_biomass_ugL"
  mb_final[i,2] <- colnames(mega2)[driver_cols[i]]
  mb_final[i,3] <- round(cor(mega2[,2],mega2[,driver_cols[i]],method = "spearman",use = "complete.obs"),2)
  
  d <- subset(mega2, mega2$Year == 2016)
  mb_final[i,4] <- round(cor(d[,2],d[,driver_cols[i]],method = "spearman",use = "complete.obs"),2)
  
  d <- subset(mega2, mega2$Year == 2017)
  mb_final[i,5] <- round(cor(d[,2],d[,driver_cols[i]],method = "spearman",use = "complete.obs"),2)
  
  d <- subset(mega2, mega2$Year == 2018)
  mb_final[i,6] <- round(cor(d[,2],d[,driver_cols[i]],method = "spearman",use = "complete.obs"),2)
  
  d <- subset(mega2, mega2$Year == 2019)
  mb_final[i,7] <- round(cor(d[,2],d[,driver_cols[i]],method = "spearman",use = "complete.obs"),2)
  
}
