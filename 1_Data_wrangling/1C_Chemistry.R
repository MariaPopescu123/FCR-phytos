#1C_Chemistry
#Author: Mary Lofton
#Date: 23JUN20

#pull in chemistry data from EDI and subset to relevant dates/drivers

#load packages
#install.packages('pacman')
pacman::p_load(tidyverse, lubridate, zoo)
rm(list=ls())

#read in sample dates and depths of phyto samples
sample_info <- read_csv("./0_Data_files/EDI_phytos/phytoplankton.csv") %>%
  select(Date, Depth_m) %>%
  distinct()
sample_info$number <- 1:100

cmax <- read_csv("./0_Data_files/FP_DistributionMetrics.csv") %>%
  select(Date, Peak_depth_m)

#read in chemistry data and isolate to Falling Creek Reservoir and correct sampling site
chem <- read_csv("./0_Data_files/chemistry.csv",
                 col_types = cols(DIC_mgL = col_double(),
                                  DC_mgL = col_double(),
                                  DN_mgL = col_double(),
                                  Flag_DIC = col_double(),
                                  Flag_DC = col_double(),
                                  Flag_DN = col_double())) %>%
  filter(Reservoir == "FCR" & Site == 50 & date(DateTime) %in% sample_info$Date) %>%
  filter(!Depth_m == 10)

#check for sampling day alignment btwn. phytos and chem
chem_dates <- data.frame(unique(chem$DateTime))
chem_dates1 <- data.frame(date(unique(chem$DateTime)))
colnames(chem_dates1) <- "Date"
chem_dates1$chem_samples <- TRUE
check <- left_join(cmax, chem_dates1, by = "Date")
check <- left_join(sample_info, chem_dates1, by = "Date")

#read in PAR metrics to get photic zone depth
par <- read_csv("./0_Data_files/Kd.csv") %>%
  mutate(Year = year(Date)) %>%
  mutate(Kd = ifelse(is.na(CTD_Kd) & is.na(YSI_Kd) & is.na(Secchi_Kd),NA,
                     ifelse(is.na(CTD_Kd) & is.na(YSI_Kd),Secchi_Kd,
                            ifelse(is.na(CTD_Kd),YSI_Kd,CTD_Kd))),
         pz_depth_m = ifelse(is.na(CTD_pz_depth_m) & is.na(YSI_pz_depth_m) & is.na(Secchi_pz_depth_m),NA,
                             ifelse(is.na(CTD_pz_depth_m) & is.na(YSI_pz_depth_m),Secchi_pz_depth_m,
                                    ifelse(is.na(CTD_pz_depth_m),YSI_pz_depth_m,CTD_pz_depth_m))),
         perc_light_thermocline = ifelse(is.na(CTD_perc_light_thermocline),YSI_perc_light_thermocline,CTD_perc_light_thermocline),
         perc_light_Cmax = ifelse(is.na(CTD_perc_light_Cmax),YSI_perc_light_Cmax,CTD_perc_light_Cmax)) %>%
  select(Year,Date,Kd,perc_light_thermocline, perc_light_Cmax,pz_depth_m) %>%
  mutate(pz_depth_m_interp = na.approx(pz_depth_m, method = "linear")) %>%
  select(Date, pz_depth_m_interp)



#two sampling days missing :-(
#2017-01-19
#2019-02-27

#two sampling days w/ duplicate chem. profiles -- using a.m. one before EM turned on
#eliminating
#2016-05-30 16:00:00
#2016-07-25 00:00:00

#subset to correct days and vars of interest
chem1 <- chem %>%
  filter(!DateTime  == chem_dates$unique.chem.DateTime.[7] & !DateTime == chem_dates$unique.chem.DateTime.[15] & !DateTime == chem_dates$unique.chem.DateTime.[68]) %>%
  mutate(DIN_ugL = NO3NO2_ugL + NH4_ugL) %>%
  select(DateTime, Depth_m, SRP_ugL, DIN_ugL, DOC_mgL)

#pull nuts at depth of Cmax and depth of max. nuts
nut_dates <- unique(chem1$DateTime)
sample_dates <- cmax[-66,] %>% #get rid of dates not in chem data
  select(Date, Peak_depth_m) %>%
  rename(Depth_m = Peak_depth_m)


final <- matrix(NA, nrow = length(nut_dates), ncol = 21)

for (i in 1:length(nut_dates)){
  chem_profile <- subset(chem1, DateTime == nut_dates[i])
  
  phyto_sample <- subset(sample_dates, Date == sample_dates$Date[i])
  
  final[i,1] <- as.character(phyto_sample$Date)
  
  Cmax <- phyto_sample$Depth_m
  
  #pull chem at Cmax
  Cmax_chem <- chem_profile[which.min(abs(chem_profile$Depth_m - Cmax)),]
  
  final[i,2] <- unlist(Cmax_chem[,2])
  final[i,3] <- unlist(Cmax_chem[,3])
  final[i,4] <- unlist(Cmax_chem[,4])
  final[i,5] <- unlist(Cmax_chem[,5])
  
  #pull chem at depth of phyto sample
  grab_sample <- subset(sample_info, Date == sample_dates$Date[i])
  grab_depth <- grab_sample$Depth_m
  
  grab_chem <- chem_profile[which.min(abs(chem_profile$Depth_m - grab_depth)),]
  
  
  final[i,6] <- unlist(grab_chem[,2])
  final[i,7] <- unlist(grab_chem[,3])
  final[i,8] <- unlist(grab_chem[,4])
  final[i,9] <- unlist(grab_chem[,5])
  
  #pull depth of max chem w/in photic zone
  pz <- subset(par, Date == sample_dates$Date[i])
  pz_chem <- subset(chem_profile, chem_profile$Depth_m <= pz$pz_depth_m_interp[1])
  
  #if no clear chem max (i.e. multiple depths with same chem value)
  #then that cell is populated w/ NA
  srp <- subset(pz_chem, SRP_ugL == max(SRP_ugL, na.rm = TRUE))
  if(length(unlist(srp[,2]))>1){
    final[i,10] <- NA
  } else {
    final[i,10] <- unlist(srp[,2])
    final[i,11] <- unlist(srp[,3])
  }
  
  din <- subset(pz_chem, DIN_ugL == max(DIN_ugL, na.rm = TRUE))
  if(length(unlist(din[,2]))>1){
    final[i,12] <- NA
  } else {
    final[i,12] <- unlist(din[,2])
    final[i,13] <- unlist(din[,4])
    
  }

  doc <- subset(pz_chem, DOC_mgL == max(DOC_mgL, na.rm = TRUE))
  if(length(unlist(doc[,2]))>1){
    final[i,14] <- NA
  } else {
    final[i,14] <- unlist(doc[,2])
    final[i,15] <- unlist(doc[,5])
  }
  
  #calculate mean and CV of nutrients across the photic zone
  final[i,16] <- mean(pz_chem$SRP_ugL, na.rm = TRUE)
  final[i,17] <- mean(pz_chem$DIN_ugL, na.rm = TRUE)
  final[i,18] <- mean(pz_chem$DOC_mgL, na.rm = TRUE)
  
  final[i,19] <- round(sd(pz_chem$SRP_ugL, na.rm = TRUE)/mean(pz_chem$SRP_ugL, na.rm = TRUE),2)
  final[i,20] <- round(sd(pz_chem$DIN_ugL, na.rm = TRUE)/mean(pz_chem$DIN_ugL, na.rm = TRUE),2)
  final[i,21] <- round(sd(pz_chem$DOC_mgL, na.rm = TRUE)/mean(pz_chem$DOC_mgL, na.rm = TRUE),2)
  

}

final <- data.frame(final)
colnames(final) <- c("Date","Chem_Depth_m","Cmax_SRP_ugL","Cmax_DIN_ugL","Cmax_DOC_mgL","Grab_Chem_Depth_m","Grab_SRP_ugL","Grab_DIN_ugL","Grab_DOC_mgL","SRPmax_depth_m","SRPmax_ugL","DINmax_depth_m","DINmax_ugL","DOCmax_depth_m","DOCmax_mgL","pz_SRP_mean","pz_DIN_mean","pz_DOC_mean","pz_SRP_CV","pz_DIN_CV","pz_DOC_CV")

#note this only has 97 rows but should still be able to left_join to final dataframe
write.csv(final, file = "./0_Data_files/chem_vars.csv",row.names = FALSE)

#visualization
chem <- read_csv("./0_Data_files/chem_vars.csv") %>%
  mutate(Year = year(Date))

ggplot(data = chem, aes(x = Date, y = Cmax_SRP_ugL))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = chem, aes(x = Date, y = Cmax_DIN_ugL))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = chem, aes(x = Date, y = Cmax_DOC_mgL))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = chem, aes(x = Date, y = SRPmax_depth_m))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()+
  scale_y_reverse()

ggplot(data = chem, aes(x = Date, y = DINmax_depth_m))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()+
  scale_y_reverse()

ggplot(data = chem, aes(x = Date, y = DOCmax_depth_m))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()+
  scale_y_reverse()

ggplot(data = chem, aes(x = Date, y = SRPmax_ugL))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = chem, aes(x = Date, y = DINmax_ugL))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = chem, aes(x = Date, y = DOCmax_mgL))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()

##distribution of epilimnetic nutrients
epi <- chem1 %>%
  filter(Depth_m <= 5)
hist(epi$SRP_ugL)
hist(epi$DIN_ugL)
hist(epi$DOC_mgL)
hist(cmax$Peak_depth_m, main = "",xlab = "Cmax depths (m)")

