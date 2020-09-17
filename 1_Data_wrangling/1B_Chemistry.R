#1B_Chemistry
#Author: Mary Lofton
#Date: 23JUN20

#pull in chemistry data from EDI and subset to relevant dates/drivers

#load packages
#install.packages('pacman')
pacman::p_load(tidyverse, lubridate)
rm(list=ls())


# #download data from EDI
# data  <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.199.6&entityid=2b3dc84ae6b12d10bd5485f1c300af13"
# destination <- "./00_Data_files"
# 
# download.file(data,destfile = "./00_Data_files/chemistry.csv", method='libcurl')

#read in sample dates and depths of phyto samples
sample_info <- read_csv("./00_Data_files/EDI_phytos/phytoplankton.csv") %>%
  select(Date, Depth_m) %>%
  distinct()
sample_info$number <- 1:100

cmax <- read_csv("./00_Data_files/FP_DistributionMetrics.csv") %>%
  select(Date, Peak_depth_m)

chem <- read_csv("./00_Data_files/chemistry.csv",
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


#two sampling days missing :-(
#2017-01-19
#2019-02-27

#two sampling days w/ duplicate chem. profiles -- using a.m. one before EM on
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

final <- matrix(NA, nrow = length(nut_dates), ncol = 10)

for (i in 1:length(nut_dates)){
  chem_profile <- subset(chem1, DateTime == nut_dates[i])
  
  phyto_sample <- subset(sample_dates, Date == sample_dates$Date[i])
  
  final[i,1] <- as.character(phyto_sample$Date)
  
  Cmax <- phyto_sample$Depth_m
  
  #pull chem at Cmax
  Cmax_chem <- chem_profile[which.min(abs(chem_profile$Depth_m - Cmax)),]
  
  final[i,2] <- unlist(Cmax_chem[,3])
  final[i,3] <- unlist(Cmax_chem[,4])
  final[i,4] <- unlist(Cmax_chem[,5])
  
  #pull depth of max chem
  #if no clear chem max (i.e. multiple depths with same chem value)
  #then that cell is populated w/ NA
  srp <- subset(chem_profile, SRP_ugL == max(SRP_ugL, na.rm = TRUE))
  if(length(unlist(srp[,2]))>1){
    final[i,5] <- NA
  } else {
    final[i,5] <- unlist(srp[,2])
    final[i,6] <- unlist(srp[,3])
  }
  
  #might have an issue here w/ ammonium buildup in hypo where no light
  din <- subset(chem_profile, DIN_ugL == max(DIN_ugL, na.rm = TRUE))
  if(length(unlist(din[,2]))>1){
    final[i,7] <- NA
  } else {
    final[i,7] <- unlist(din[,2])
    final[i,8] <- unlist(din[,4])
    
  }

  doc <- subset(chem_profile, DOC_mgL == max(DOC_mgL, na.rm = TRUE))
  if(length(unlist(doc[,2]))>1){
    final[i,9] <- NA
  } else {
    final[i,9] <- unlist(doc[,2])
    final[i,10] <- unlist(doc[,5])
  }

}

final <- data.frame(final)
colnames(final) <- c("Date","Cmax_SRP_ugL","Cmax_DIN_ugL","Cmax_DOC_mgL","SRPmax_depth_m","SRPmax_ugL","DINmax_depth_m","DINmax_ugL","DOCmax_depth_m","DOCmax_mgL")

#note this only has 98 rows but should still be able to left_join to final dataframe
write.csv(final, file = "./00_Data_files/chem_vars.csv",row.names = FALSE)

#visualization
chem <- read_csv("./00_Data_files/chem_vars.csv") %>%
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

