#1B_Chemistry
#Author: Mary Lofton
#Date: 23JUN20

#pull in chemistry data from EDI and subset to relevant dates/drivers

#load packages
#install.packages('pacman')
pacman::p_load(tidyverse, lubridate)
rm(list=ls())


#download data from EDI
data  <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.199.6&entityid=2b3dc84ae6b12d10bd5485f1c300af13"
destination <- "./00_Data_files"

download.file(data,destfile = "./00_Data_files/chemistry.csv", method='libcurl')

#read in sample dates and depths of phyto samples
sample_info <- read_csv("./00_Data_files/EDI_phytos/phytoplankton.csv") %>%
  select(Date, Depth_m) %>%
  distinct()
sample_info$number <- 1:100

chem <- read_csv("./00_Data_files/chemistry.csv",
                 col_types = cols(DIC_mgL = col_double(),
                                  DC_mgL = col_double(),
                                  DN_mgL = col_double(),
                                  Flag_DIC = col_double(),
                                  Flag_DC = col_double(),
                                  Flag_DN = col_double())) %>%
  filter(Reservoir == "FCR" & Site == 50 & date(DateTime) %in% sample_info$Date) 

#check for sampling day alignment btwn. phytos and chem
chem_dates <- data.frame(unique(chem$DateTime))
chem_dates$number <- 1:100
check <- left_join(sample_info, chem_dates, by = "number")

#two sampling days missing :-(
#2017-01-19
#2019-02-27

#two sampling days w/ duplicate chem. profiles -- using a.m. one before EM on
#eliminating
#2016-05-30 16:00:00
#2016-07-25 00:00:00

#subset to correct days and vars of interest
chem1 <- chem %>%
  filter(!DateTime  == chem_dates$unique.chem.DateTime.[7] & !DateTime == chem_dates$unique.chem.DateTime.[15]) %>%
  mutate(DIN_ugL = NO3NO2_ugL + NH4_ugL) %>%
  select(DateTime, Depth_m, SRP_ugL, DIN_ugL, DOC_mgL)

#pull nuts at depth of Cmax and depth of max. nuts
nut_dates <- unique(chem1$DateTime)
sample_dates <- sample_info[-c(24,68),] %>% #get rid of dates not in chem data
  select(Date, Depth_m)

final <- matrix(NA, nrow = length(nut_dates), ncol = 7)

for (i in 1:length(nut_dates)){
  chem_profile <- subset(chem1, DateTime == nut_dates[i])
  
  phyto_sample <- subset(sample_dates, Date == sample_dates$Date[i])
  
  final[i,1] <- as.character(phyto_sample$Date)
  
  Cmax <- phyto_sample$Depth_m
  
  #pull chem at Cmax
  Cmax_chem <- chem1[which.min(abs(chem1$Depth_m - Cmax)),]
  
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
  }
  
  #might have an issue here w/ ammonium buildup in hypo where no light
  din <- subset(chem_profile, DIN_ugL == max(DIN_ugL, na.rm = TRUE))
  if(length(unlist(din[,2]))>1){
    final[i,6] <- NA
  } else {
    final[i,6] <- unlist(din[,2])
  }

  doc <- subset(chem_profile, DOC_mgL == max(DOC_mgL, na.rm = TRUE))
  if(length(unlist(doc[,2]))>1){
    final[i,7] <- NA
  } else {
    final[i,7] <- unlist(doc[,2])
  }

}

final <- data.frame(final)
colnames(final) <- c("Date","Cmax_SRP_ugL","Cmax_DIN_ugL","Cmax_DOC_ugL","SRPmax_depth_m","DINmax_depth_m","DOCmax_depth_m")

#note this only has 98 rows but should still be able to left_join to final dataframe
write.csv(final, file = "./00_Data_files/chem_vars.csv",row.names = FALSE)
