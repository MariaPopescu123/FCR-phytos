##Pull an FP sample from a specific site, time, depth
##Author: Mary Lofton
##Date: 27SEP18

#load packages
pacman::p_load(tidyverse, lubridate)

#clear environment
rm(list=ls())

# #download data from EDI
# data  <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.272.4&entityid=e7e3e6e513985a602d9a5f22687d4efc"
# destination <- "./00_Data_files"
# 
# download.file(data,destfile = "./00_Data_files/FP.csv", method='libcurl')

#read in sample dates and depths of phyto samples
sample_info <- read_csv("./00_Data_files/EDI_phytos/phytoplankton.csv") %>%
  select(Date, Depth_m) %>%
  distinct()
sample_info$number <- 1:100

#load fp data to get biomass value
fp <- read_csv("./00_Data_files/FP.csv")%>%
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


#check to make sure don't have duplicate casts on same day
check <- fp_sample %>%
  distinct(Date,Hour)
#the "duplicate" on 07-03 is ok because the cast was taken over 2 hours

#calculate depth distribution characteristics for each profile
dates <- unique(fp_sample$Date)

final <- matrix(NA,nrow = 98, ncol = 5)

for (i in 1:length(dates)){
  #source dependent scripts
  source('./0_Function_library/Script_S2_GND_fit_functions.R')
  
  profile <- fp_sample %>%
    filter(Date == dates[i])
  
  final[i,1]<- as.character(dates[i])
  final[i,2]<- unlist(profile[which.max(profile$TotalConc_ugL),"TotalConc_ugL"])
  final[i,3]<- unlist(profile[which.max(profile$TotalConc_ugL),"Depth_m"])
  final[i,4]<- max(profile$TotalConc_ugL, na.rm = TRUE) - median(profile$TotalConc_ugL, na.rm = TRUE)
  
  #calculate DCM fit
  DCM.fit = fit.GND(profile$Depth_m, profile$TotalConc_ugL)
  DCM.depth  = DCM.fit$par[3]
  DCM.std    = DCM.fit$par[4]
  # calculate DCM top depth, will set to 0 if this is above surface of lake
  DCM.breadth.top = c()
  if (DCM.depth - DCM.std < 0) {
    DCM.breadth.top = 0
  } else if (DCM.depth - DCM.std >= 0) {
    DCM.breadth.top = DCM.depth - DCM.std
  }
  # calculate the DCM bottom depth, will set to max. depth of profile if below
  DCM.breadth.bottom = c() # 
  if (DCM.depth + DCM.std > max(profile$Depth_m)) {
    DCM.breadth.bottom = max(profile$Depth_m)
  } else {
    DCM.breadth.bottom = DCM.depth + DCM.std
  }
  # calculate peak width
  peak.width = DCM.breadth.bottom - DCM.breadth.top
  
  final[i,5] <- peak.width
  
}

final <- data.frame(final)

colnames(final) <- c("Date","Max_biomass_ugL","Peak_depth_m","Peak_magnitude_ugL","Peak_width_m")

#tweaks to match up dates with phytoplankton dates
#need to NOTE in Methods that in cases where profiles were not available from
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

#note there are only 98 rows here but should be ok when left-join to phyto samples
write.csv(final, "./00_Data_files/FP_DistributionMetrics.csv", row.names = FALSE)

#read data back in
final <- read_csv("./00_Data_files/FP_DistributionMetrics.csv")

##PLOTTING
final1 <- final %>%
  mutate(Date = as.Date(Date),
         Max_biomass_ugL = as.numeric(Max_biomass_ugL),
         Peak_depth_m = as.numeric(Peak_depth_m),
         Peak_magnitude_ugL = as.numeric(Peak_magnitude_ugL),
         Peak_width_m = as.numeric(Peak_width_m)) %>%
  gather(Max_biomass_ugL:Peak_width_m, key = "Dist_var",value = "value") %>%
  mutate(Year = year(Date))

vars <- ggplot(data = final1, aes(x = Date, y = value, group = Dist_var, color = Dist_var))+
  facet_wrap(~Dist_var, nrow = 4, scales = "free_y")+
  geom_point(size = 2)+
  geom_line()+
  theme_classic()
vars
ggsave(vars, filename = "./Exploratory_viz/FP_dist_vars.png",
       height = 5, width = 6, units = "in", device = "png")

final1$Dist_var <- factor(final1$Dist_var, levels = c("Max_biomass_ugL","Peak_depth_m","Peak_magnitude_ugL","Peak_width_m"),
                   ordered = TRUE, labels=c(expression(paste("Maximum biomass ","(",mu,g,~L^-1,")")),
                                            expression(paste("Peak depth (m)")), 
                                            expression(paste("Peak magnitude ","(",mu,g,~L^-1,")")),
                                            expression(paste("Peak width (m)"))))
final1$Year <- as.factor(final1$Year)

vars1 <- ggplot(data = final1, aes(x = value, group = Year, color = Year, fill = Year))+
  geom_density(alpha = 0.5)+
  facet_wrap(~Dist_var, nrow = 2, scales = "free",labeller = label_parsed)+
  xlab("")+
  theme_classic()
vars1
ggsave(vars1, filename = "C:/Users/Mary Lofton/Dropbox/Ch_2/Exploratory_viz/FP_dist_vars_hist.png",
       height = 3, width = 10, units = "in", device = "png")

final2 <- final %>%
  mutate(Year = as.factor(year(Date)))

vars2 <- ggplot(data = final2, aes(x = Date, y = Peak_depth_m, group = Year, color = Year, fill = Year))+
  geom_point()+
  geom_line()+
  facet_wrap(vars(Year), scales = "free_x")+
  theme_classic()+
  scale_y_reverse()+
  ylab("Peak depth (m)")
vars2
ggsave(plot = vars2, filename = "C:/Users/Mary Lofton/Dropbox/Ch_2/Exploratory_viz/peakdepth.png",
       device = "png",height = 3, width = 5, units = "in")

vars3 <- ggplot(data = final2, aes(x = Date, y = Peak_magnitude_ugL, group = Year, color = Year, fill = Year))+
  geom_point()+
  geom_line()+
  facet_wrap(vars(Year), scales = "free_x")+
  theme_classic()+
  ylab(expression(paste("Peak magnitude ","(",mu,g,~L^-1,")")))
vars3
ggsave(plot = vars3, filename = "C:/Users/Mary Lofton/Dropbox/Ch_2/Exploratory_viz/peakmagnitude.png",
       device = "png",height = 3, width = 5, units = "in")

vars4 <- ggplot(data = final2, aes(x = Date, y = Peak_width_m, group = Year, color = Year, fill = Year))+
  geom_point()+
  geom_line()+
  facet_wrap(vars(Year), scales = "free_x")+
  theme_classic()+
  ylab(expression(paste("Peak width (m)")))
vars4
ggsave(plot = vars4, filename = "C:/Users/Mary Lofton/Dropbox/Ch_2/Exploratory_viz/peakwidth.png",
       device = "png",height = 3, width = 5, units = "in")

yrs <- unique(final1$Year)
dist_vars <- unique(final1$Dist_var)

png(file = "./Exploratory_viz/FP_pacf.png",width = 16, height = 16,
    units = "cm",res = 300)
par(mfrow = c(4,4), mgp = c(2,0.5,0),mar = c(4,3,3,1))

for (j in 1:length(yrs)){
  for (k in 1:length(dist_vars)){
    
    mydata <- final1 %>%
      filter(Year == yrs[j],
             Dist_var == dist_vars[k])
    
    myacf <- acf(mydata$value, 
                type = "partial",
                plot = FALSE)
    plot(myacf,main = "")
    title(c(yrs[j],dist_vars[k]),line = 1)
    
  }
  
}

dev.off()


