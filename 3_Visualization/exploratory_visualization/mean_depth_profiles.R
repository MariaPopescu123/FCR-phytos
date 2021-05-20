#Fig. 1: Conceptual fig about gradients over depth
#Author: Mary Lofton
#Date: 01MAR21

#load packages
#install.packages('pacman')
pacman::p_load(tidyverse, lubridate, rLakeAnalyzer, data.table, cowplot)
rm(list=ls())

#read in sample dates and depths of phyto samples
sample_info <- read_csv("./00_Data_files/EDI_phytos/phytoplankton.csv") %>%
  select(Date, Depth_m) %>%
  distinct()
sample_info$number <- 1:100

#read in FP data so can match temp profiles to hour FP profiles were taken
replacement_dates <- as.Date(c("2016-07-12","2018-05-24","2019-07-03","2019-07-11","2019-07-18","2019-10-22"))

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
fp_sample <- fp_sample %>%
  distinct()
fp_sample <- fp_sample[-81,]

#read in CTD data and limit to temperature, DO, pH
ctd <- fread("./00_Data_files/CTD.csv")
ctd <- tibble(ctd) %>%
  select(Reservoir, Site, Date, Depth_m, Temp_C, DO_mgL, pH, PAR_umolm2s) %>%
  filter(Reservoir == "FCR" & Site == 50 & date(Date) %in% sample_info$Date) %>%
  mutate(Hour = hour(Date),
         Date = date(Date))

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

#data munging for rLakeAnalyzer
depths = c(0.1, 0.8, 1.6, 2.8, 3.8, 5.0, 6.2, 8.0, 9.0, 9.3)
df.final<-data.frame()

for (i in 1:length(depths)){
  
  ctd_layer<-final %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
  
  # Bind each of the data layers together.
  df.final = bind_rows(df.final, ctd_layer)
  
}

# Re-arrange the data frame by date
ctd_new <- arrange(df.final, Date)

# Round each extracted depth to the nearest 10th. 
ctd_new$Depth_m <- round(as.numeric(ctd_new$Depth_m), digits = 0.5) 

#fix a couple of wonky depths
ctd_new$Depth_m[59:60] <- 9.0
ctd_new <- ctd_new[-c(39:40,60,298:300),]

#cleaning casts where measurements don't start until deeper in water column
ctd_new <- ctd_new %>%
  filter(!(Depth_m >0.3 & Depth_m < 0.8))

ctd_check <- ctd_new %>%
  filter(Depth_m < 0.8)
hist(ctd_check$Depth_m)

ctd_new <- ctd_new %>%
  mutate(Depth_m = ifelse(Depth_m < 0.8, 0.1, 
                          ifelse(Depth_m == 9.2, 9.3, Depth_m))) 

ctd_new <- ctd_new %>%
  mutate(Year = year(Date))

my.cs.data <- read_csv("./2_Data_analysis/CS_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_DO_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") %>%
  select(-MonthDay) %>%
  filter(!is.na(BV_TOTAL))

ctd_new <- ctd_new %>%
  filter(Date %in% my.cs.data$Date) %>%
  mutate(Year = as.character(Year))

mean_td <- ctd_new %>%
  group_by(Depth_m)%>%
  summarize(mean_temp = mean(Temp_C, na.rm = TRUE),
            q2.5 = quantile(Temp_C,0.025, na.rm = TRUE),
            q97.5 = quantile(Temp_C, 0.975, na.rm = TRUE)) %>%
  gather(mean_temp:q97.5, key = "summary.stat", value = "temp_c")

p.ctd <- ggplot(data = mean_td, aes(x = temp_c, y = Depth_m, group = summary.stat, linetype = summary.stat))+
  geom_path(size = 1.5, col = "cornflowerblue")+
  scale_y_reverse()+
  theme_classic()+
  xlab("Water temperature (degrees Celsius)")+
  ylab("Depth (m)")
p.ctd

p.ctd1 <- ggplot(data = subset(mean_td, summary.stat == "mean_temp"), aes(x = temp_c, y = Depth_m, group = summary.stat, linetype = summary.stat))+
  geom_path(size = 1.5, col = "cornflowerblue")+
  scale_y_reverse()+
  theme_classic()+
  xlab("Water temperature (degrees Celsius)")+
  ylab("Depth (m)")
p.ctd1
ggsave(p.ctd1, filename = "C:/Users/Mary Lofton/Documents/Presentations/mean_temp_FCR.png",
       device = "png", height = 2.5, width = 3.5, units = "in")

mean_td_temp <- ctd_new[1:10,]
mean_td_temp$Date <- as.Date("2021-01-01")
mean_td_temp$Depth_m <- mean_td$Depth_m
mean_td_temp$Temp_C <- mean_td$mean_temp
mean_td_temp$Year <- "mean"
mean_td_temp[,6:8] <- NA

ctd_new <- bind_rows(ctd_new, mean_td_temp)

p <- ggplot(data = subset(ctd_new, ctd_new$Date == "2016-07-18"), aes(x = Temp_C, y = Depth_m, group = Date, linetype = as.factor(Date)))+
  geom_path(size = 1.5, col = "cornflowerblue")+
  scale_y_reverse()+
  theme_classic() +
  ylab("Depth (m)")+
  xlab("Water temperature (degrees Celsius)")
p 
ggsave(p, filename = "C:/Users/Mary Lofton/Documents/Presentations/td_pre_FCR.png",
       device = "png", height = 2.5, width = 3.5, units = "in")

p <- ggplot(data = subset(ctd_new, ctd_new$Date == "2016-07-18"| ctd_new$Date == "2016-08-01"), aes(x = Temp_C, y = Depth_m, group = Date, linetype = as.factor(Date)))+
  geom_path(size = 1.5, col = "cornflowerblue")+
  scale_y_reverse()+
  theme_classic() +
  ylab("Depth (m)")+
  xlab("Water temperature (degrees Celsius)")
p 
ggsave(p, filename = "C:/Users/Mary Lofton/Documents/Presentations/td_lowering_FCR.png",
       device = "png", height = 2.5, width = 3.5, units = "in")


##getting layers for FP profiles
#data munging for rLakeAnalyzer
depths = c(0.1, 0.8, 1.6, 2.8, 3.8, 5.0, 6.2, 8.0, 9.0, 9.3)
df.final.fp<-data.frame()

for (i in 1:length(depths)){
  
  fp_layer<-fp_sample %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
  
  # Bind each of the data layers together.
  df.final.fp = bind_rows(df.final.fp, fp_layer)
  
}

# Re-arrange the data frame by date
fp_new <- arrange(df.final.fp, Date)

# Round each extracted depth to the nearest 10th. 
fp_new$Depth_m <- round(as.numeric(fp_new$Depth_m), digits = 0.5) 

#cleaning casts where measurements don't start until deeper in water column

fp_new <- fp_new %>%
  mutate(Depth_m = ifelse(Depth_m < 0.6, 0.5, 
                          ifelse(Depth_m >= 0.6 & Depth_m <=1, 0.8,
                                 ifelse(Depth_m >= 1.4 & Depth_m <= 1.8,1.6,
                                        ifelse(Depth_m >= 2.6 & Depth_m <= 3.0,2.8,
                                               ifelse(Depth_m >= 3.6 & Depth_m <= 4.0,3.8,
                                                      ifelse(Depth_m >=4.8 & Depth_m <= 5.2, 5.0,
                                                             ifelse(Depth_m >= 6.0 & Depth_m <= 6.4, 6.2,
                                                                    ifelse(Depth_m >= 7.8 & Depth_m <= 8.2, 8.0,
                                                                           ifelse(Depth_m >= 8.8 & Depth_m <= 9.2, 9.0,
                                                                                  ifelse(Depth_m > 9.2 & Depth_m <= 9.5, 9.3, NA)))))))))))  %>%
  mutate(Year = year(Date)) %>%
  filter(Date %in% my.cs.data$Date)

fp_new <- fp_new %>%
  filter(!is.na(Depth_m))

p <- ggplot(data = subset(fp_new, fp_new$Year == "2016"), aes(x = TotalConc_ugL, y = Depth_m, group = Date))+
  geom_path()+
  scale_y_reverse()+
  theme_classic()
p 

mean_fp <- fp_new %>%
  group_by(Depth_m)%>%
  summarize(mean_conc = mean(TotalConc_ugL, na.rm = TRUE),
            q2.5 = quantile(TotalConc_ugL,0.025, na.rm = TRUE),
            q97.5 = quantile(TotalConc_ugL, 0.975, na.rm = TRUE)) %>%
  gather(mean_conc:q97.5, key = "summary.stat", value = "ugL")

p.fp <- ggplot(data = mean_fp, aes(x = ugL, y = Depth_m, group = summary.stat, linetype = summary.stat))+
  geom_path(size = 1.5)+
  scale_y_reverse()+
  theme_classic()
p.fp


##chemistry
#read in sample dates and depths of phyto samples
chem <- read_csv("./00_Data_files/chemistry.csv",
                 col_types = cols(DIC_mgL = col_double(),
                                  DC_mgL = col_double(),
                                  DN_mgL = col_double(),
                                  Flag_DIC = col_double(),
                                  Flag_DC = col_double(),
                                  Flag_DN = col_double())) %>%
  filter(Reservoir == "FCR" & Site == 50 & date(DateTime) %in% my.cs.data$Date) %>%
  filter(Depth_m %in% c(0.1, 0.8, 1.6, 2.8, 3.8, 5.0, 6.2, 8.0, 9.0, 9.3))

#check for sampling day alignment btwn. phytos and chem
chem_dates <- data.frame(unique(chem$DateTime))
chem_dates1 <- data.frame(date(unique(chem$DateTime)))
colnames(chem_dates1) <- "Date"
chem_dates1$chem_samples <- TRUE
check <- left_join(cmax, chem_dates1, by = "Date")

#subset to correct days and vars of interest
chem1 <- chem %>%
  filter(!DateTime  == chem_dates$unique.chem.DateTime.[7] & !DateTime == chem_dates$unique.chem.DateTime.[15] & !DateTime == chem_dates$unique.chem.DateTime.[68]) %>%
  mutate(DIN_ugL = NO3NO2_ugL + NH4_ugL) %>%
  select(DateTime, Depth_m, SRP_ugL, DIN_ugL, DOC_mgL) %>%
  mutate(Date = date(DateTime), Year = year(DateTime)) %>%
  select(-DateTime)

p <- ggplot(data = chem1, aes(x = DIN_ugL, y = Depth_m, group = Date, color = as.factor(Year)))+
  geom_path()+
  scale_y_reverse()+
  theme_classic()
p 

p <- ggplot(data = subset(chem1, chem1$Year == "2016"), aes(x = SRP_ugL, y = Depth_m, group = Date))+
  geom_path(color = "darkorange")+
  scale_y_reverse()+
  theme_classic()+
  ylab("Depth (m)")+
  xlab("Soluble reactive phosphorus (ug/L)")
p
ggsave(p, filename = "C:/Users/Mary Lofton/Documents/Presentations/SRP_profiles_2016_FCR.png",
       device = "png", height = 2.5, width = 3.5, units = "in")

p 

p <- ggplot(data = chem1, aes(x = DOC_mgL, y = Depth_m, group = Date, color = as.factor(Year)))+
  geom_path()+
  scale_y_reverse()+
  theme_classic()
p 

mean_chem <- chem1 %>%
  group_by(Depth_m)%>%
  summarize(mean_DIN = mean(DIN_ugL, na.rm = TRUE),
            q2.5_DIN = quantile(DIN_ugL,0.025, na.rm = TRUE),
            q97.5_DIN = quantile(DIN_ugL, 0.975, na.rm = TRUE),
            mean_SRP = mean(SRP_ugL, na.rm = TRUE),
            q2.5_SRP = quantile(SRP_ugL,0.025, na.rm = TRUE),
            q97.5_SRP = quantile(SRP_ugL, 0.975, na.rm = TRUE),
            mean_DOC = mean(DOC_mgL, na.rm = TRUE),
            q2.5_DOC = quantile(DOC_mgL,0.025, na.rm = TRUE),
            q97.5_DOC = quantile(DOC_mgL, 0.975, na.rm = TRUE)) %>%
  gather(mean_DIN:q97.5_DOC, key = "summary.stat", value = "conc")
plotchem <- mean_chem[grep("mean",mean_chem$summary.stat),]

p.nuts <- ggplot(data = plotchem, aes(x = conc, y = Depth_m, group = summary.stat, color = summary.stat))+
  geom_path(size = 1.5, col = "royalblue")+
  facet_wrap(vars(summary.stat), scales = "free_x")+
  scale_y_reverse()+
  theme_classic()+
  ylab("Depth (m)")+
  xlab("Nutrient concentration")
p.nuts

DIN <- mean_chem %>%
  filter(summary.stat %in% c("mean_DIN","q2.5_DIN","q97.5_DIN"))
p.DIN <- ggplot(data = DIN, aes(x = conc, y = Depth_m, group = summary.stat, linetype = summary.stat))+
  geom_path(size = 1.5)+
  scale_y_reverse()+
  theme_classic()
p.DIN

p.DIN1 <- ggplot(data = subset(DIN, DIN$summary.stat == "mean_DIN"), aes(x = conc, y = Depth_m, group = summary.stat, linetype = summary.stat))+
  geom_path(size = 1.5, col = "royalblue")+
  scale_y_reverse()+
  theme_classic()+
  ylab("Depth (m)")+
  xlab("Dissolved inorganic nitrogen (ug/L)")
p.DIN1
ggsave(p.DIN1, filename = "C:/Users/Mary Lofton/Documents/Presentations/mean_DIN_FCR.png",
       device = "png", height = 2.5, width = 3.5, units = "in")


SRP <- mean_chem %>%
  filter(summary.stat %in% c("mean_SRP","q2.5_SRP","q97.5_SRP"))
p.SRP <- ggplot(data = SRP, aes(x = conc, y = Depth_m, group = summary.stat, linetype = summary.stat))+
  geom_path(size = 1.5)+
  scale_y_reverse()+
  theme_classic()
p.SRP

p.SRP1 <- ggplot(data = subset(SRP, SRP$summary.stat == "mean_SRP"), aes(x = conc, y = Depth_m, group = summary.stat, linetype = summary.stat))+
  geom_path(size = 1.5, col = "darkorange")+
  scale_y_reverse()+
  theme_classic()+
  ylab("Depth (m)")+
  xlab("Soluble reactive phosphorus (ug/L)")
p.SRP1
ggsave(p.SRP1, filename = "C:/Users/Mary Lofton/Documents/Presentations/mean_SRP_FCR.png",
       device = "png", height = 2.5, width = 3.5, units = "in")


DOC <- mean_chem %>%
  filter(summary.stat %in% c("mean_DOC","q2.5_DOC","q97.5_DOC"))
p.DOC <- ggplot(data = DOC, aes(x = conc, y = Depth_m, group = summary.stat, linetype = summary.stat))+
  geom_path(size = 1.5)+
  scale_y_reverse()+
  theme_classic()
p.DOC

p.DOC1 <- ggplot(data = subset(DOC, DOC$summary.stat == "mean_DOC"), aes(x = conc, y = Depth_m, group = summary.stat, linetype = summary.stat))+
  geom_path(size = 1.5, col = "black")+
  scale_y_reverse()+
  theme_classic()+
  ylab("Depth (m)")+
  xlab("Dissolved organic carbon (mg/L)")
p.DOC1
ggsave(p.DOC1, filename = "C:/Users/Mary Lofton/Documents/Presentations/mean_DOC_FCR.png",
       device = "png", height = 2.5, width = 3.5, units = "in")




##PAR data
#data munging for rLakeAnalyzer
depths = c(-0.1, 0.1, 0.8, 1.6, 2.8, 3.8, 5.0, 6.2, 8.0, 9.0, 9.3)
df.final.par<-data.frame()

for (i in 1:length(depths)){
  
  par_layer<-final %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
  
  # Bind each of the data layers together.
  df.final.par = bind_rows(df.final.par, par_layer)
  
}

# Re-arrange the data frame by date
par_new <- arrange(df.final.par, Date) %>%
  filter(!is.na(PAR_umolm2s))

datez <- unique(par_new$Date)

par_newer <- data.frame()

for (i in 1:length(datez)){
  
  today <- subset(par_new, par_new$Date == datez[i])
  
  today$PAR_perc <- today$PAR_umolm2s/today$PAR_umolm2s[1]*100
  
  # Bind each of the data layers together.
  par_newer = bind_rows(par_newer, today)
  
}

# Round each extracted depth to the nearest 10th. 
par_newer$Depth_m <- round(as.numeric(par_newer$Depth_m), digits = 0.5) 
unique(par_newer$Depth_m)

# #fix a couple of wonky depths
# ctd_new$Depth_m[59:60] <- 9.0
# ctd_new <- ctd_new[-c(39:40,60,298:300),]

#cleaning casts where measurements don't start until deeper in water column
par_newer <- par_newer  %>%
  mutate(Depth_m = ifelse(Depth_m == 6.3, 6.2, 
                          ifelse(Depth_m == 9.2, 9.3, Depth_m))) 

par_newer <- par_newer %>%
  mutate(Year = year(Date)) %>%
  filter(Date %in% my.cs.data$Date) 

check <- par_newer %>%
  filter(PAR_perc > 100)

par_newer <- par_newer[-c(12:22,133:143),]

p <- ggplot(data = par_newer, aes(x = PAR_perc, y = Depth_m, group = Date, color = as.factor(Year)))+
  geom_path()+
  scale_y_reverse()+
  theme_classic()
p 

mean_par <- par_newer %>%
  group_by(Depth_m)%>%
  summarize(mean_par = mean(PAR_perc, na.rm = TRUE),
            q2.5 = quantile(PAR_perc,0.025, na.rm = TRUE),
            q97.5 = quantile(PAR_perc, 0.975, na.rm = TRUE)) %>%
  gather(mean_par:q97.5, key = "summary.stat", value = "par")

p.par <- ggplot(data = mean_par, aes(x = par, y = Depth_m, group = summary.stat, linetype = summary.stat))+
  geom_path(size = 1.5)+
  scale_y_reverse()+
  theme_classic()
p.par

p.par1 <- ggplot(data = subset(mean_par, summary.stat == "mean_par"), aes(x = par, y = Depth_m, group = summary.stat, linetype = summary.stat))+
  geom_path(size = 1.5, col = "goldenrod")+
  scale_y_reverse()+
  theme_classic()+
  xlab("% surface light")+
  ylab("Depth (m)")
p.par1

ggsave(p.par1, filename = "C:/Users/Mary Lofton/Documents/Presentations/mean_par_FCR.png",
       device = "png", height = 2.5, width = 3.5, units = "in")

fig1 <- plot_grid(p.ctd, p.par, p.DIN, p.SRP, p.DOC, p.fp, nrow = 2, ncol = 3, 
                  align = "hv")
fig1
ggsave(fig1, filename = "./3_Visualization/Fig1.tif",height = 6, width = 12,
       units = "in", dpi = 300, dev = "tiff")
