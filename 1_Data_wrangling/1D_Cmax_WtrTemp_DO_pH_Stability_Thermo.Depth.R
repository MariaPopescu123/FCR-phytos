#1D_WtrTemp_DO_pH_Stability_Thermo.Depth
#Author: Mary Lofton
#Date: 15SEP20

#pull in wtr temp data from EDI and subset to relevant dates/drivers

####SET-UP####
#load packages
#install.packages('pacman')
pacman::p_load(tidyverse, lubridate, rLakeAnalyzer, data.table)
rm(list=ls())

#read in sample dates and depths of phyto samples
sample_info <- read_csv("./0_Data_files/phytoplankton.csv") %>%
  select(Date, Depth_m) %>%
  distinct()
sample_info$number <- 1:100

#read in FP data so can match temp profiles to hour FP profiles were taken
replacement_dates <- as.Date(c("2016-07-12","2018-05-24","2019-07-03","2019-07-11","2019-07-18","2019-10-22"))

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
fp_sample <- fp_sample %>%
  distinct(Date,Hour)
fp_sample <- fp_sample[-81,]

#read in met data for lake number and Wedderburn number
#read in and wrangle data
met <- fread("./0_Data_files/Met_final_2015_2020.csv", fill = TRUE, blank.lines.skip = FALSE,select = c("DateTime",
                                                                                                        "Rain_Total_mm",                                                        
                                                                                                        "WindSpeed_Average_m_s",
                                                                                                        "Flag_Rain_Total_mm",
                                                                                                        "Flag_WindSpeed_Average_m_s"))%>%
  mutate(Date = date(DateTime)) %>%
  filter(Flag_Rain_Total_mm == 0 & Flag_WindSpeed_Average_m_s == 0) %>%
  group_by(Date)  %>%
  summarize(daily_precip = sum(Rain_Total_mm, na.rm = TRUE),
            daily_wind = mean(WindSpeed_Average_m_s, na.rm = TRUE))

met2 <- met %>%
  filter(Date %in% sample_info$Date)

wnd <- met2 %>%
  rename(datetime = Date, wnd = daily_wind) %>%
  select(datetime, wnd)
write.table(wnd, "./0_Data_files/wnd.wnd",sep='\t',row.names = FALSE)
wnd <- load.ts("./0_Data_files/wnd.wnd")

####CTD####
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
  select(Date, Depth_m, Temp_C)

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
                          ifelse(Depth_m == 9.2, 9.3, Depth_m))) %>%
  spread(key = "Depth_m", value = "Temp_C") %>%
  rename(datetime = Date)

colnames(ctd_new)[2:11] <- paste("wtr", colnames(ctd_new)[2:11], sep = "_")

wtr_ctd <- ctd_new %>%
  mutate(datetime = format(datetime, format = "%Y-%m-%d %H:%M:%S"))
write.table(wtr_ctd, "./0_Data_files/wtr_ctd.wtr",sep='\t',row.names = FALSE)
wtr_ctd <- load.ts("./0_Data_files/wtr_ctd.wtr")

bathy <- load.bathy("./0_Data_files/FCR.bth")

#calculate Schmidt stability
ss_ctd <- ts.schmidt.stability(wtr = wtr_ctd, bathy = bathy, na.rm = TRUE)
ss_ctd[ss_ctd<0] <- 0

ss_ctd1 <- ss_ctd %>%
  mutate(Year = year(datetime))

ggplot(data = ss_ctd1, aes(x = datetime, y = schmidt.stability))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()

#calculate thermocline depth
td_ctd <- ts.thermo.depth(wtr = wtr_ctd, na.rm = TRUE, seasonal = TRUE) 

td_ctd1 <- td_ctd %>%
  mutate(Year = year(datetime))

ggplot(data = td_ctd1, aes(x = datetime, y = thermo.depth))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  scale_y_reverse()+
  theme_classic()

#calculate buoyancy frequency
bf_ctd <- ts.buoyancy.freq(wtr = wtr_ctd, na.rm = TRUE, seasonal = TRUE) 

bf_ctd1 <- bf_ctd %>%
  mutate(Year = year(datetime))

ggplot(data = bf_ctd1, aes(x = datetime, y = n2))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  scale_y_reverse()+
  theme_classic()

#calculate lake number
wnd_ctd <- wnd %>%
  filter(datetime %in% unique(wtr_ctd$datetime)) 

#define internal functions for lake number
drop.datetime = function(data, error=FALSE){
  datetime.pattern = "(datetime|timestamp|time|date)"
  
  header = names(data)
  dt_indx = grep(datetime.pattern, header, ignore.case=TRUE)
  
  if(length(dt_indx) < 1){
    if(error){
      stop('Unable to find a datetime column. Datetime column was supplied.')
    }else{
      warning('Unable to find a datetime column. Assuming no datetime column was supplied.')
      return(data)
    }
    
  }else if(length(dt_indx) > 1){
    stop('datetime column ambiguity. You can only have one column of datetime.')
  }
  
  return(data[,-dt_indx, drop=FALSE])
}

get.datetime = function(data, error=FALSE){
  datetime.pattern = "(datetime|timestamp|time|date)"
  
  header = names(data)
  dt_indx = grep(datetime.pattern, header, ignore.case=TRUE)
  
  if(length(dt_indx) < 1){
    if(error){
      stop('Unable to find a datetime column.')
    }else{
      warning('Unable to find a datetime column, attempting to ignore.')
      return(NULL)
    }
  }else if(length(dt_indx) > 1){
    stop('datetime column ambiguity. You can only have one column of datetime.')
  }
  
  return(data[,dt_indx])
}

#define lake number function with modifications to handle missing data
ts.lake.number <- function(wtr, wnd, wnd.height, bathy, seasonal=TRUE){
  
  # Make sure data frames match by date/time. 
  all.data = merge(wtr, wnd, by='datetime')
  
  cols = ncol(all.data)
  wtr = all.data[,-cols]
  wnd = all.data[,c(1, cols)]
  
  n = nrow(wtr)
  l.n = rep(NA, n)
  
  wtr.mat = as.matrix(wtr[,-1])
  dimnames(wtr.mat) <- NULL
  
  for(i in 1:n){
    
    depths = get.offsets(wtr)
    
    if(any(is.na(wtr.mat[i,])) || is.na(wnd[i,2])){
      
      wtr.now <- wtr.mat[i,!is.na(wtr.mat[i,])]
      depths.now <- depths[!is.na(wtr.mat[i,])]
      
    } else {
      
      wtr.now <- wtr.mat[i,]
      depths.now <- depths
      
    }
    
    m.d = meta.depths(wtr.now, depths.now, seasonal=seasonal)
    if(any(is.na(m.d))){
      next
    }
    
    epi.dens = layer.density(0.1, m.d[1], wtr.now, depths.now, bathy$areas, bathy$depths)
    hypo.dens = layer.density(m.d[2], max(depths.now), wtr.now, depths.now, bathy$areas, bathy$depths)
    
    
    uS = uStar(wnd[i,2], wnd.height, epi.dens)
    
    St = schmidt.stability(wtr.now, depths.now, bathy$areas, bathy$depths)
    
    #thermo.depth <- function(wtr, depths, Smin = 0.1){\
    l.n[i] = lake.number(bathy$areas, bathy$depths, uS, St, m.d[1], m.d[2], hypo.dens)
  }
  
  output = data.frame(datetime=get.datetime(wtr), lake.number=l.n)
  
  return(output)
}


ln_ctd <- ts.lake.number(wtr = wtr_ctd, wnd = wnd_ctd, wnd.height = 3, bathy = bathy, seasonal = TRUE) 

ln_ctd1 <- ln_ctd %>%
  mutate(Year = year(datetime))

ggplot(data = ln_ctd1, aes(x = datetime, y = lake.number))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()

#calculate Wedderburn number

#define edited function to handle missing data
ts.wedderburn.number <- function(wtr, wnd, wnd.height, bathy, Ao, seasonal=TRUE){
  
  # Make sure data frames match by date/time. 
  all.data = merge(wtr, wnd, by='datetime')
  
  cols = ncol(all.data)
  wtr = all.data[,-cols]
  wnd = all.data[,c(1, cols)]
  
  n = nrow(wtr)
  w.n = rep(NA, n)
  
  wtr.mat = as.matrix(wtr[,-1])
  dimnames(wtr.mat) <- NULL
  
  for(i in 1:n){
    
    depths = get.offsets(wtr)
    
    #check we have all the data necessary
    if(any(is.na(wtr.mat[i,])) || is.na(wnd[i,2])){
      
      wtr.now <- wtr.mat[i,!is.na(wtr.mat[i,])]
      depths.now <- depths[!is.na(wtr.mat[i,])]
      
    } else {
      
      wtr.now <- wtr.mat[i,]
      depths.now <- depths
      
    }
    
    m.d = meta.depths(wtr.now, depths.now, seasonal=seasonal)
    if(any(is.na(m.d))){
      next
    }
    
    #Need epi and hypo density for wedderburn.number calc
    epi.dens = layer.density(0.1, m.d[1], wtr.now, depths.now, bathy$areas, bathy$depths)
    hypo.dens = layer.density(m.d[2], max(depths.now), wtr.now, depths.now, bathy$areas, bathy$depths)
    
    
    uS = uStar(wnd[i,2], wnd.height, epi.dens)
    
    #St = schmidt.stability(wtr.mat[i,], depths, bathy$areas, bathy$depths)
    
    #thermo.depth <- function(wtr, depths, Smin = 0.1){\
    #l.n[i] = lake.number(bathy$areas, bathy$depths, uS, St, m.d[1], m.d[2], hypo.dens)
    
    w.n[i] = wedderburn.number(hypo.dens - epi.dens, m.d[1], uS, Ao, hypo.dens)
    
  }
  
  output = data.frame(datetime=get.datetime(wtr), wedderburn.number=w.n)
  
  return(output)
}

wn_ctd <- ts.wedderburn.number(wtr = wtr_ctd, wnd = wnd_ctd, wnd.height = 3, bathy = bathy, Ao = 119000, seasonal = TRUE) 

wn_ctd1 <- wn_ctd %>%
  mutate(Year = year(datetime))

ggplot(data = wn_ctd1, aes(x = datetime, y = wedderburn.number))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()

#calculate water density

#calculate metalimnion depths
meta_ctd <- ts.meta.depths(wtr = wtr_ctd, na.rm = TRUE, seasonal = TRUE)
mean(meta_ctd$top, na.rm = TRUE)
mean(meta_ctd$bottom, na.rm = TRUE)

#find water temperature at depth of Cmax
cmax <- read_csv("./0_Data_files/FP_DistributionMetrics.csv")

depths = cmax$Peak_depth_m
dates  = cmax$Date

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

df.final<-data.frame()

for (i in 1:length(dates)){
  
  profile <- final %>%
    filter(Date == dates[i])
  
  ctd_layer<-profile %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
  
  # Bind each of the data layers together.
  df.final = bind_rows(df.final, ctd_layer)
  
}

grab_depths = sample_info$Depth_m
grab_dates = sample_info$Date

final.grab <- ctd[0,]

ctd_dates <- unique(ctd$Date)
for (i in 1:length(ctd_dates)){
  ctd_sample <- subset(ctd, ctd$Date == ctd_dates[i])
  fp_hour <- subset(fp_sample, fp_sample$Date == ctd_dates[i])
  ctd_profile <- ctd_sample[ctd_sample[, "Hour"] == closest(ctd_sample$Hour,fp_hour$Hour[1]),]
  final.grab <- bind_rows(final.grab, ctd_profile)
}

final.grab <- final.grab %>%
  select(Date, Depth_m, Temp_C)

df.final.grab<-data.frame()

for (i in 1:length(grab_dates)){
  
  profile <- final.grab %>%
    filter(Date == grab_dates[i])
  
  ctd_layer<-profile %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - grab_depths[i])))
  
  # Bind each of the data layers together.
  df.final.grab = bind_rows(df.final.grab, ctd_layer)
  
}
colnames(df.final.grab) <-c("Date","Grab_T_Depth_m","Grab_Depth_Temp_C")

check <- left_join(df.final.grab,sample_info,by = "Date")


colnames(ss_ctd)[1] <- "Date"
colnames(td_ctd)[1] <- "Date"
colnames(bf_ctd)[1] <- "Date"
colnames(ln_ctd)[1] <- "Date"
colnames(wn_ctd)[1] <- "Date"


ctd_all <- left_join(df.final, ss_ctd, by = "Date")
ctd_all1 <- left_join(ctd_all, td_ctd, by = "Date")
ctd_all2 <- left_join(ctd_all1, bf_ctd, by = "Date")
ctd_all3 <- left_join(ctd_all2, ln_ctd, by = "Date")
ctd_all4 <- left_join(ctd_all3, wn_ctd, by = "Date")
ctd_final <- left_join(ctd_all4,df.final.grab, by = "Date")


####YSI####
#read in YSI data and limit to temperature
ysi <- read_csv("./0_Data_files/YSI.csv") %>%
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
  select(Date, Depth_m, Temp_C)

#data munging for rLakeAnalyzer
depths = c(0.1, 0.8, 1.6, 2.8, 3.8, 5.0, 6.2, 8.0, 9.0, 9.3)
df.final<-data.frame()

for (i in 1:length(depths)){
  
  ysi_layer<-final %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
  
  # Bind each of the data layers together.
  df.final = bind_rows(df.final, ysi_layer)
  
}

# Re-arrange the data frame by date
ysi_new <- arrange(df.final, Date)

# Round each extracted depth to the nearest 10th. 
ysi_new$Depth_m <- round(as.numeric(ysi_new$Depth_m), digits = 0.5) 

# Adjust days when we measured at 0.7 and 2.9
ysi_new <- ysi_new %>%
  mutate(Depth_m = ifelse(Depth_m == 2.9, 2.8, 
                          ifelse(Depth_m == 0.7, 0.8, Depth_m))) 

ysi_new <- ysi_new[complete.cases(ysi_new), ]

ysi_new <- ysi_new %>%
  distinct() %>%
  filter(!Date == "2019-06-10") %>%
  spread(key = "Depth_m", value = "Temp_C") %>%
  filter(!Date == "2019-06-10") %>%
  rename(datetime = Date)

ysi_new <- ysi_new[-c(36, 64:67),]

colnames(ysi_new)[2:12] <- paste("wtr", colnames(ysi_new)[2:12], sep = "_")

wtr_ysi <- ysi_new %>%
  mutate(datetime = format(datetime, format = "%Y-%m-%d %H:%M:%S"))
wtr_ysi[37,2] <- NA
write.table(wtr_ysi, "./0_Data_files/wtr_ysi.wtr",sep='\t',row.names = FALSE)
wtr_ysi <- load.ts("./0_Data_files/wtr_ysi.wtr")
bathy <- load.bathy("./0_Data_files/FCR.bth")

#calculate Schmidt stability
ss_ysi <- ts.schmidt.stability(wtr = wtr_ysi, bathy = bathy, na.rm = TRUE)
ss_ysi[ss_ysi<0] <- 0

ss_ysi1 <- ss_ysi %>%
  mutate(Year = year(datetime))

ggplot(data = ss_ysi1, aes(x = datetime, y = schmidt.stability))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()

#calculate thermocline depth
td_ysi <- ts.thermo.depth(wtr = wtr_ysi, na.rm = TRUE, seasonal = TRUE) 

td_ysi1 <- td_ysi %>%
  mutate(Year = year(datetime))

ggplot(data = td_ysi1, aes(x = datetime, y = thermo.depth))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  scale_y_reverse()+
  theme_classic()

#calculate buoyancy frequency
bf_ysi <- ts.buoyancy.freq(wtr = wtr_ysi, na.rm = TRUE, seasonal = TRUE) 

bf_ysi1 <- bf_ysi %>%
  mutate(Year = year(datetime))

ggplot(data = bf_ysi1, aes(x = datetime, y = n2))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  scale_y_reverse()+
  theme_classic()

#calculate lake number
wnd_ysi <- wnd %>%
  filter(datetime %in% unique(wtr_ysi$datetime)) 

ln_ysi <- ts.lake.number(wtr = wtr_ysi, wnd = wnd_ysi, wnd.height = 3, bathy = bathy, seasonal = TRUE) 

ln_ysi1 <- ln_ysi %>%
  mutate(Year = year(datetime))

ggplot(data = ln_ysi1, aes(x = datetime, y = lake.number))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()

#calculate wedderburn number
wn_ysi <- ts.wedderburn.number(wtr = wtr_ysi, wnd = wnd_ysi, wnd.height = 3, bathy = bathy, Ao = 119000, seasonal = TRUE) 

wn_ysi1 <- wn_ysi %>%
  mutate(Year = year(datetime))

ggplot(data = wn_ysi1, aes(x = datetime, y = wedderburn.number))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()

#find water temperature at depth of Cmax
cmax <- read_csv("./0_Data_files/FP_DistributionMetrics.csv")

depths = cmax$Peak_depth_m
dates  = cmax$Date

final <- ysi[0,]

ysi_dates <- unique(ysi$Date)
for (i in 1:length(ysi_dates)){
  ysi_sample <- subset(ysi, ysi$Date == ysi_dates[i])
  fp_hour <- subset(fp_sample, fp_sample$Date == ysi_dates[i])
  ysi_profile <- ysi_sample[ysi_sample[, "Hour"] == closest(ysi_sample$Hour,fp_hour$Hour[1]),]
  final <- bind_rows(final, ysi_profile)
}

final <- final %>%
  select(Date, Depth_m, Temp_C, DO_mgL, pH)

df.final<-data.frame()

for (i in 1:length(dates)){
  
  profile <- final %>%
    filter(Date == dates[i])
  
  ysi_layer<-profile %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
  
  # Bind each of the data layers together.
  df.final = bind_rows(df.final, ysi_layer)
  
}

final.grab <- ysi[0,]

ysi_dates <- unique(ysi$Date)
for (i in 1:length(ysi_dates)){
  ysi_sample <- subset(ysi, ysi$Date == ysi_dates[i])
  fp_hour <- subset(fp_sample, fp_sample$Date == ysi_dates[i])
  ysi_profile <- ysi_sample[ysi_sample[, "Hour"] == closest(ysi_sample$Hour,fp_hour$Hour[1]),]
  final.grab <- bind_rows(final.grab, ysi_profile)
}

final.grab <- final.grab %>%
  select(Date, Depth_m, Temp_C)

df.final.grab<-data.frame()

for (i in 1:length(grab_dates)){
  
  profile <- final.grab %>%
    filter(Date == grab_dates[i])
  
  ysi_layer<-profile %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - grab_depths[i])))
  
  # Bind each of the data layers together.
  df.final.grab = bind_rows(df.final.grab, ysi_layer)
  
}

colnames(df.final.grab) <-c("Date","Grab_T_Depth_m","Grab_Depth_Temp_C")

check <- left_join(df.final.grab,sample_info,by = "Date")

colnames(ss_ysi)[1] <- "Date"
colnames(td_ysi)[1] <- "Date"
colnames(bf_ysi)[1] <- "Date"
colnames(ln_ysi)[1] <- "Date"
colnames(wn_ysi)[1] <- "Date"


ysi_all <- left_join(df.final, ss_ysi, by = "Date")
ysi_all1 <- left_join(ysi_all, td_ysi, by = "Date")
ysi_all2 <- left_join(ysi_all1, bf_ysi, by = "Date")
ysi_all3 <- left_join(ysi_all2, ln_ysi, by = "Date")
ysi_all4 <- left_join(ysi_all3, wn_ysi, by = "Date")
ysi_final <- left_join(ysi_all4, df.final.grab, by = "Date")


ysi_final <- ysi_final[-c(66:69, 71),]

####SCC THERMISTORS####
#read in SCC data and limit to temperature
scc <- fread("./0_Data_files/SCC.csv")
scc <- tibble(scc) %>%
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

wtr_scc <- scc_new %>%
  mutate(datetime = format(datetime, format = "%Y-%m-%d %H:%M:%S"))
write.table(wtr_scc, "./0_Data_files/wtr_scc.wtr",sep='\t',row.names = FALSE)
wtr_scc <- load.ts("./0_Data_files/wtr_scc.wtr")
bathy <- load.bathy("./0_Data_files/FCR.bth")

#calculate Schmidt stability
ss_scc <- ts.schmidt.stability(wtr = wtr_scc, bathy = bathy, na.rm = TRUE)
ss_scc[ss_scc<0] <- 0

ss_scc1 <- ss_scc %>%
  mutate(Year = year(datetime))

ggplot(data = ss_scc1, aes(x = datetime, y = schmidt.stability))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()

#calculate thermocline depth
td_scc <- ts.thermo.depth(wtr = wtr_scc, na.rm = TRUE, seasonal = TRUE) 

td_scc1 <- td_scc %>%
  mutate(Year = year(datetime))

ggplot(data = td_scc1, aes(x = datetime, y = thermo.depth))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  scale_y_reverse()+
  theme_classic()

#calculate buoyancy frequency
bf_scc <- ts.buoyancy.freq(wtr = wtr_scc, na.rm = TRUE, seasonal = TRUE) 

bf_scc1 <- bf_scc %>%
  mutate(Year = year(datetime))

ggplot(data = bf_scc1, aes(x = datetime, y = n2))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  scale_y_reverse()+
  theme_classic()

#calculate lake number
wnd_scc <- wnd %>%
  filter(datetime %in% unique(wtr_scc$datetime)) 

ln_scc <- ts.lake.number(wtr = wtr_scc, wnd = wnd_scc, wnd.height = 3, bathy = bathy, seasonal = TRUE) 

ln_scc1 <- ln_scc %>%
  mutate(Year = year(datetime))

ggplot(data = ln_scc1, aes(x = datetime, y = lake.number))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()

#calculate wedderburn number
wn_scc <- ts.wedderburn.number(wtr = wtr_scc, wnd = wnd_scc, wnd.height = 3, bathy = bathy, Ao = 119000, seasonal = TRUE) 

wn_scc1 <- wn_scc %>%
  mutate(Year = year(datetime))

ggplot(data = wn_scc1, aes(x = datetime, y = wedderburn.number))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()


#find water temperature at depth of Cmax
cmax <- read_csv("./0_Data_files/FP_DistributionMetrics.csv")

depths = cmax$Peak_depth_m
dates  = cmax$Date

df.final<-data.frame()

scc_long <- wtr_scc %>%
  gather(wtr_0.1:wtr_9, key = "Depth",value = "Temp_C")%>%
  mutate(Depth_m = fun1(strsplit(Depth, split = "_"),2)) %>%
  select(-Depth)

scc_long <- arrange(scc_long, datetime)


for (i in 1:length(dates)){
  
  profile <- scc_long %>%
    filter(datetime == dates[i])
  
  scc_layer<-profile %>% group_by(datetime) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
  
  # Bind each of the data layers together.
  df.final = bind_rows(df.final, scc_layer)
  
}

df.final.grab<-data.frame()


for (i in 1:length(grab_dates)){
  
  profile <- scc_long %>%
    filter(datetime == grab_dates[i])
  
  scc_layer<-profile %>% group_by(datetime) %>% slice(which.min(abs(as.numeric(Depth_m) - grab_depths[i])))
  
  # Bind each of the data layers together.
  df.final.grab = bind_rows(df.final.grab, scc_layer)
  
}

colnames(df.final.grab) <-c("Date","Grab_Depth_Temp_C","Grab_T_Depth_m")

check <- left_join(df.final.grab,sample_info,by = "Date")

colnames(df.final)[1] <- "Date"
colnames(ss_scc)[1] <- "Date"
colnames(td_scc)[1] <- "Date"
colnames(bf_scc)[1] <- "Date"
colnames(ln_scc)[1] <- "Date"
colnames(wn_scc)[1] <- "Date"


scc_all <- left_join(df.final, ss_scc, by = "Date")
scc_all1 <- left_join(scc_all, td_scc, by = "Date")
scc_all2 <- left_join(scc_all1, bf_scc, by = "Date")
scc_all3 <- left_join(scc_all2, ln_scc, by = "Date")
scc_all4 <- left_join(scc_all3, wn_scc, by = "Date")
scc_final <- left_join(scc_all4, df.final.grab, by = "Date")

####FINAL COLLATION OF DATA####
#join data frames from three data sources (CTD, YSI, SCC thermistor string) into one final data fram
Tmetrics <- data.frame(sample_info$Date)
colnames(Tmetrics) <- "Date"

Tmetrics1 <- left_join(Tmetrics, ctd_final, by = "Date")  
colnames(Tmetrics1)[2:12] <- paste("CTD", colnames(Tmetrics1)[2:12], sep = "_")
Tmetrics2 <- left_join(Tmetrics1, ysi_final, by = "Date")
colnames(Tmetrics2)[13:23] <- paste("YSI", colnames(Tmetrics2)[13:23], sep = "_")
Tmetrics3 <- left_join(Tmetrics2, scc_final, by = "Date")
colnames(Tmetrics3)[24:32] <- paste("SCC", colnames(Tmetrics3)[24:32], sep = "_")

Tmetrics4 <- Tmetrics3

write.csv(Tmetrics4, "./0_Data_files/WtrTemp_Stability_DO_pH.csv",row.names = FALSE)

####WRITE PLOTS OF TEMP PROFILE AND THERMOCLINE FOR REVIEWER 2####
ctd_profiles <- wtr_ctd %>%
  gather(wtr_0.1:wtr_9.3, key = "depth", value = "temp") %>%
  mutate(depth = fun1(strsplit(depth, split = "_"),2)) %>%
  arrange(datetime, depth)

ysi_profiles <- wtr_ysi %>%
  gather(wtr_0.1:wtr_9.5, key = "depth", value = "temp") %>%
  mutate(depth = fun1(strsplit(depth, split = "_"),2)) %>%
  arrange(datetime, depth)

scc_profiles <- wtr_scc %>%
  gather(wtr_0.1:wtr_9, key = "depth", value = "temp") %>%
  mutate(depth = fun1(strsplit(depth, split = "_"),2)) %>%
  arrange(datetime, depth)

td.data <- Tmetrics4 %>%
  mutate(thermo.depth = ifelse((is.na(CTD_thermo.depth) & is.na(YSI_thermo.depth)),SCC_thermo.depth,
                               ifelse(is.na(CTD_thermo.depth),YSI_thermo.depth,CTD_thermo.depth)),
         sensor = ifelse((is.na(CTD_thermo.depth) & is.na(YSI_thermo.depth)),"SCC",
                         ifelse(is.na(CTD_thermo.depth),"YSI","CTD"))) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") %>%
  select(Date, thermo.depth, sensor) %>%
  filter(complete.cases(.))

for(i in 1:length(td.data$Date)){
  if(td.data[i,"sensor"] == "CTD"){
    temp_profile <- ctd_profiles %>%
      filter(datetime == td.data$Date[i])
  } else if(td.data[i,"sensor"] == "YSI"){
    temp_profile <- ysi_profiles %>%
      filter(datetime == td.data$Date[i])
  } else {
    temp_profile <- scc_profiles %>%
      filter(datetime == td.data$Date[i])
  }
  
  temp_profile <- tibble(temp_profile) %>%
    mutate(depth = as.numeric(depth))
  
  p1 <- ggplot(data = temp_profile, aes(x = temp, y = depth))+
    geom_path()+
    geom_point()+
    scale_y_reverse(c(0,10))+
    xlab("Water temperature (degrees Celsuis)")+
    ylab("Depth (m)")+
    xlim(c(5,30))+
    ggtitle(temp_profile$datetime[1])+
    geom_hline(yintercept = td.data[i,2], color = "red")+
    theme_classic()
  
  my.filename = paste0("C:/Users/Mary Lofton/Dropbox/Ch_2/FWB_revision_files/Temp_profile_plots/",temp_profile$datetime[1],".png")
  ppi = 300
  png(file = my.filename, width = 3*ppi, height = 3*ppi, res = ppi)
  print(p1)
  dev.off()
}

