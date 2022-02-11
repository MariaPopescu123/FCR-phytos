#3B_EM_manipulation_data_viz
#Author: Mary Lofton
#Date: 17DEC20

#load packages and clear environment
pacman::p_load(tidyverse, lubridate, cowplot, grid, data.table, viridis)
rm(list=ls())

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
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") %>%
  mutate(Year = year(Date)) %>%
  filter(Year %in% c(2016:2019))

met3 <- met %>%
  mutate(Year = year(Date)) %>%
  group_by(Year) %>%
  summarize(annual_precip = sum(daily_precip))

rain_storms <- met2 %>%
  filter(daily_precip > quantile(daily_precip, 0.95)) %>%
  mutate(rain_storm = "TRUE")
quantile(met2$daily_precip, 0.95)

wind_storms <- met2 %>%
  filter(daily_wind > quantile(daily_wind, 0.95)) %>%
  mutate(wind_storm = "TRUE")
quantile(met2$daily_wind, 0.95)

all_storms <- full_join(rain_storms, wind_storms)



#get data
my.fp.data <- read_csv("./2_Data_analysis/FP_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_DO_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") 
colnames(my.fp.data)
my.fp.data <- my.fp.data[,c(2,1,3:48)]

my.cs.data <- read_csv("./2_Data_analysis/CS_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_DO_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") %>%
  select(Date,shannon:BV_TOTAL) %>%
  filter(!is.na(BV_TOTAL)) 

mydata <- left_join(my.cs.data,my.fp.data, by = "Date") %>%
  mutate(Year = as.factor(Year), MonthDay = format(Date, format="%m-%d"),
         Peak_depth_m = ifelse(Peak_depth_m > 9, NA, Peak_depth_m)) %>%
  select(Year, Date, MonthDay, thermo.depth, schmidt.stability, Peak_depth_m, n2,SRPmax_depth_m)

my.cols <- viridis(12)
barplot(1:12, col = viridis(12))
#use 5-6 and 11-12

p2016 <- ggplot(data = subset(mydata, mydata$Year == 2016), aes(x = Date, y = thermo.depth))+
  geom_line(color = my.cols[5], size = 1.5)+
  geom_point(color = my.cols[5], size = 3)+
  scale_y_reverse()+
  ylim(6,0)+
  geom_vline(xintercept = as.numeric(as.Date("2016-05-05")))+
  geom_vline(xintercept = as.numeric(as.Date("2016-05-30")), lty = 2)+
  geom_vline(xintercept = as.numeric(as.Date("2016-06-26")), lty = 2)+
  geom_vline(xintercept = as.numeric(as.Date("2016-07-25")), lty = 2)+
  ylab("Thermocline depth (m)")+
  xlab("")+
  ggtitle("C. 2016")+
  theme_classic()+
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 16),plot.title = element_text(size = 18, face = "bold"))

p2017 <- ggplot(data = subset(mydata, mydata$Year == 2017), aes(x = Date, y = thermo.depth))+
  geom_line(color = my.cols[5], size = 1.5)+
  geom_point(color = my.cols[5], size = 3)+
  scale_y_reverse()+
  ylim(6,0)+
  geom_vline(xintercept = as.numeric(as.Date("2017-05-30")), lty = 2)+
  geom_vline(xintercept = as.numeric(as.Date("2017-07-10")), lty = 2)+
  ylab("Thermocline depth (m)")+
  xlab("")+
  ggtitle("F. 2017")+
  theme_classic()+
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 16),plot.title = element_text(size = 18, face = "bold"))

p2018 <- ggplot(data = subset(mydata, mydata$Year == 2018), aes(x = Date, y = thermo.depth))+
  geom_line(color = my.cols[11], size = 1.5)+
  geom_point(color = my.cols[11], size = 3)+
  scale_y_reverse()+
  ylim(6,0)+
  ylab("Thermocline depth (m)")+
  xlab("")+
  ggtitle("I. 2018")+
  theme_classic()+
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 16),plot.title = element_text(size = 18, face = "bold"))

p2019 <- ggplot(data = subset(mydata, mydata$Year == 2019), aes(x = Date, y = thermo.depth))+
  geom_line(color = my.cols[11], size = 1.5)+
  geom_point(color = my.cols[11], size = 3)+
  scale_y_reverse()+
  ylim(6,0)+
  ylab("Thermocline depth (m)")+
  geom_vline(xintercept = as.numeric(as.Date("2019-06-08")))+
  xlab("")+
  ggtitle("L. 2019")+
  theme_classic()+
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 16),plot.title = element_text(size = 18, face = "bold"))


# figs3 <- plot_grid(p2016,p2017,p2018,p2019,nrow = 2, ncol = 2, align = "hv")
# ggsave(figs3, filename = "./3_Visualization/FigS3.tif",height = 6, width = 9,
#        units = "in", dpi = 300, dev = "tiff")

###SI fig for wind/rain data
quantile(met2$daily_precip, 0.95)
quantile(met2$daily_wind, 0.95)


pp2016 <- ggplot(data = subset(met2, met2$Year == 2016 & met2$Date >= "2016-05-02" & met2$Date <= "2016-09-20"), aes(x = Date, y = daily_precip))+
  geom_line(color = my.cols[5], size = 1.5)+
  geom_point(color = my.cols[5], size = 3)+
  geom_vline(xintercept = as.numeric(as.Date("2016-05-05")))+
  geom_hline(yintercept = 22.1, lty = "dashed")+
  geom_vline(xintercept = as.numeric(as.Date("2016-05-30")), lty = 2)+
  geom_vline(xintercept = as.numeric(as.Date("2016-06-26")), lty = 2)+
  geom_vline(xintercept = as.numeric(as.Date("2016-07-25")), lty = 2)+
  ylab("Precipitation (mm)")+
  xlab("")+
  ggtitle("A. 2016")+
  theme_classic()+
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 16),plot.title = element_text(size = 18, face = "bold"))
#pp2016
pp2017 <- ggplot(data = subset(met2, met2$Year == 2017 & met2$Date >= "2017-05-15" & met2$Date <= "2017-09-04"), aes(x = Date, y = daily_precip))+
  geom_line(color = my.cols[5], size = 1.5)+
  geom_point(color = my.cols[5], size = 3)+
  geom_hline(yintercept = 22.1, lty = "dashed")+
  geom_vline(xintercept = as.numeric(as.Date("2017-05-30")), lty = 2)+
  geom_vline(xintercept = as.numeric(as.Date("2017-07-10")), lty = 2)+
  ylab("Precipitation (mm)")+
  xlab("")+
  ggtitle("D. 2017")+
  theme_classic()+
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 16),plot.title = element_text(size = 18, face = "bold"))
#pp2017
pp2018 <- ggplot(data = subset(met2, met2$Year == 2018& met2$Date >= "2018-05-07" & met2$Date <= "2018-09-10"), aes(x = Date, y = daily_precip))+
  geom_line(color = my.cols[11], size = 1.5)+
  geom_point(color = my.cols[11], size = 3)+
  geom_hline(yintercept = 22.1, lty = "dashed")+
  ylab("Precipitation (mm)")+
  xlab("")+
  ggtitle("G. 2018")+
  theme_classic()+
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 16),plot.title = element_text(size = 18, face = "bold"))
#pp2018
pp2019 <- ggplot(data = subset(met2, met2$Year == 2019& met2$Date >= "2019-05-06" & met2$Date <= "2019-09-11"), aes(x = Date, y = daily_precip))+
  geom_line(color = my.cols[11], size = 1.5)+
  geom_point(color = my.cols[11], size = 3)+
  geom_vline(xintercept = as.numeric(as.Date("2019-06-08")))+
  geom_hline(yintercept = 22.1, lty = "dashed")+
  ylab("Precipitation (mm)")+
  xlab("")+
  ggtitle("J. 2019")+
  theme_classic()+
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 16),plot.title = element_text(size = 18, face = "bold"))
#pp2019

pw2016 <- ggplot(data = subset(met2, met2$Year == 2016& met2$Date >= "2016-05-02" & met2$Date <= "2016-09-20"), aes(x = Date, y = daily_wind))+
  geom_line(color = my.cols[5], size = 1.5)+
  geom_point(color = my.cols[5], size = 3)+
  geom_vline(xintercept = as.numeric(as.Date("2016-05-05")))+
  geom_hline(yintercept = 2.7, lty = "dashed")+
  geom_vline(xintercept = as.numeric(as.Date("2016-05-30")), lty = 2)+
  geom_vline(xintercept = as.numeric(as.Date("2016-06-26")), lty = 2)+
  geom_vline(xintercept = as.numeric(as.Date("2016-07-25")), lty = 2)+
  ylab(expression(paste("Daily mean windspeed  ","(",m~s^-1,")")))+
  xlab("")+
  ggtitle("B. 2016")+
  theme_classic()+
  theme(axis.title = element_text(size = 14),axis.text = element_text(size = 16),plot.title = element_text(size = 18, face = "bold"))
#pw2016
pw2017 <- ggplot(data = subset(met2, met2$Year == 2017& met2$Date >= "2017-05-15" & met2$Date <= "2017-09-04"), aes(x = Date, y = daily_wind))+
  geom_line(color = my.cols[5], size = 1.5)+
  geom_point(color = my.cols[5], size = 3)+
  geom_hline(yintercept = 2.7, lty = "dashed")+
  geom_vline(xintercept = as.numeric(as.Date("2017-05-30")), lty = 2)+
  geom_vline(xintercept = as.numeric(as.Date("2017-07-10")), lty = 2)+
  ylab(expression(paste("Daily mean windspeed  ","(",m~s^-1,")")))+
  xlab("")+
  ggtitle("E. 2017")+
  theme_classic()+
  theme(axis.title = element_text(size = 14),axis.text = element_text(size = 16),plot.title = element_text(size = 18, face = "bold"))
#pw2017
pw2018 <- ggplot(data = subset(met2, met2$Year == 2018& met2$Date >= "2018-05-07" & met2$Date <= "2018-09-10"), aes(x = Date, y = daily_wind))+
  geom_line(color = my.cols[11], size = 1.5)+
  geom_point(color = my.cols[11], size = 3)+
  geom_hline(yintercept = 2.7, lty = "dashed")+
  ylab(expression(paste("Daily mean windspeed  ","(",m~s^-1,")")))+
  xlab("")+
  ggtitle("H. 2018")+
  theme_classic()+
  theme(axis.title = element_text(size = 14),axis.text = element_text(size = 16),plot.title = element_text(size = 18, face = "bold"))
#pw2018
pw2019 <- ggplot(data = subset(met2, met2$Year == 2019 & met2$Date >= "2019-05-06" & met2$Date <= "2019-09-11"), aes(x = Date, y = daily_wind))+
  geom_line(color = my.cols[11], size = 1.5)+
  geom_point(color = my.cols[11], size = 3)+
  geom_vline(xintercept = as.numeric(as.Date("2019-06-08")))+
  geom_hline(yintercept = 2.7, lty = "dashed")+
  ylab(expression(paste("Daily mean windspeed  ","(",m~s^-1,")")))+
  xlab("")+
  ggtitle("K. 2019")+
  theme_classic()+
  theme(axis.title = element_text(size = 14),axis.text = element_text(size = 16),plot.title = element_text(size = 18, face = "bold"))
#pw2019

# figs3 <- plot_grid(pp2016,pw2016,pp2017,pw2017,pp2018,pw2018,pp2019,pw2019,nrow = 4, ncol = 2, align = "hv")
# ggsave(figs3, filename = "./3_Visualization/met_data.tif",height = 12, width = 9,
#        units = "in", dpi = 300, dev = "tiff")

fig3_revision <- plot_grid(pp2016,pw2016,p2016, pp2017,pw2017, p2017, pp2018,pw2018,p2018, pp2019,pw2019, p2019, nrow = 4, ncol = 3, align = "hv")
ggsave(fig3_revision, filename = "./3_Visualization/Fig3_revision.tif",height = 12, width = 13.5,
       units = "in", dpi = 300, dev = "tiff")

