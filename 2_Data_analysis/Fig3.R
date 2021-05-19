#3B_EM_manipulation_data_viz
#Author: Mary Lofton
#Date: 17DEC20

pacman::p_load(tidyverse, lubridate, cowplot)
rm(list=ls())


#get data
my.fp.data <- read_csv("./2_Data_analysis/FP_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_DO_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") 
colnames(my.fp.data)
my.fp.data <- my.fp.data[,c(2,1,3:34)]

my.cs.data <- read_csv("./2_Data_analysis/CS_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_DO_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") %>%
  select(Date,shannon:BV_TOTAL) %>%
  filter(!is.na(BV_TOTAL)) 

mydata <- left_join(my.cs.data,my.fp.data, by = "Date") %>%
  mutate(Year = as.factor(Year), MonthDay = format(Date, format="%m-%d"),
         Peak_depth_m = ifelse(Peak_depth_m > 9, NA, Peak_depth_m)) %>%
  select(Year, Date, MonthDay, thermo.depth, schmidt.stability, Peak_depth_m)

axis_dates <- mydata[1:3,] 
axis_dates$Year[1:3]<-NA
axis_dates$Date[1:3]<- c("2016-05-01","2016-06-01","2016-09-01")
axis_dates$MonthDay<- c("05-01","06-01","09-01")
axis_dates[,4:6] <- NA
distdates <- rbind(mydata, axis_dates) 

p1 <- ggplot(data = distdates, aes(x = MonthDay, y = thermo.depth, group = Year, color = Year))+
  geom_line(size = 1.5,)+
  geom_point(size = 3)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  scale_color_discrete(na.translate = FALSE)+
  ylab("Thermocline depth (m)")+
  xlab("Date (mm-dd)")+
  theme_classic()+
  ggtitle("A")
p1

p2 <- ggplot(data = distdates, aes(x = MonthDay, y = schmidt.stability, group = Year, color = Year))+
  geom_line(size = 1.5,)+
  geom_point(size = 3)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  scale_color_discrete(na.translate = FALSE)+
  ylab("Schmidt stability (J/m2)")+
  xlab("Date (mm-dd)")+
  theme_classic()+
  ggtitle("C")
p2

p3 <- ggplot(data = distdates, aes(x = MonthDay, y = Peak_depth_m, group = Year, color = Year))+
  geom_line(size = 1.5,)+
  geom_point(size = 3)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  scale_color_discrete(na.translate = FALSE)+
  ylab("Depth of peak biomass (m)")+
  xlab("Date (mm-dd)")+
  theme_classic()+
  ggtitle("E")
p3

p4 <- ggplot(data = distdates, aes(x = thermo.depth, group = Year, color = Year, fill = Year))+
  geom_density(alpha = 0.5)+
  xlab("Thermocline depth (m)")+
  ylab("Density")+
  theme_classic()+
  ggtitle("B")
p4

p5 <- ggplot(data = distdates, aes(x = schmidt.stability, group = Year, color = Year, fill = Year))+
  geom_density(alpha = 0.5)+
  xlab("Schmidt stability (J/m2)")+
  ylab("Density")+
  theme_classic()+
  ggtitle("D")
p5

p6 <- ggplot(data = distdates, aes(x = Peak_depth_m, group = Year, color = Year, fill = Year))+
  geom_density(alpha = 0.5)+
  xlab("Depth of peak biomass (m)")+
  ylab("Density")+
  theme_classic()+
  ggtitle("F")
p6

fig3 <- plot_grid(p1,p4,p2,p5,p3,p6,nrow = 3, ncol = 2, align = "hv")
ggsave(fig3, filename = "./3_Visualization/Fig3.tif",height = 9, width = 9,
       units = "in", dpi = 300, dev = "tiff")
