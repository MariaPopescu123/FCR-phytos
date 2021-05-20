#3B_EM_manipulation_data_viz
#Author: Mary Lofton
#Date: 17DEC20

pacman::p_load(tidyverse, lubridate, cowplot, grid)
rm(list=ls())


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
  select(Year, Date, MonthDay, Max_biomass_ugL, Peak_width_m)

axis_dates <- mydata[1:3,] 
axis_dates$Year[1:3]<-NA
axis_dates$Date[1:3]<- c("2016-05-01","2016-06-01","2016-09-01")
axis_dates$MonthDay<- c("05-01","06-01","09-01")
axis_dates[,4:5] <- NA
distdates <- rbind(mydata, axis_dates) %>%
  mutate(EM2 = as.factor(ifelse(Year %in% c(2016:2017),"mix","no mix")))

mb <- mean(distdates$Max_biomass_ugL, na.rm = TRUE)
sd_bio <- sd(distdates$Max_biomass_ugL, na.rm = TRUE)

mb + 2*sd_bio
#65 ugL

p1 <- ggplot(data = distdates, aes(x = MonthDay, y = Max_biomass_ugL, group = Year, color = Year))+
  geom_line(size = 1.5,)+
  geom_point(size = 3)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  scale_color_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  ylab(expression(paste("Biomass (",mu,g,~L^-1,")")))+
  xlab("Date (mm-dd)")+
  theme_classic()+  
  theme(plot.title = element_text(face = "bold"))+
  geom_hline(yintercept = 65, size = 1.5, lty = "dashed")+
  ggtitle("A")
p1

p2 <- ggplot(data = distdates, aes(x = MonthDay, y = Peak_width_m, group = Year, color = Year))+
  geom_line(size = 1.5,)+
  geom_point(size = 3)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  scale_color_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  ylab("Peak width (m)")+
  xlab("Date (mm-dd)")+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  ggtitle("B")
p2


figs4 <- plot_grid(p1,p2,nrow = 2, ncol = 1, align = "hv")
ggsave(figs4, filename = "./3_Visualization/FigS4.tif",height = 6, width = 6,
       units = "in", dpi = 300, dev = "tiff")

