#FigS3
#Author: Mary Lofton
#Date: 17DEC20
#Updated Feb. 2022

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
  select(Year, Date, MonthDay, thermo.depth, schmidt.stability,n2, SRPmax_depth_m, pz_DOC_mean,DOCmax_mgL,Cmax_DOC_mgL,Grab_DOC_mgL,Peak_depth_m,BV_Desmids,rel_abund_Desmids)

axis_dates <- mydata[1:3,] 
axis_dates$Year[1:3]<-NA
axis_dates$Date[1:3]<- c("2016-05-01","2016-06-01","2016-09-01")
axis_dates$MonthDay<- c("05-01","06-01","09-01")
axis_dates[,4:7] <- NA
distdates <- rbind(mydata, axis_dates) %>%
  mutate(EM2 = as.factor(ifelse(Year %in% c(2016:2017),"mix","no mix")))

p1 <- ggplot(data = distdates, aes(x = MonthDay, y = thermo.depth, group = Year, color = Year))+
  geom_line(size = 1.5,)+
  geom_point(size = 3)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  scale_color_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  ylab("Thermocline depth (m)")+
  xlab("Date (mm-dd)")+
  theme_classic()+  
  theme(legend.position = "none",plot.title = element_text(face = "bold"))+
  scale_y_reverse()+
  ggtitle("A")
p1

p2 <- ggplot(data = distdates, aes(x = MonthDay, y = schmidt.stability, group = Year, color = Year))+
  geom_line(size = 1.5,)+
  geom_point(size = 3)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  scale_color_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  ylab(expression(paste("Schmidt stability ","(",J,~m^-2,")")))+
  xlab("Date (mm-dd)")+
  theme_classic()+
  theme(legend.position = "none",plot.title = element_text(face = "bold"))+
  ggtitle("B")
p2

p3 <- ggplot(data = distdates, aes(x = MonthDay, y = n2, group = Year, color = Year))+
  geom_line(size = 1.5,)+
  geom_point(size = 3)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  scale_color_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  ylab(expression(paste("Buoyancy frequency ","(",~s^-1,")")))+
  scale_y_continuous(breaks = c(0.002,0.004,0.006,0.008,0.010,0.012,0.014))+
  xlab("Date (mm-dd)")+
  theme_classic()+
  theme(legend.position = "none",plot.title = element_text(face = "bold"))+
  ggtitle("C")
p3

p4 <- ggplot(data = distdates, aes(x = MonthDay, y = SRPmax_depth_m, group = Year, color = Year))+
  #geom_line(size = 1.5)+
  geom_point(size = 3)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  scale_color_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  ylab(expression(paste("Depth of maximum SRP (m)")))+
  xlab("Date (mm-dd)")+
  theme_classic()+
  theme(legend.position = "none",plot.title = element_text(face = "bold"))+
  scale_y_reverse()+
  ggtitle("D")
p4

p5 <- ggplot(data = distdates, aes(x = MonthDay, y = pz_DOC_mean, group = Year, color = Year))+
  geom_line(size = 1.5)+
  geom_point(size = 3)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  scale_color_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  ylab(expression(paste("Mean photic zone DOC (",mg,~L^-1,")")))+
  xlab("Date (mm-dd)")+
  theme_classic()+
  theme(legend.position = "none",plot.title = element_text(face = "bold"))+
  ggtitle("E")
p5

p6 <- ggplot(data = distdates, aes(x = MonthDay, y = DOCmax_mgL, group = Year, color = Year))+
  geom_line(size = 1.5)+
  geom_point(size = 3)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  scale_color_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  ylab(expression(paste("Maximum photic zone DOC (",mg,~L^-1,")")))+
  xlab("Date (mm-dd)")+
  theme_classic()+
  theme(legend.position = "none",plot.title = element_text(face = "bold"))+
  ggtitle("F")
p6

p7 <- ggplot(data = distdates, aes(x = MonthDay, y = Cmax_DOC_mgL, group = Year, color = Year))+
  geom_line(size = 1.5)+
  geom_point(size = 3)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  scale_color_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  ylab(expression(paste("DOC at peak depth (",mg,~L^-1,")")))+
  xlab("Date (mm-dd)")+
  theme_classic()+
  theme(legend.position = "none",plot.title = element_text(face = "bold"))+
  ggtitle("G")
p7

p8 <- ggplot(data = distdates, aes(x = MonthDay, y = Grab_DOC_mgL, group = Year, color = Year))+
  geom_line(size = 1.5)+
  geom_point(size = 3)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  scale_color_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  ylab(expression(paste("DOC at depth sample (",mg,~L^-1,")")))+
  xlab("Date (mm-dd)")+
  theme_classic()+
  theme(legend.position = "none",plot.title = element_text(face = "bold"))+
  ggtitle("H")
p8

p9 <- ggplot(data = distdates, aes(x = MonthDay, y = Peak_depth_m, group = Year, color = Year))+
  geom_line(size = 1.5,)+
  geom_point(size = 3)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  scale_color_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  ylab("Depth of peak biomass (m)")+
  xlab("Date (mm-dd)")+
  theme_classic()+
  theme(legend.position = "none",plot.title = element_text(face = "bold"))+
  scale_y_reverse()+
  ggtitle("I")
p9

p10 <- ggplot(data = distdates, aes(x = MonthDay, y = log(BV_Desmids), group = Year, color = Year))+
  geom_line(size = 1.5)+
  geom_point(size = 3)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  scale_color_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  ylab(expression(paste("log(Desmid biovolume ","(",mu,m^3,~mL^-1,"))")))+
  xlab("Date (mm-dd)")+
  theme_classic()+
  theme(legend.position = "none",plot.title = element_text(face = "bold"))+
  ggtitle("J")
p10

p11 <- ggplot(data = distdates, aes(x = MonthDay, y = rel_abund_Desmids, group = Year, color = Year))+
  geom_line(size = 1.5)+
  geom_point(size = 3)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  scale_color_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  ylab("relative abundance of Desmids")+
  xlab("Date (mm-dd)")+
  theme_classic()+
  theme(legend.position = "none",plot.title = element_text(face = "bold"))+
  ggtitle("K")
p11

pleg_line <- ggplot(data = distdates, aes(x = MonthDay, y = Peak_depth_m, group = Year, color = Year))+
  geom_line(size = 1.5,)+
  geom_point(size = 3)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  scale_color_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  ylab("Jaccard dissimilarity index")+
  xlab("Date (mm-dd)")+
  theme_classic()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_blank(), #change legend title font size
        legend.text = element_text(size=16))#change legend text font size   ggtitle("G")
pleg_line

line_legend <- cowplot::get_legend(pleg_line)

grid.newpage()
grid.draw(line_legend)


figs3 <- plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,line_legend,nrow = 4, ncol = 3, align = "hv")
fig3
ggsave(figs3, filename = "./3_Visualization/FigS3_v2.tif",height = 11, width = 8.5,
       units = "in", dpi = 300, dev = "tiff")

