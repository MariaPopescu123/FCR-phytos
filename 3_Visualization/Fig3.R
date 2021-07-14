#Fig3
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
  select(Year, Date, MonthDay, thermo.depth, schmidt.stability, Peak_depth_m, n2,SRPmax_depth_m)

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
  ggtitle("C")
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
  ggtitle("E")
p3

p4 <- ggplot(data = distdates, aes(x = MonthDay, y = SRPmax_depth_m, group = Year, color = Year))+
  #geom_line(size = 1.5)+
  geom_point(size = 3)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  scale_color_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  ylab(expression(paste("Depth of maximum SRP (",mu,g,~L^-1,")")))+
  xlab("Date (mm-dd)")+
  theme_classic()+
  theme(legend.position = "none",plot.title = element_text(face = "bold"))+
  scale_y_reverse()+
  ggtitle("G")
p4

p5 <- ggplot(data = distdates, aes(x = MonthDay, y = Peak_depth_m, group = Year, color = Year))+
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
p5

boxplot_dates <- distdates %>%
  filter(!is.na(Year)) #%>%
  # mutate(td_groups = ifelse(Year %in% c(2016:2017),"a","b")) %>%
  # mutate(td_y = 2) %>%
  # mutate(ss_groups = ifelse(Year %in% c(2016:2017),"a","b")) %>%
  # mutate(ss_y = 50) %>%
  # mutate(pd_groups = ifelse(Year == 2016,"a",
  #                           ifelse(Year == 2017,"ab",
  #                                  ifelse(Year == 2018,"c","bc")))) %>%
  # mutate(pd_y = 0) %>%
  # mutate(jac_groups = ifelse(Year == 2016,"ab",
  #                           ifelse(Year == 2017,"a","b"))) %>%
  # mutate(jac_y = 0.8)%>%
  # mutate(dSRP_groups = ifelse(Year == 2016,"a",
  #                             ifelse(Year == 2017,"ab",
  #                                    ifelse(Year == 2018,"bc","c")))) %>%
  # mutate(dSRP_y = 0) %>%
  # mutate(grabdin_groups = ifelse(Year == 2016,"a",
  #                                ifelse(Year %in% c(2017:2018),"ab","b"))) %>%
  # mutate(grabdin_y = 200) %>%
  # mutate(cvdin_groups = "a") %>%
  # mutate(cvdin_y = 1.6)

p6 <- ggplot(data = boxplot_dates, aes(x = Year, y = thermo.depth, group = Year, fill = Year))+
  geom_boxplot()+
  ylab("Thermocline depth (m)")+
  scale_fill_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  xlab("Year")+
  theme_classic()+
  theme(legend.position = "none",plot.title = element_text(face = "bold"))+
  #geom_text(data = boxplot_dates, aes(y = td_y-0.1, x = Year,label = td_groups), size = 4)+
  scale_y_reverse()+
  ggtitle("B")
p6

p7 <- ggplot(data = boxplot_dates, aes(x = Year, y = schmidt.stability, group = Year, fill = Year))+
  geom_boxplot()+
  ylab(expression(paste("Schmidt stability ","(",J,~m^-2,")")))+
  scale_fill_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  xlab("Year")+
  theme_classic()+
  #geom_text(data = boxplot_dates, aes(y = ss_y+2, x = Year,label = ss_groups), size = 4)+
  theme(legend.position = "none",plot.title = element_text(face = "bold"))+
  ggtitle("D")
p7

p8 <- ggplot(data = boxplot_dates, aes(x = Year, y = n2, group = Year, fill = Year))+
  geom_boxplot()+
  ylab(expression(paste("Buoyancy frequency ","(",~s^-1,")")))+
  scale_fill_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  xlab("Year")+
  theme_classic()+
  scale_y_continuous(breaks = c(0.002,0.004,0.006,0.008,0.010,0.012,0.014))+
  #geom_text(data = boxplot_dates, aes(y = ss_y+2, x = Year,label = ss_groups), size = 4)+
  theme(legend.position = "none",plot.title = element_text(face = "bold"))+
  ggtitle("F")
p8

p9 <- ggplot(data = boxplot_dates, aes(x = Year, y = SRPmax_depth_m, group = Year, fill = Year))+
  geom_boxplot()+
  ylab(expression(paste("Depth of maximum SRP (",mu,g,~L^-1,")")))+
  scale_fill_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  xlab("Year")+
  theme_classic()+
  scale_y_reverse()+
  #geom_text(data = boxplot_dates, aes(y = dSRP_y-0.3, x = Year,label = dSRP_groups), size = 4)+
  theme(legend.position = "none",plot.title = element_text(face = "bold"))+
  ggtitle("H")
p9

p10 <- ggplot(data = boxplot_dates, aes(x = Year, y = Peak_depth_m, group = Year, fill = Year))+
  geom_boxplot()+
  ylab("Depth of peak biomass (m)")+
  scale_fill_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  xlab("Year")+
  theme_classic()+
  #geom_text(data = boxplot_dates, aes(y = pd_y-0.1, x = Year,label = pd_groups), size = 4)+
  theme(legend.position = "none",plot.title = element_text(face = "bold"))+
  scale_y_reverse()+
  ggtitle("J")
p10

pleg_box <- ggplot(data = boxplot_dates, aes(x = Year, y = Peak_depth_m, group = Year, fill = Year))+
  geom_boxplot()+
  ylab("Jaccard dissimilarity index")+
  scale_fill_manual(na.translate = FALSE, values = c("coral","coral3","cadetblue3","cadetblue4"))+
  xlab("Year")+
  theme_classic()+
  #geom_text(data = boxplot_dates, aes(y = jac_y, x = Year,label = jac_groups), size = 4)+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_blank(), #change legend title font size
        legend.text = element_text(size=16),
        legend.position="bottom")
pleg_box

boxplot_legend <- cowplot::get_legend(pleg_box)

grid.newpage()
grid.draw(boxplot_legend)

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
        legend.text = element_text(size=16),#change legend text font size   ggtitle("G")
        legend.position="bottom")
pleg_line

line_legend <- cowplot::get_legend(pleg_line)

grid.newpage()
grid.draw(line_legend)


fig3 <- plot_grid(p1,p6,p2,p7,p3,p8,p4,p9,p5,p10,line_legend,boxplot_legend,nrow = 6, ncol = 2, align = "hv",
                  rel_heights = c(4,4,4,4,4,1))
fig3
ggsave(fig3, filename = "./3_Visualization/Fig3_v2.tif",height = 12, width = 8.5,
       units = "in", dpi = 300, dev = "tiff")

