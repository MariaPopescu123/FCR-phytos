#3B_EM_manipulation_data_viz
#Author: Mary Lofton
#Date: 17DEC20

pacman::p_load(tidyverse, lubridate, cowplot, grid, data.table)
rm(list=ls())



#get data
mydata <- read_csv( "./00_Data_files/WRT.csv") %>%
  mutate(Year = year(Date))


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

my.cols <- gg_color_hue(4)

p2016 <- ggplot(data = subset(mydata, mydata$Year == 2016), aes(x = Date, y = WRT_day))+
  geom_line(color = my.cols[1], size = 1.5)+
  geom_point(color = my.cols[1], size = 3)+
  ylab("Water residence time (days)")+
  xlab("")+
  ggtitle("A. 2016")+
  theme_classic()+
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 16),plot.title = element_text(size = 18))
p2016
p2017 <- ggplot(data = subset(mydata, mydata$Year == 2017), aes(x = Date, y = WRT_day))+
  geom_line(color = my.cols[2], size = 1.5)+
  geom_point(color = my.cols[2], size = 3)+
  ylab("Water residence time (days)")+
  xlab("")+
  ggtitle("B. 2017")+
  theme_classic()+
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 16),plot.title = element_text(size = 18))
p2017
p2018 <- ggplot(data = subset(mydata, mydata$Year == 2018), aes(x = Date, y = WRT_day))+
  geom_line(color = my.cols[3], size = 1.5)+
  geom_point(color = my.cols[3], size = 3)+
  ylab("Water residence time (days)")+
  xlab("")+
  ggtitle("C. 2018")+
  theme_classic()+
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 16),plot.title = element_text(size = 18))
p2018
p2019 <- ggplot(data = subset(mydata, mydata$Year == 2019), aes(x = Date, y = WRT_day))+
  geom_line(color = my.cols[4], size = 1.5)+
  geom_point(color = my.cols[4], size = 3)+
  ylab("Water residence time (days)")+
  xlab("")+
  ggtitle("D. 2019")+
  theme_classic()+
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 16),plot.title = element_text(size = 18))
p2019

figs3 <- plot_grid(p2016,p2017,p2018,p2019,nrow = 2, ncol = 2, align = "hv")
ggsave(figs3, filename = "./3_Visualization/FigS6.tif",height = 6, width = 9,
       units = "in", dpi = 300, dev = "tiff")

