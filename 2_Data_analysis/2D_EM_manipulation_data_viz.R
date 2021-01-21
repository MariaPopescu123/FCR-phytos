#3B_EM_manipulation_data_viz
#Author: Mary Lofton
#Date: 17DEC20

pacman::p_load(tidyverse, lubridate)

mydata <- read_csv("./2_Data_analysis/megamatrix.csv") %>%
  mutate(MonthDay = format(Date, format="%m-%d"),
         Year = as.factor(Year)) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20")

p1 <- ggplot(data = mydata, aes(x = MonthDay, y = EM, group = Year, color = Year))+
  geom_line(size = 2)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  theme_classic()
p1

p2 <- ggplot(data = mydata, aes(x = MonthDay, y = thermo.depth, group = Year, color = Year))+
  geom_line(size = 2)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  theme_classic()
p2

my.aov <- aov(thermo.depth ~ Year, data = mydata)
summary(my.aov)
my.tukey <- TukeyHSD(my.aov)
my.tukey

p3 <- ggplot(data = mydata, aes(x = MonthDay, y = Peak_depth_m, group = Year, color = Year))+
  geom_line(size = 2)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  theme_classic()
p3

my.aov <- aov(Peak_depth_m ~ Year, data = mydata)
summary(my.aov)
my.tukey <- TukeyHSD(my.aov)
my.tukey

em <- subset(mydata, Year %in% c(2016,2017))
no_em <- subset(mydata, Year %in% c(2018,2019))
mean(em$Peak_depth_m, na.rm = TRUE)
mean(no_em$Peak_depth_m, na.rm = TRUE)

p4 <- ggplot(data = mydata, aes(x = MonthDay, y = schmidt.stability, group = Year, color = Year))+
  geom_line(size = 2)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  theme_classic()
p4

my.aov <- aov(schmidt.stability ~ Year, data = mydata)
summary(my.aov)
my.tukey <- TukeyHSD(my.aov)
my.tukey

p5 <- ggplot(data = mydata, aes(x = MonthDay, y = perc_light_thermocline, group = Year, color = Year))+
  geom_line(size = 2)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  theme_classic()
p5

my.aov <- aov(perc_light_thermocline ~ Year, data = mydata)
summary(my.aov)
my.tukey <- TukeyHSD(my.aov)
my.tukey

p6 <- ggplot(data = mydata, aes(x = MonthDay, y = DOCmax_mgL, group = Year, color = Year))+
  geom_line(size = 2)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  theme_classic()
p6

my.aov <- aov(DOCmax_depth_m ~ Year, data = mydata)
summary(my.aov)
my.tukey <- TukeyHSD(my.aov)
my.tukey

##THIS NEEDS TO BE EDITED NOW THAT HAVE PERC LIGHT AT THERMOCLINE
for (i in c(6:14,16:20)){
  var <- colnames(mydata[,i])
pa <- ggplot(data = mydata, aes_string(x = var, group = "Year", color = "Year", fill = "Year"))+
  geom_density(alpha = 0.5)+
  theme_classic()
print(pa)

print(var)
my.var <- unlist(mydata[,i])
my.aov <- aov(my.var ~ mydata$Year)
print(summary(my.aov))
my.tukey <- TukeyHSD(my.aov)
print(my.tukey)
}
#DOCmax_mgL not significant, but looks like maybe fairly different
#between EM vs. non-EM years?
#Ditto for Cmax_DOC_mgL but 2019-2017 are similar so likely more
#responsive to HOx operation??
#next step look @ ks tests b/c even if means similar, distributions look
#super different

mydata1 <- mydata %>%
  filter(abs(Peak_depth_m - Phyto_Depth_m) <= 1.1) %>%
  gather(Peak_depth_m, Phyto_Depth_m, key = "sample_type",value = "m")
ggplot(data = mydata1, aes(x = Date, y = m, group = sample_type, color = sample_type))+
  geom_point(size = 2)+
  theme_classic()

#THIS NEEDS TO BE EDITED NOW THAT HAVE PERC LIGHT AT THERMOCLINE
for (i in c(20:46)){
  var <- colnames(mydata1[,i])
  pb <- ggplot(data = mydata1, aes_string(x = var, group = "Year", color = "Year", fill = "Year"))+
    geom_density(alpha = 0.5)+
    theme_classic()
  print(pb)
  
  print(var)
  my.var <- unlist(mydata1[,i])
  my.aov <- aov(my.var ~ mydata1$Year)
  print(summary(my.aov))
  my.tukey <- TukeyHSD(my.aov)
  print(my.tukey)
}


