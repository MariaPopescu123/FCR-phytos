#3B_EM_manipulation_data_viz
#Author: Mary Lofton
#Date: 17DEC20

pacman::p_load(tidyverse, lubridate)
rm(list=ls())


#get data
my.fp.data <- read_csv("./2_Data_analysis/FP_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") 
my.fp.data <- my.fp.data[,c(2,1,3:23)]

my.cs.data <- read_csv("./2_Data_analysis/CS_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") %>%
  select(Date,shannon:BV_TOTAL)

mydata <- left_join(my.fp.data,my.cs.data, by = "Date") %>%
  mutate(Year = as.factor(Year),
         Peak_depth_m = ifelse(Peak_depth_m > 9, NA, Peak_depth_m)) 

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

t.test(thermo.depth ~ EM, data = mydata)


p3 <- ggplot(data = mydata, aes(x = MonthDay, y = Peak_depth_m, group = Year, color = Year))+
  geom_line(size = 2)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  theme_classic()
p3

my.aov <- aov(Peak_depth_m ~ Year, data = mydata)
summary(my.aov)
my.tukey <- TukeyHSD(my.aov)
my.tukey

t.test(Peak_depth_m ~ EM, data = mydata)


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

t.test(schmidt.stability ~ EM, data = mydata)


p5 <- ggplot(data = mydata, aes(x = MonthDay, y = rel_abund_Cyanobacteria, group = Year, color = Year))+
  geom_line(size = 2)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  theme_classic()
p5

p6 <- ggplot(data = mydata, aes(x = MonthDay, y = rel_abund_Cryptophytes, group = Year, color = Year))+
  geom_line(size = 2)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  theme_classic()
p6

p7 <- ggplot(data = mydata, aes(x = MonthDay, y = Max_biomass_ugL, group = Year, color = Year))+
  geom_line(size = 2)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  theme_classic()
p7

p8 <- ggplot(data = mydata, aes(x = MonthDay, y = richness, group = Year, color = Year))+
  geom_line(size = 2)+
  scale_x_discrete(breaks = c("05-01","06-01","07-01","08-01","09-01","10-01"))+
  theme_classic()
p8

##THIS NEEDS TO BE EDITED NOW THAT HAVE PERC LIGHT AT THERMOCLINE
for (i in c(5:22,24:32)){
  var <- colnames(mydata[,i])
pa <- ggplot(data = mydata, aes_string(x = var, group = "Year", color = "Year", fill = "Year"))+
  geom_density(alpha = 0.5)+
  theme_classic()
print(pa)

print(var)
my.var <- unlist(mydata[,i])
# my.aov <- aov(my.var ~ mydata$Year)
# print(summary(my.aov))
# my.tukey <- TukeyHSD(my.aov)
# print(my.tukey)
# 
# print(t.test(my.var ~ EM, data = mydata))
print(t.test(my.var ~ HOx, data = mydata))

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

##
ggplot(data = mydata, aes(x = Date, y = DOCmax_mgL))+
  geom_line()+
  theme_classic()
ggplot(data = mydata, aes(x = Date, y = BV_Cyanobacteria))+
  geom_line()+
  theme_classic()
ggplot(data = mydata, aes(x = Date, y = BV_TOTAL))+
  geom_line()+
  geom_point()+
  theme_classic()
par(mfrow = c(2,2))
ggplot(data = subset(mydata, mydata$Year == 2017), aes(x = Date, y = Cmax_DOC_mgL))+
  geom_line()+
  geom_point()+
  theme_classic()
ggplot(data = subset(mydata, mydata$Year == 2019), aes(x = Date, y = Cmax_DOC_mgL))+
  geom_line()+
  geom_point()+
  theme_classic()
ggplot(data = subset(mydata, mydata$Year == 2017), aes(x = Date, y = DOCmax_mgL))+
  geom_line()+
  geom_point()+
  theme_classic()
ggplot(data = subset(mydata, mydata$Year == 2019), aes(x = Date, y = DOCmax_mgL))+
  geom_line()+
  geom_point()+
  theme_classic()
ggplot(data = subset(mydata, mydata$Year == 2017), aes(x = Date, y = BV_Cyanobacteria))+
  geom_line()+
  geom_point()+
  theme_classic()
ggplot(data = subset(mydata, mydata$Year == 2019), aes(x = Date, y = BV_Cyanobacteria))+
  geom_line()+
  geom_point()+
  theme_classic()

for (i in 3:18){
  my.var <- mydata[,i]
  #print(colnames(my.var))
  #print(cor(my.var, mydata$Peak_depth_m, method = "spearman", use = "complete.obs"))
  #print(cor(my.var, mydata$BV_Cyanobacteria, method = "spearman", use = "complete.obs"))
  print(cor(my.var, mydata$BV_Cryptophytes, method = "spearman", use = "complete.obs"))
  #print(cor(my.var, mydata$DOCmax_mgL, method = "spearman", use = "complete.obs"))
  #print(cor(my.var, mydata$DINmax_ugL, method = "spearman", use = "complete.obs"))
}

for (i in c(19:22, 24:32)){
  my.var <- mydata[,i]
  #print(colnames(my.var))
  #print(cor(my.var, mydata$Peak_depth_m, method = "spearman", use = "complete.obs"))
  #print(cor(my.var, mydata$BV_Cyanobacteria, method = "spearman", use = "complete.obs"))
  print(cor(my.var, mydata$EM, method = "spearman", use = "complete.obs"))
  #print(cor(my.var, mydata$DOCmax_mgL, method = "spearman", use = "complete.obs"))
  #print(cor(my.var, mydata$DINmax_ugL, method = "spearman", use = "complete.obs"))
}

ggplot(data = subset(mydata, month(mydata$Date) != 1), aes(x = Date, y = (thermo.depth - Peak_depth_m)))+
  geom_point()+
  theme_classic()+
  geom_smooth(method = "lm")+
  geom_hline(yintercept = 0)+
  facet_wrap(facets = vars(Year), nrow = 2, scales = "free_x")

allcor <- cor(mydata[,c(3:22, 24:32)], method = "spearman", use = "complete.obs")
write.csv(allcor, "./2_Data_analysis/spearman_results_all.csv",row.names = TRUE)
allcor[allcor<0.5]<- NA
allcor[allcor == 1] <- NA
hist(allcor)
allcor[lower.tri(allcor)] = ""
