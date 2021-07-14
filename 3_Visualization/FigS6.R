#FigS6
#Author: Mary Lofton
#Date: 11DEC2017
#Updated: 01FEB2021

####SET-UP####
#load packages
pacman::p_load(tidyverse, lubridate, vegan, Hmisc, picante, indicspecies, vegan3d)
rm(list=ls())

#set functions
#get colors for plotting
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#read in environmental data
my.cs.data <- read_csv("./2_Data_analysis/CS_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_DO_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") %>%
  select(-MonthDay) %>%
  filter(!is.na(BV_TOTAL)) 

colnames(my.cs.data)
env <- my.cs.data[,c(2:6,7:9,16:21, 23:26, 33)]
env$Month <- month(my.cs.data$Date)
env <- data.frame(scale(env, center = TRUE, scale = TRUE)) 

#read in community data
phytos <- read_csv("./00_Data_files/EDI_phytos/phytoplankton.csv") %>%
  select(Date, Genus, BV_um3mL) %>%
  spread(key = Genus, value = BV_um3mL)

# any non-detects filled with 0
phytos[is.na(phytos)]<- 0

# create data frame with different EM periods for plotting later on
em <- my.cs.data[,c(1,2,4,5,6)] %>%
  mutate(Month = month(Date))

phytos1 <- left_join(em, phytos, by = "Date") %>%
  mutate(MonthDay = format(Date, format="%m-%d"))

# prepare phyto data for NMDS input
phytos2 <- phytos1 %>%
  select(-Date,-EM1,-EM2,-EM3,-Month,-Year,-MonthDay) 

phytos3 <- as.matrix(phytos2)

#### NMDS FOR 2016 ONLY ####

#select 2016 data
p2016.data <- phytos1 %>%
  filter(year(Date) == 2016)

#prep data for NMDS input
p2016 <- p2016.data %>%
  select(-Year,-Date,-EM1,-EM2,-EM3,-Month,-MonthDay)
p2016 <- as.matrix(p2016)



#run NMDS for all months
set.seed(1)
Q <- metaMDS(p2016, distance='bray', k=2, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
# plot(Q, display=c('sites'), choices=c(1,2), type='p')



#coded by month
my.cols <- gg_color_hue(5)
png(filename = "./3_Visualization/Fig6.png",width = 8,height = 8,units = "in",res = 300)
par(mfrow = c(2,2),mar = c(3,3,0,0),mgp = c(1.5,0.5,0))
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')

#May
points(Q$points[p2016.data$Month==5 & p2016.data$EM3 == 0,1], Q$points[p2016.data$Month==5 & p2016.data$EM3 == 0,2], pch=1,col=my.cols[1],bg=my.cols[1], cex = 2)
points(Q$points[p2016.data$Month==5 & p2016.data$EM3 == 1,1], Q$points[p2016.data$Month==5 & p2016.data$EM3 == 1,2], pch=10,col=my.cols[1],bg=my.cols[1], cex = 2)

#June
points(Q$points[p2016.data$Month==6 & p2016.data$EM3 == 1,1], Q$points[p2016.data$Month==6 & p2016.data$EM3 == 1,2], pch=10,col=my.cols[2],bg=my.cols[2], cex = 2)
points(Q$points[p2016.data$Month==6 & p2016.data$EM3 == 2,1], Q$points[p2016.data$Month==6 & p2016.data$EM3 == 2,2], pch=13,col=my.cols[2],bg=my.cols[2], cex = 2)

#July
points(Q$points[p2016.data$Month==7 & p2016.data$EM3 == 2,1], Q$points[p2016.data$Month==7 & p2016.data$EM3 == 2,2], pch=13,col=my.cols[3],bg=my.cols[3], cex = 2)
points(Q$points[p2016.data$Month==7 & p2016.data$EM3 == 3,1], Q$points[p2016.data$Month==7 & p2016.data$EM3 == 3,2], pch=19,col=my.cols[3],bg=my.cols[3], cex = 2)

#August
points(Q$points[p2016.data$Month==8 & p2016.data$EM3 == 3,1], Q$points[p2016.data$Month==8 & p2016.data$EM3 == 3,2], pch=19,col=my.cols[4],bg=my.cols[4], cex = 2)

#September
points(Q$points[p2016.data$Month==9 & p2016.data$EM3 == 3,1], Q$points[p2016.data$Month==9 & p2016.data$EM3 == 3,2], pch=19,col=my.cols[5],bg=my.cols[5], cex = 2)

#ordiarrows (Q, groups = p2016.data$Year, order.by = p2016.data$Date, startmark = "+", label = FALSE, length = .1)
#text(Y.sel[,1], Y.sel[,2], rownames(Y.sel), col="black", cex = 0.8)

legend("bottomleft",legend = c("pre-EM","post-EM 1","post-EM 2","post-EM 3"), pch = c(1,10,13,19), bty = "n", pt.cex = 2)
#legend("bottomright",legend = c("May","June","July","August","September"), pch = 22, pt.bg = my.cols, bty = "n", pt.cex = 2)
legend("topright",legend = c("k = 0.13"),bty = "n")
legend(-1.8,2.2,legend = c("A. 2016"),bty = "n", cex = 1.5)

#### 2017 ####

#select just 2017 data
p2017.data <- phytos1 %>%
  filter(year(Date) == 2017)

#prep data for NMDS input
p2017 <- p2017.data %>%
  select(-Year,-Date,-EM1,-EM2,-EM3,-Month,-MonthDay)
p2017 <- as.matrix(p2017)


#run NMDS for 2017
set.seed(1)
Q <- metaMDS(p2017, distance='bray', k=2, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
# plot(Q, display=c('sites'), choices=c(1,2), type='p')



#coded by month
my.cols <- gg_color_hue(5)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')

#May
points(Q$points[p2017.data$Month==5 & p2017.data$EM3 == 0,1], Q$points[p2017.data$Month==5 & p2017.data$EM3 == 0,2], pch=1,col=my.cols[1],bg=my.cols[1], cex = 2)
points(Q$points[p2017.data$Month==5 & p2017.data$EM3 == 1,1], Q$points[p2017.data$Month==5 & p2017.data$EM3 == 1,2], pch=10,col=my.cols[1],bg=my.cols[1], cex = 2)

#June
points(Q$points[p2017.data$Month==6 & p2017.data$EM3 == 1,1], Q$points[p2017.data$Month==6 & p2017.data$EM3 == 1,2], pch=10,col=my.cols[2],bg=my.cols[2], cex = 2)

#July
points(Q$points[p2017.data$Month==7 & p2017.data$EM3 == 1,1], Q$points[p2017.data$Month==7 & p2017.data$EM3 == 1,2], pch=10,col=my.cols[3],bg=my.cols[3], cex = 2)
points(Q$points[p2017.data$Month==7 & p2017.data$EM3 == 2,1], Q$points[p2017.data$Month==7 & p2017.data$EM3 == 2,2], pch=13,col=my.cols[3],bg=my.cols[3], cex = 2)

#August
points(Q$points[p2017.data$Month==8 & p2017.data$EM3 == 2,1], Q$points[p2017.data$Month==8 & p2017.data$EM3 == 2,2], pch=13,col=my.cols[4],bg=my.cols[4], cex = 2)

#September
points(Q$points[p2017.data$Month==9 & p2017.data$EM3 == 2,1], Q$points[p2017.data$Month==9 & p2017.data$EM3 == 2,2], pch=13,col=my.cols[5],bg=my.cols[5], cex = 2)

#ordiarrows (Q, groups = p2017.data$Year, order.by = p2017.data$Date, startmark = "+", label = FALSE, length = .1)
#text(Y.sel[,1], Y.sel[,2], rownames(Y.sel), col="black", cex = 0.8)

#legend("bottomleft",legend = c("pre-EM","post-EM 1","post-EM 2"), pch = c(1,10,13), bty = "n", pt.cex = 2)
#legend("bottomright",legend = c("May","June","July","August","September"), pch = 22, pt.bg = my.cols, bty = "n", pt.cex = 2)
legend("topright",legend = c("k = 0.08"),bty = "n")
legend(-1.6,1.45,legend = c("B. 2017"),bty = "n", cex = 1.5)

#### 2018 ####

#subset just 2018 data
p2018.data <- phytos1 %>%
  filter(year(Date) == 2018)

#prep data for NMDS input
p2018 <- p2018.data %>%
  select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay,-Year)
p2018 <- as.matrix(p2018)


#run NMDS for all months
set.seed(1)
Q <- metaMDS(p2018, distance='bray', k=2, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
# plot(Q, display=c('sites'), choices=c(1,2), type='p')



#coded by month
my.cols <- gg_color_hue(5)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')

#May
points(Q$points[p2018.data$Month==5 ,1], Q$points[p2018.data$Month==5 ,2], pch=1,col=my.cols[1], cex = 2)

#June
points(Q$points[p2018.data$Month==6 ,1], Q$points[p2018.data$Month==6 ,2], pch=1,col=my.cols[2], cex = 2)

#July
points(Q$points[p2018.data$Month==7 ,1], Q$points[p2018.data$Month==7 ,2], pch=1,col=my.cols[3], cex = 2)

#August
points(Q$points[p2018.data$Month==8 ,1], Q$points[p2018.data$Month==8 ,2], pch=1,col=my.cols[4], cex = 2)

#September
points(Q$points[p2018.data$Month==9 ,1], Q$points[p2018.data$Month==9 ,2], pch=1,col=my.cols[5], cex = 2)


legend("bottomleft",legend = c("May","June","July","August","September"), pch = 22, pt.bg = my.cols, bty = "n", pt.cex = 2)
legend("topright",legend = c("k = 0.12"),bty = "n")
legend(-2.6,1.5,legend = c("C. 2018"),bty = "n", cex = 1.5)



#### 2019 ####

#subset 2019 data only
p2019.data <- phytos1 %>%
  filter(year(Date) == 2019) %>%
  mutate(Storm = ifelse(Date <= "2019-06-08",1,2))

#prep data for NMDS input
p2019 <- p2019.data %>%
  select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay,-Year,-Storm)
p2019 <- as.matrix(p2019)


#run NMDS for all months
set.seed(1)
Q <- metaMDS(p2019, distance='bray', k=2, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
# plot(Q, display=c('sites'), choices=c(1,2), type='p')



#coded by month
my.cols <- gg_color_hue(5)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')

#May
points(Q$points[p2019.data$Month==5 ,1], Q$points[p2019.data$Month==5 ,2], pch=3,col=my.cols[1], cex = 2)

#June
points(Q$points[p2019.data$Month==6 & p2019.data$Storm == 1,1], Q$points[p2019.data$Month==6 & p2019.data$Storm == 1,2], pch=3,col=my.cols[2], cex = 2)
points(Q$points[p2019.data$Month==6 & p2019.data$Storm == 2,1], Q$points[p2019.data$Month==6 & p2019.data$Storm == 2,2], pch=8,col=my.cols[2], cex = 2)

#July
points(Q$points[p2019.data$Month==7 ,1], Q$points[p2019.data$Month==7 ,2], pch=8,col=my.cols[3], cex = 2)

#August
points(Q$points[p2019.data$Month==8 ,1], Q$points[p2019.data$Month==8 ,2], pch=8,col=my.cols[4], cex = 2)

#September
points(Q$points[p2019.data$Month==9 ,1], Q$points[p2019.data$Month==9 ,2], pch=8,col=my.cols[5], cex = 2)


#legend("bottomright",legend = c("May","June","July","August","September"), pch = 22, pt.bg = my.cols, bty = "n", pt.cex = 2)
legend("bottomleft",legend = c("pre-storm","post-storm"), pch = c(3,8), bty = "n", pt.cex = 2)
legend("topright",legend = c("k = 0.04"),bty = "n")
legend(-2.25,1.75,legend = c("D. 2019"),bty = "n", cex = 1.5)

dev.off()

