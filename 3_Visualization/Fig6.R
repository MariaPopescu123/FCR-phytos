#Fig6
#Author: Mary Lofton
#Date: 23FEB21

####SET-UP####
#load packages
pacman::p_load(tidyverse, lubridate, vegan, Hmisc, picante, indicspecies, vegan3d, cowplot, viridis)
rm(list=ls())

#set functions
#get colors for plotting
gg_color_hue <- function(n) {
  hues = seq(15, 450, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#select variables from envfit

#Function: select.envfit - Setting r2 cutoff values to display in an ordination.r.select<-0.3 # correlation threshold, see function below
#__FUNCTION: select.envfit__#
# function (select.envit) filters the resulting list of function (envfit) based on their p values. This allows to display only significant values in the final plot.
# just run this
select.envfit<-function(fit, name.list){ #needs two sorts of input: fit= result of envfit, r.select= numeric, correlation minimum threshold
  for (i in 1:length(row.names(fit$vectors$arrows))) { #run for-loop through the entire length of the column r in object fit$vectors$r starting at i=1
    myvars <- row.names(fit$vectors$arrows)
    if (myvars[i] %in% name.list) { #Check wether r<r.select, i.e. if the correlation is weaker than the threshold value. Change this Parameter for r-based selection
      fit$vectors$arrows[i,]=NA #If the above statement is TRUE, i.e. r is smaller than r.select, then the coordinates of the vectors are set to NA, so they cannot be displayed
      i=i+1 #increase the running parameter i from 1 to 2, i.e. check the next value in the column until every value has been checked
    } #close if-loop
  } #close for-loop
  return(fit) #return fit as the result of the function
} #close the function


#read in environmental data
my.cs.data <- read_csv("./2_Data_analysis/CS_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_DO_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") %>%
  select(-MonthDay) %>%
  filter(!is.na(BV_TOTAL)) 

colnames(my.cs.data)
env <- my.cs.data[,c(2,5:6,11:13,27,29:32,39,40:41,43)]
env$Month <- month(my.cs.data$Date)
env <- data.frame(scale(env, center = TRUE, scale = TRUE)) 

#read in fp data
my.fp.data <- read_csv("./2_Data_analysis/FP_megamatrix.csv")

#read in community data
phytos <- read_csv("./0_Data_files/phytoplankton.csv") %>%
  select(Date, Genus, BV_um3mL) %>%
  spread(key = Genus, value = BV_um3mL)

# any non-detects filled with 0
phytos[is.na(phytos)]<- 0

# create data frame with different EM periods for plotting later on
em <- my.cs.data[,c(1,2,4,5,6)] %>%
  mutate(Month = month(Date))

phytos1 <- left_join(em, phytos, by = "Date") %>%
  mutate(MonthDay = format(Date, format="%m-%d"))

fp <- left_join(phytos1,my.fp.data, by = "Date") %>%
  select(Date, Peak_depth_m)

# prepare phyto data for NMDS input
phytos2 <- phytos1 %>%
  select(-Date,-EM1,-EM2,-EM3,-Month,-Year,-MonthDay) 

phytos3 <- as.matrix(phytos2)


####RUN THE NMDS FOR ALL YEARS####


#based on results of scree plot, the best choice for k is 3

#run NMDS for all months
set.seed(1)
Q <- metaMDS(phytos3, distance='bray', k=3, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
plot(Q, display=c('sites'), choices=c(1,2), type='p')
plot(Q, display=c('sites'), choices=c(1,3), type='p')
# Q$species
# Q$points

# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
# significance level < 0.01
env_all_12 <- env %>%
  select(Year, EM2, EM3, Temp_C_grab, n2,lake.number, wedderburn.number, Max_biomass_ugL, Peak_depth_m, Month)

env_all_13 <- env %>%
  select(EM3, Temp_C_grab, wedderburn.number, perc_light_grab, Max_biomass_ugL, Peak_depth_m, Month)

env_all <- env %>%
  select(Year, EM2, EM3, Temp_C_grab, n2,lake.number, wedderburn.number, perc_light_grab, Max_biomass_ugL, Peak_depth_m, Month)

#run final correlations w/ significant environmental variables
set.seed(1)
en12 = envfit(Q, env_all_12, permutations = 9999, na.rm = TRUE, choices = c(1,2))
en12
row.names(en12$vectors$arrows) <- c("Yr","D1/PD","D2","T","N2","LN","WN","Biom","","Mo")

set.seed(1)
en13 = envfit(Q, env_all_13, permutations = 9999, na.rm = TRUE, choices = c(1,3))
en13
row.names(en13$vectors$arrows) <- c("D2","T","","%PAR","Biom","PD","Mo")

en12_notEM <- select.envfit(fit = en12, name.list = c("D1/PD","D2","N2",""))
en12_EM <- select.envfit(fit = en12, name.list = c("Yr","T","LN","WN","Biom","Mo"))

en13_notEM <- select.envfit(fit = en13, name.list = c("D2","PD"))
en13_EM <- select.envfit(fit = en13, name.list = c("T","","%PAR","Biom","Mo"))


#first and second axes by EM2
tiff(filename = "./3_Visualization/Fig6.tif",width = 6.5,height = 9.75,units = "in",res = 300)
par(mfrow = c(3,2),mar = c(3,3,1.5,0),mgp = c(1.5,0.5,0))
plot(Q, display=c('sites'),choices=c(1,2), type='n')
my.cols <- viridis(5)
#2016
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 5,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 5,2], pch=10,col = my.cols[1],bg=my.cols[1], cex = 2)
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 6,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 6,2], pch=10,col = my.cols[2],bg=my.cols[2], cex = 2)
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 7,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 7,2], pch=10,col = my.cols[3],bg=my.cols[3], cex = 2)
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 8,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 8,2], pch=10,col = my.cols[4],bg=my.cols[4], cex = 2)
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 9,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 9,2], pch=10,col = my.cols[5],bg=my.cols[5], cex = 2)

#2017
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 5,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 5,2], pch=8,col = my.cols[1],bg=my.cols[1], cex = 2)
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 6,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 6,2], pch=8,col = my.cols[2],bg=my.cols[2], cex = 2)
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 7,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 7,2], pch=8,col = my.cols[3],bg=my.cols[3], cex = 2)
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 8,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 8,2], pch=8,col = my.cols[4],bg=my.cols[4], cex = 2)
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 9,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 9,2], pch=8,col = my.cols[5],bg=my.cols[5], cex = 2)

#2018
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 5,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 5,2], pch=23,col = "black",bg=my.cols[1],cex = 2)
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 6,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 6,2], pch=23,col = "black",bg=my.cols[2], cex = 2)
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 7,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 7,2], pch=23,col = "black",bg=my.cols[3], cex = 2)
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 8,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 8,2], pch=23,col = "black",bg=my.cols[4], cex = 2)
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 9,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 9,2], pch=23,col = "black",bg=my.cols[5], cex = 2)

#2019
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 5,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 5,2], pch=25,col = "black",bg=my.cols[1],cex = 2)
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 6,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 6,2], pch=25,col = "black",bg=my.cols[2], cex = 2)
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 7,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 7,2], pch=25,col = "black",bg=my.cols[3], cex = 2)
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 8,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 8,2], pch=25,col = "black",bg=my.cols[4], cex = 2)
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 9,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 9,2], pch=25,col = "black",bg=my.cols[5], cex = 2)

plot(en12_notEM, col = "gray48")
plot(en12_EM, col = "black")
legend("bottomleft",legend = c("2016","2017","2018","2019"), pch = c(10,8,23,25), pt.bg = c(NA,NA,"black","black"), bty = "n", pt.cex = 2)
#legend(x = -2, y = -0.2,legend = c("May","June","July","August","September"), pch = 22, pt.bg = my.cols, bty = "n", pt.cex = 2)
legend("topright",legend = c("k = 0.11"),bty = "n")
title("A. all summers",adj = 0, line = 0.5)

plot(Q, display=c('sites'),choices=c(1,3), type='n')
my.cols <- viridis(5)
#2016
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 5,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 5,3], pch=10,col = my.cols[1],bg=my.cols[1], cex = 2)
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 6,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 6,3], pch=10,col = my.cols[2],bg=my.cols[2], cex = 2)
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 7,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 7,3], pch=10,col = my.cols[3],bg=my.cols[3], cex = 2)
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 8,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 8,3], pch=10,col = my.cols[4],bg=my.cols[4], cex = 2)
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 9,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2016 & phytos1$Month == 9,3], pch=10,col = my.cols[5],bg=my.cols[5], cex = 2)

#2017
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 5,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 5,3], pch=8,col = my.cols[1],bg=my.cols[1], cex = 2)
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 6,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 6,3], pch=8,col = my.cols[2],bg=my.cols[2], cex = 2)
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 7,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 7,3], pch=8,col = my.cols[3],bg=my.cols[3], cex = 2)
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 8,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 8,3], pch=8,col = my.cols[4],bg=my.cols[4], cex = 2)
points(Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 9,1], Q$points[phytos1$EM2==1 & phytos1$Year == 2017 & phytos1$Month == 9,3], pch=8,col = my.cols[5],bg=my.cols[5], cex = 2)

#2018
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 5,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 5,3], pch=23,col = "black",bg=my.cols[1],cex = 2)
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 6,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 6,3], pch=23,col = "black",bg=my.cols[2], cex = 2)
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 7,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 7,3], pch=23,col = "black",bg=my.cols[3], cex = 2)
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 8,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 8,3], pch=23,col = "black",bg=my.cols[4], cex = 2)
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 9,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2018 & phytos1$Month == 9,3], pch=23,col = "black",bg=my.cols[5], cex = 2)

#2019
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 5,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 5,3], pch=25,col = "black",bg=my.cols[1],cex = 2)
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 6,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 6,3], pch=25,col = "black",bg=my.cols[2], cex = 2)
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 7,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 7,3], pch=25,col = "black",bg=my.cols[3], cex = 2)
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 8,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 8,3], pch=25,col = "black",bg=my.cols[4], cex = 2)
points(Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 9,1], Q$points[phytos1$EM2==0 & phytos1$Year == 2019 & phytos1$Month == 9,3], pch=25,col = "black",bg=my.cols[5], cex = 2)

plot(en13_notEM, col = "gray48")
plot(en13_EM, col = "black")
text(1, 0.45, "WN",col = "gray48")
legend("bottomleft",legend = c("2016","2017","2018","2019"), pch = c(10,8,23,25), pt.bg = c(NA,NA,"black","black"), bty = "n", pt.cex = 2)
#legend("bottomleft",legend = c("2016 (EM)","2017 (EM)","2018 (no EM)","2019 (no EM)"), pch = c(10,8,23,25), pt.bg = c(NA,NA,"black","black"), bty = "n", pt.cex = 2)
#legend("bottomright",legend = c("May","June","July","August","September"), pch = 22, pt.bg = my.cols, bty = "n", pt.cex = 2)
legend("topright",legend = c("k = 0.11"),bty = "n")
title("B. all summers",adj = 0, line = 0.5)

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
phytos <- read_csv("./0_Data_files/phytoplankton.csv") %>%
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
my.cols <- viridis(5)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
title("C. 2016",adj = 0, line = 0.5)


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
my.cols <- viridis(5)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
title("D. 2017",adj = 0, line = 0.5)


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

legend("bottomleft",legend = c("pre-EM","post-EM 1","post-EM 2"), pch = c(1,10,13), bty = "n", pt.cex = 2)
#legend("bottomright",legend = c("May","June","July","August","September"), pch = 22, pt.bg = my.cols, bty = "n", pt.cex = 2)
legend("topright",legend = c("k = 0.08"),bty = "n")

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
my.cols <- viridis(5)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
title("E. 2018",adj = 0, line = 0.5)


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
my.cols <- viridis(5)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
title("F. 2019",adj = 0, line = 0.5)

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

dev.off()



