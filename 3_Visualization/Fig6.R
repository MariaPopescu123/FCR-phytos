#Fig. 4
#Author: Mary Lofton
#Date: 23FEB21

####SET-UP####
#load packages
pacman::p_load(tidyverse, lubridate, vegan, Hmisc, picante, indicspecies, vegan3d, cowplot)
rm(list=ls())

#set functions
#get colors for plotting
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
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
env <- my.cs.data[,c(2,5:6,11:13,27,29:30,38,39:40,42)]
env$Month <- month(my.cs.data$Date)
env <- data.frame(scale(env, center = TRUE, scale = TRUE)) 

#read in fp data
my.fp.data <- read_csv("./2_Data_analysis/FP_megamatrix.csv")

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

fp <- left_join(phytos1,my.fp.data, by = "Date") %>%
  select(Date, Peak_depth_m)

# prepare phyto data for NMDS input
phytos2 <- phytos1 %>%
  select(-Date,-EM1,-EM2,-EM3,-Month,-Year,-MonthDay) 

phytos3 <- as.matrix(phytos2)


####RUN THE NMDS FOR ALL YEARS####

# #calculate Hellinger distance - FOR COMMUNITY DATA
# WHY WOULD I DO THIS??
# hellinger_vars <- decostand(phytos[,-c(1:3)], method = "hellinger")


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
  select(Year, EM2, EM3, Temp_C_grab, n2,Max_biomass_ugL, Peak_depth_m, Month)

env_all_13 <- env %>%
  select(EM3, Temp_C_grab, perc_light_grab, Max_biomass_ugL, Peak_depth_m, Month)

env_all <- env %>%
  select(Year, EM2, EM3, Temp_C_grab, n2,perc_light_grab, Max_biomass_ugL, Peak_depth_m, Month)

#run final correlations w/ significant environmental variables
set.seed(1)
en12 = envfit(Q, env_all_12, permutations = 9999, na.rm = TRUE, choices = c(1,2))
en12
row.names(en12$vectors$arrows) <- c("Yr","DEEP1/PD","DEEP2","T","N2","Biom","","Mo")

set.seed(1)
en13 = envfit(Q, env_all_13, permutations = 9999, na.rm = TRUE, choices = c(1,3))
en13
row.names(en13$vectors$arrows) <- c("DEEP2","T","%PAR","Biom","PD","Mo")

en12_notEM <- select.envfit(fit = en12, name.list = c("DEEP1/PD","EM2","N2",""))
en12_EM <- select.envfit(fit = en12, name.list = c("Yr","T","Biom","Mo"))

en13_notEM <- select.envfit(fit = en13, name.list = c("DEEP2","PD"))
en13_EM <- select.envfit(fit = en13, name.list = c("T","%PAR","Biom","Mo"))


#first and second axes by EM2
tiff(filename = "./3_Visualization/Fig4.tif",width = 12.4,height = 6.2,units = "in",res = 300)
par(mfrow = c(1,2))
plot(Q, display=c('sites'),choices=c(1,2), type='n')
my.cols = gg_color_hue(5)
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
legend(x = -2, y = -0.2,legend = c("May","June","July","August","September"), pch = 22, pt.bg = my.cols, bty = "n", pt.cex = 2)
legend("topright",legend = c("k = 0.11"),bty = "n")
title("A",adj = 0, line = 0.5)

plot(Q, display=c('sites'),choices=c(1,3), type='n')
my.cols = gg_color_hue(5)
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
#legend("bottomleft",legend = c("2016 (EM)","2017 (EM)","2018 (no EM)","2019 (no EM)"), pch = c(10,8,23,25), pt.bg = c(NA,NA,"black","black"), bty = "n", pt.cex = 2)
#legend("bottomright",legend = c("May","June","July","August","September"), pch = 22, pt.bg = my.cols, bty = "n", pt.cex = 2)
#legend("topright",legend = c("k = 0.11"),bty = "n")
title("B",adj = 0, line = 0.5)
dev.off()

