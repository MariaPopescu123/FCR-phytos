#Title: NMDS for FCR phytos
#Author: Mary Lofton
#Date: 11DEC2017
#Updated: 01FEB2021

####SET-UP####
#load packages
pacman::p_load(tidyverse, lubridate, vegan, Hmisc, picante)
rm(list=ls())

#read in environmental data
my.cs.data <- read_csv("./2_Data_analysis/CS_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") %>%
  select(-MonthDay) %>%
  filter(!is.na(BV_TOTAL))

colnames(my.cs.data)
env <- my.cs.data[,c(2:5,7:20)]
env <- data.frame(scale(env, center = TRUE, scale = TRUE)) 

#read in community data
phytos <- read_csv("./00_Data_files/EDI_phytos/phytoplankton.csv") %>%
  select(Date, Genus, BV_um3mL) %>%
  spread(key = Genus, value = BV_um3mL)

# any non-detects filled with 0
phytos[is.na(phytos)]<- 0

# create data frame with different EM periods for plotting later on
em <- my.cs.data[,c(1,4,5,6)] %>%
  mutate(Month = month(Date))

phytos1 <- left_join(em, phytos, by = "Date")

# prepare phyto data for NMDS input
phytos2 <- phytos1 %>%
  select(-Date,-EM1,-EM2,-EM3,-Month)

phytos3 <- as.matrix(phytos2)


####RUN THE NMDS FOR ALL YEARS####

# #calculate Hellinger distance - FOR COMMUNITY DATA
# WHY WOULD I DO THIS??
# hellinger_vars <- decostand(phytos[,-c(1:3)], method = "hellinger")

##scree plot
stress <- rep(NA,6)
k <- c(1:6)
for (i in 1:length(k)){
  Q <- metaMDS(phytos3, distance='bray', k=i, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
  stress[i] <- Q$stress
}
plot(k, stress, type='b', xlab='Dimensions', ylab='Stress', pch=21, bg='grey75', cex=1.3)
legend("topright",legend = c("all years"),bty = 'n')

#based on results of scree plot, the best choice for k is 3

#run NMDS for all months
Q <- metaMDS(phytos3, distance='bray', k=3, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
plot(Q, display=c('sites'), choices=c(1,2), type='p')
plot(Q, display=c('sites'), choices=c(1,3), type='p')

#run initial correlations with environmental variables
en_0 = envfit(Q, env, permutations = 999, na.rm = TRUE)
en_0

# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
env_all <- env %>%
  select(Year, EM2, DINmax_ugL, DOCmax_depth_m, DOCmax_mgL,
         thermo.depth, Kd, perc_light_thermocline)

#run final correlations w/ significant environmental variables
en = envfit(Q, env_all, permutations = 999, na.rm = TRUE)
en

#get colors for plotting
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
my.cols <- gg_color_hue(14)


#plot results with EM periods and relationships w/ environmental variables

#first and second axes
png(filename = "./3_Visualization/NMDS_all_1_2_EM3.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[phytos1$EM3=="pre-mix 2016",1], Q$points[phytos1$EM3=="pre-mix 2016",2], pch=21,bg=my.cols[1])
points(Q$points[phytos1$EM3=="post-EM1 2016",1], Q$points[phytos1$EM3=="post-EM1 2016",2], pch=21,bg=my.cols[2])
points(Q$points[phytos1$EM3=="pre-EM2 2016",1], Q$points[phytos1$EM3=="pre-EM2 2016",2], pch=21,bg=my.cols[3])
points(Q$points[phytos1$EM3=="post-EM2 2016",1], Q$points[phytos1$EM3=="post-EM2 2016",2], pch=21,bg=my.cols[4])
points(Q$points[phytos1$EM3=="pre-EM3 2016",1], Q$points[phytos1$EM3=="pre-EM3 2016",2], pch=21,bg=my.cols[5])
points(Q$points[phytos1$EM3=="post-EM3 2016",1], Q$points[phytos1$EM3=="post-EM3 2016",2], pch=21,bg=my.cols[6])
points(Q$points[phytos1$EM3=="post-mix 2016",1], Q$points[phytos1$EM3=="post-mix 2016",2], pch=21,bg=my.cols[7])
points(Q$points[phytos1$EM3=="pre-mix 2017",1], Q$points[phytos1$EM3=="pre-mix 2017",2], pch=21,bg=my.cols[8])
points(Q$points[phytos1$EM3=="post-EM1 2017",1], Q$points[phytos1$EM3=="post-EM1 2017",2], pch=21,bg=my.cols[9])
points(Q$points[phytos1$EM3=="pre-EM2 2017",1], Q$points[phytos1$EM3=="pre-EM2 2017",2], pch=21,bg=my.cols[10])
points(Q$points[phytos1$EM3=="post-EM2 2017",1], Q$points[phytos1$EM3=="post-EM2 2017",2], pch=21,bg=my.cols[11])
points(Q$points[phytos1$EM3=="post-mix 2017",1], Q$points[phytos1$EM3=="post-mix 2017",2], pch=21,bg=my.cols[12])
points(Q$points[phytos1$EM3=="2018",1], Q$points[phytos1$EM3=="2018",2], pch=21,bg=my.cols[13])
points(Q$points[phytos1$EM3=="2019",1], Q$points[phytos1$EM3=="2019",2], pch=21,bg=my.cols[14])
plot(en)
legend("topright",legend = c("pre-mix 2016","post-EM1 2016","pre-EM2 2016",
                             "post-EM2 2016","pre-EM3 2016","post-EM3 2016",
                             "post-mix 2016","pre-mix 2017","post-EM1 2017",
                             "pre-EM2 2017","post-EM2 2017","post-mix 2017",
                             "2018","2019"), pch = 21, pt.bg = my.cols, bty = "n")
legend("bottomright",legend = c("k = 0.11"),bty = "n")
dev.off()

#first two axes in black and white
png(filename = "./3_Visualization/NMDS_all_1_2_EM1.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[phytos1$EM1==1,1], Q$points[phytos1$EM1==1,2], pch=21,bg="black")
points(Q$points[phytos1$EM1==0,1], Q$points[phytos1$EM1==0,2], pch=21,bg="white")
plot(en)
legend("topright",legend = c("EM + 2 wks post","no EM"), pch = 21, pt.bg = c("black","white"), bty = "n")
legend("bottomright",legend = c("k = 0.11"),bty = "n")
dev.off()

#first two axes in black and white w/ EM coded by year
png(filename = "./3_Visualization/NMDS_all_1_2_EM2.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[phytos1$EM2==1,1], Q$points[phytos1$EM2==1,2], pch=21,bg="black")
points(Q$points[phytos1$EM2==0,1], Q$points[phytos1$EM2==0,2], pch=21,bg="white")
plot(en)
legend("topright",legend = c("EM years","years w/ no EM"), pch = 21, pt.bg = c("black","white"), bty = "n")
legend("bottomright",legend = c("k = 0.11"),bty = "n")
dev.off()

#first and third axes
png(filename = "./3_Visualization/NMDS_all_1_3_EM3.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,3), type='n')
points(Q$points[phytos1$EM3=="pre-mix 2016",1], Q$points[phytos1$EM3=="pre-mix 2016",3], pch=21,bg=my.cols[1])
points(Q$points[phytos1$EM3=="post-EM1 2016",1], Q$points[phytos1$EM3=="post-EM1 2016",3], pch=21,bg=my.cols[2])
points(Q$points[phytos1$EM3=="pre-EM2 2016",1], Q$points[phytos1$EM3=="pre-EM2 2016",3], pch=21,bg=my.cols[3])
points(Q$points[phytos1$EM3=="post-EM2 2016",1], Q$points[phytos1$EM3=="post-EM2 2016",3], pch=21,bg=my.cols[4])
points(Q$points[phytos1$EM3=="pre-EM3 2016",1], Q$points[phytos1$EM3=="pre-EM3 2016",3], pch=21,bg=my.cols[5])
points(Q$points[phytos1$EM3=="post-EM3 2016",1], Q$points[phytos1$EM3=="post-EM3 2016",3], pch=21,bg=my.cols[6])
points(Q$points[phytos1$EM3=="post-mix 2016",1], Q$points[phytos1$EM3=="post-mix 2016",3], pch=21,bg=my.cols[7])
points(Q$points[phytos1$EM3=="pre-mix 2017",1], Q$points[phytos1$EM3=="pre-mix 2017",3], pch=21,bg=my.cols[8])
points(Q$points[phytos1$EM3=="post-EM1 2017",1], Q$points[phytos1$EM3=="post-EM1 2017",3], pch=21,bg=my.cols[9])
points(Q$points[phytos1$EM3=="pre-EM2 2017",1], Q$points[phytos1$EM3=="pre-EM2 2017",3], pch=21,bg=my.cols[10])
points(Q$points[phytos1$EM3=="post-EM2 2017",1], Q$points[phytos1$EM3=="post-EM2 2017",3], pch=21,bg=my.cols[11])
points(Q$points[phytos1$EM3=="post-mix 2017",1], Q$points[phytos1$EM3=="post-mix 2017",3], pch=21,bg=my.cols[12])
points(Q$points[phytos1$EM3=="2018",1], Q$points[phytos1$EM3=="2018",3], pch=21,bg=my.cols[13])
points(Q$points[phytos1$EM3=="2019",1], Q$points[phytos1$EM3=="2019",3], pch=21,bg=my.cols[14])
plot(en)
legend("topright",legend = c("pre-mix 2016","post-EM1 2016","pre-EM2 2016",
                             "post-EM2 2016","pre-EM3 2016","post-EM3 2016",
                             "post-mix 2016","pre-mix 2017","post-EM1 2017",
                             "pre-EM2 2017","post-EM2 2017","post-mix 2017",
                             "2018","2019"), pch = 21, pt.bg = my.cols, bty = "n")
legend("bottomright",legend = c("k = 0.11"),bty = "n")
dev.off()

#first and third axes in black and white
png(filename = "./3_Visualization/NMDS_all_1_3_EM1.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,3), type='n')
points(Q$points[phytos1$EM1==1,1], Q$points[phytos1$EM1==1,3], pch=21,bg="black")
points(Q$points[phytos1$EM1==0,1], Q$points[phytos1$EM1==0,3], pch=21,bg="white")
plot(en)
legend("topright",legend = c("EM + 2 wks post","no EM"), pch = 21, pt.bg = c("black","white"), bty = "n")
legend("bottomright",legend = c("k = 0.11"),bty = "n")
dev.off()

#first and third axes in black and white w/ EM coded by year
png(filename = "./3_Visualization/NMDS_all_1_3_EM2.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,3), type='n')
points(Q$points[phytos1$EM2==1,1], Q$points[phytos1$EM2==1,3], pch=21,bg="black")
points(Q$points[phytos1$EM2==0,1], Q$points[phytos1$EM2==0,3], pch=21,bg="white")
plot(en)
legend("topright",legend = c("EM years","years w/ no EM"), pch = 21, pt.bg = c("black","white"), bty = "n")
legend("bottomright",legend = c("k = 0.11"),bty = "n")
dev.off()


#plot results by month including relationships w/ environmental variables

#first and second axes
my.cols <- gg_color_hue(5)
png(filename = "./3_Visualization/NMDS_all_1_2_month.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[phytos1$Month==5,1], Q$points[phytos1$Month==5,2], pch=21,bg=my.cols[1])
points(Q$points[phytos1$Month==6,1], Q$points[phytos1$Month==6,2], pch=21,bg=my.cols[2])
points(Q$points[phytos1$Month==7,1], Q$points[phytos1$Month==7,2], pch=21,bg=my.cols[3])
points(Q$points[phytos1$Month==8,1], Q$points[phytos1$Month==8,2], pch=21,bg=my.cols[4])
points(Q$points[phytos1$Month==9,1], Q$points[phytos1$Month==9,2], pch=21,bg=my.cols[5])
plot(en)
legend("topright", legend=c('May','June', 'July','Aug','Sept'), pch=21, pt.bg=c(my.cols),bty = "n") 
legend("bottomright",legend = c("k = 0.11"),bty = "n")
dev.off()

#first and third axes
png(filename = "./3_Visualization/NMDS_all_1_3_month.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,3), type='n')
points(Q$points[phytos1$Month==5,1], Q$points[phytos1$Month==5,3], pch=21,bg=my.cols[1])
points(Q$points[phytos1$Month==6,1], Q$points[phytos1$Month==6,3], pch=21,bg=my.cols[2])
points(Q$points[phytos1$Month==7,1], Q$points[phytos1$Month==7,3], pch=21,bg=my.cols[3])
points(Q$points[phytos1$Month==8,1], Q$points[phytos1$Month==8,3], pch=21,bg=my.cols[4])
points(Q$points[phytos1$Month==9,1], Q$points[phytos1$Month==9,3], pch=21,bg=my.cols[5])
plot(en)
legend("topright", legend=c('May','June', 'July','Aug','Sept'), pch=21, pt.bg=c(my.cols),bty = "n") 
legend("bottomright",legend = c("k = 0.11"),bty = "n")
dev.off()

#plot w/ genera overlain - this is just for preliminary interpretation, 
#these will be difficult to read
png(filename = "./3_Visualization/NMDS_all_1_2_genera.png",width = 9,height = 9,units = "in",res = 300)
fig <- ordiplot(Q, choices = c(1,2))
points(Q$points[phytos1$EM1==1,1], Q$points[phytos1$EM1==1,2], pch=21,bg='black')
text(fig, "species", col="blue", cex=0.6)
legend("topright",legend = c("EM + 2 wks post","no EM"), pch = 21, pt.bg = c("black","white"), bty = "n")
dev.off()

png(filename = "./3_Visualization/NMDS_all_1_3_genera.png",width = 9,height = 9,units = "in",res = 300)
fig <- ordiplot(Q, choices = c(1,3))
points(Q$points[phytos1$EM1==1,1], Q$points[phytos1$EM1==1,3], pch=21,bg='black')
text(fig, "species", col="blue", cex=0.6)
legend("topright",legend = c("EM + 2 wks post","no EM"), pch = 21, pt.bg = c("black","white"), bty = "n")
dev.off()


#### NMDS FOR 2016 ONLY ####

#select 2016 data
p2016.data <- phytos1 %>%
  filter(year(Date) == 2016)

#prep data for NMDS input
p2016 <- p2016.data %>%
  select(-Date,-EM1,-EM2,-EM3,-Month)
p2016 <- as.matrix(p2016)

#select out environmental variables that are constant in 2016
env2016 <- env[1:20,] %>%
  select(-Year, -HOx, -EM2)

##scree plot
#  Basic outline for making your own Scree plot
stress <- rep(NA,6)
k <- c(1:6)
for (i in 1:length(k)){
  Q <- metaMDS(p2016, distance='bray', k=i, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
  stress[i] <- Q$stress
}
plot(k, stress, type='b', xlab='Dimensions', ylab='Stress', pch=21, bg='grey75', cex=1.3)
legend("topright",legend = c("2016"),bty = 'n')

#based on scree plot results, the best value for k is 2

#run NMDS for all months
Q <- metaMDS(p2016, distance='bray', k=2, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
plot(Q, display=c('sites'), choices=c(1,2), type='p')

#get initial environmental relationships
en_0 = envfit(Q, env2016, permutations = 999, na.rm = TRUE)
en_0

# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
env_2016 <- env2016 %>%
  select(DINmax_depth_m, Temp_C, 
         thermo.depth, perc_light_thermocline)

#run final correlations w/ significant environmental variables
en = envfit(Q, env_2016, permutations = 999, na.rm = TRUE)
en

#get colors for plots
my.cols <- gg_color_hue(7)

#plot results
png(filename = "./3_Visualization/NMDS_2016_1_2_EM3.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[p2016.data$EM3=="pre-mix 2016",1], Q$points[p2016.data$EM3=="pre-mix 2016",2], pch=21,bg=my.cols[1])
points(Q$points[p2016.data$EM3=="post-EM1 2016",1], Q$points[p2016.data$EM3=="post-EM1 2016",2], pch=21,bg=my.cols[2])
points(Q$points[p2016.data$EM3=="pre-EM2 2016",1], Q$points[p2016.data$EM3=="pre-EM2 2016",2], pch=21,bg=my.cols[3])
points(Q$points[p2016.data$EM3=="post-EM2 2016",1], Q$points[p2016.data$EM3=="post-EM2 2016",2], pch=21,bg=my.cols[4])
points(Q$points[p2016.data$EM3=="pre-EM3 2016",1], Q$points[p2016.data$EM3=="pre-EM3 2016",2], pch=21,bg=my.cols[5])
points(Q$points[p2016.data$EM3=="post-EM3 2016",1], Q$points[p2016.data$EM3=="post-EM3 2016",2], pch=21,bg=my.cols[6])
points(Q$points[p2016.data$EM3=="post-mix 2016",1], Q$points[p2016.data$EM3=="post-mix 2016",2], pch=21,bg=my.cols[7])
plot(en)
legend("topright",legend = c("pre-mix 2016","post-EM1 2016","pre-EM2 2016",
                             "post-EM2 2016","pre-EM3 2016","post-EM3 2016",
                             "post-mix 2016"), pch = 21, pt.bg = my.cols, bty = "n")
legend("bottomright",legend = c("k = 0.13"),bty = "n")
dev.off()

#in black and white
png(filename = "./3_Visualization/NMDS_2016_1_2_EM1.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[p2016.data$EM1==1,1], Q$points[p2016.data$EM1==1,2], pch=21,bg="black")
points(Q$points[p2016.data$EM1==0,1], Q$points[p2016.data$EM1==0,2], pch=21,bg="white")
plot(en)
legend("topright",legend = c("EM + 2 wks post","no EM"), pch = 21, pt.bg = c("black","white"), bty = "n")
legend("bottomright",legend = c("k = 0.13"),bty = "n")
dev.off()

#coded by month
my.cols <- gg_color_hue(5)
png(filename = "./3_Visualization/NMDS_2016_1_2_month.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[p2016.data$Month==5,1], Q$points[p2016.data$Month==5,2], pch=21,bg=my.cols[1])
points(Q$points[p2016.data$Month==6,1], Q$points[p2016.data$Month==6,2], pch=21,bg=my.cols[2])
points(Q$points[p2016.data$Month==7,1], Q$points[p2016.data$Month==7,2], pch=21,bg=my.cols[3])
points(Q$points[p2016.data$Month==8,1], Q$points[p2016.data$Month==8,2], pch=21,bg=my.cols[4])
points(Q$points[p2016.data$Month==9,1], Q$points[p2016.data$Month==9,2], pch=21,bg=my.cols[5])
plot(en)
legend("topright", legend=c('May','June', 'July','Aug','Sept'), pch=21, pt.bg=c(my.cols),bty = "n") 
legend("bottomright",legend = c("k = 0.13"),bty = "n")
dev.off()

#plot w/ genera overlain - this is just for preliminary interpretation, 
#these will be difficult to read
png(filename = "./3_Visualization/NMDS_2016_1_2_genera.png",width = 9,height = 9,units = "in",res = 300)
fig <- ordiplot(Q, choices = c(1,2))
points(Q$points[p2016.data$EM1==1,1], Q$points[p2016.data$EM1==1,2], pch=21,bg='black')
text(fig, "species", col="blue", cex=0.6)
legend("topright",legend = c("EM + 2 wks post","no EM"), pch = 21, pt.bg = c("black","white"), bty = "n")
dev.off()


#### 2017 ####

#select just 2017 data
p2017.data <- phytos1 %>%
  filter(year(Date) == 2017)

#prep data for NMDS input
p2017 <- p2017.data %>%
  select(-Date,-EM1,-EM2,-EM3,-Month)
p2017 <- as.matrix(p2017)

#select out environmental variables that are constant in 2017
env2017 <- env[21:35,] %>%
  select(-Year, -HOx, -EM2)

##scree plot
#  Basic outline for making your own Scree plot
stress <- rep(NA,6)
k <- c(1:6)
for (i in 1:length(k)){
  Q <- metaMDS(p2017, distance='bray', k=i, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
  stress[i] <- Q$stress
}
plot(k, stress, type='b', xlab='Dimensions', ylab='Stress', pch=21, bg='grey75', cex=1.3)
legend("topright",legend = c("2017"),bty = 'n')

#based on scree plot the best value for k is 2

#run NMDS for 2017
Q <- metaMDS(p2017, distance='bray', k=2, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
plot(Q, display=c('sites'), choices=c(1,2), type='p')

#get initial environmental relationships
en_0 = envfit(Q, env2017, permutations = 999, na.rm = TRUE)
en_0

# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
env_2017 <- env2017 %>%
  select(EM1, Cmax_DOC_mgL, DINmax_depth_m, DINmax_ugL)

#run final correlations w/ significant environmental variables
en = envfit(Q, env_2017, permutations = 999, na.rm = TRUE)
en

my.cols <- gg_color_hue(5)
#plot results
png(filename = "./3_Visualization/NMDS_2017_1_2_EM3.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[p2017.data$EM3=="pre-mix 2017",1], Q$points[p2017.data$EM3=="pre-mix 2017",2], pch=21,bg=my.cols[1])
points(Q$points[p2017.data$EM3=="post-EM1 2017",1], Q$points[p2017.data$EM3=="post-EM1 2017",2], pch=21,bg=my.cols[2])
points(Q$points[p2017.data$EM3=="pre-EM2 2017",1], Q$points[p2017.data$EM3=="pre-EM2 2017",2], pch=21,bg=my.cols[3])
points(Q$points[p2017.data$EM3=="post-EM2 2017",1], Q$points[p2017.data$EM3=="post-EM2 2017",2], pch=21,bg=my.cols[4])
points(Q$points[p2017.data$EM3=="post-mix 2017",1], Q$points[p2017.data$EM3=="post-mix 2017",2], pch=21,bg=my.cols[5])
plot(en)
legend("topright",legend = c("pre-mix 2017","post-EM1 2017","pre-EM2 2017",
                             "post-EM2 2017",
                             "post-mix 2017"), pch = 21, pt.bg = my.cols, bty = "n")
legend("bottomright",legend = c("k = 0.08"),bty = "n")
dev.off()

#in black and white
png(filename = "./3_Visualization/NMDS_2017_1_2_EM1.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[p2017.data$EM1==1,1], Q$points[p2017.data$EM1==1,2], pch=21,bg="black")
points(Q$points[p2017.data$EM1==0,1], Q$points[p2017.data$EM1==0,2], pch=21,bg="white")
plot(en)
legend("topright",legend = c("EM + 2 wks post","no EM"), pch = 21, pt.bg = c("black","white"), bty = "n")
legend("bottomright",legend = c("k = 0.13"),bty = "n")
dev.off()

#coded by month
my.cols <- gg_color_hue(5)
png(filename = "./3_Visualization/NMDS_2017_1_2_month.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[p2017.data$Month==5,1], Q$points[p2017.data$Month==5,2], pch=21,bg=my.cols[1])
points(Q$points[p2017.data$Month==6,1], Q$points[p2017.data$Month==6,2], pch=21,bg=my.cols[2])
points(Q$points[p2017.data$Month==7,1], Q$points[p2017.data$Month==7,2], pch=21,bg=my.cols[3])
points(Q$points[p2017.data$Month==8,1], Q$points[p2017.data$Month==8,2], pch=21,bg=my.cols[4])
points(Q$points[p2017.data$Month==9,1], Q$points[p2017.data$Month==9,2], pch=21,bg=my.cols[5])
plot(en)
legend("topright", legend=c('May','June', 'July','Aug','Sept'), pch=21, pt.bg=c(my.cols),bty = "n") 
legend("bottomright",legend = c("k = 0.08"),bty = "n")
dev.off()

#plot w/ genera overlain - this is just for preliminary interpretation, 
#these will be difficult to read
png(filename = "./3_Visualization/NMDS_2017_1_2_genera.png",width = 9,height = 9,units = "in",res = 300)
fig <- ordiplot(Q, choices = c(1,2))
points(Q$points[p2017.data$EM1==1,1], Q$points[p2017.data$EM1==1,2], pch=21,bg='black')
text(fig, "species", col="blue", cex=0.6)
legend("topright",legend = c("EM + 2 wks post","no EM"), pch = 21, pt.bg = c("black","white"), bty = "n")
dev.off()


#### 2018 ####

#subset just 2018 data
p2018.data <- phytos1 %>%
  filter(year(Date) == 2018)

#prep data for NMDS input
p2018 <- p2018.data %>%
  select(-Date,-EM1,-EM2,-EM3,-Month)
p2018 <- as.matrix(p2018)

#select out environmental variables that are constant in 2018
env2018 <- env[36:50,] %>%
  select(-Year, -EM1, -EM2)

##scree plot
#  Basic outline for making your own Scree plot
stress <- rep(NA,6)
k <- c(1:6)
for (i in 1:length(k)){
  Q <- metaMDS(p2018, distance='bray', k=i, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
  stress[i] <- Q$stress
}
plot(k, stress, type='b', xlab='Dimensions', ylab='Stress', pch=21, bg='grey75', cex=1.3)
legend("topright",legend = c("2018"),bty = 'n')

#based on the scree plot the best value for k is 2

#run NMDS for all months
Q <- metaMDS(p2018, distance='bray', k=2, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
plot(Q, display=c('sites'), choices=c(1,2), type='p')

#get initial environmental relationships
en_0 = envfit(Q, env2018, permutations = 999, na.rm = TRUE)
en_0

# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
env_2018 <- env2018 %>%
  select(DOCmax_depth_m)

#run final correlations w/ significant environmental variables
en = envfit(Q, env_2018, permutations = 999, na.rm = TRUE)
en

#no EM plot because there was no EM in 2018

my.cols <- gg_color_hue(5)
#plot results by month
png(filename = "./3_Visualization/NMDS_2018_1_2_month.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[p2018.data$Month==5,1], Q$points[p2018.data$Month==5,2], pch=21,bg=my.cols[1])
points(Q$points[p2018.data$Month==6,1], Q$points[p2018.data$Month==6,2], pch=21,bg=my.cols[2])
points(Q$points[p2018.data$Month==7,1], Q$points[p2018.data$Month==7,2], pch=21,bg=my.cols[3])
points(Q$points[p2018.data$Month==8,1], Q$points[p2018.data$Month==8,2], pch=21,bg=my.cols[4])
points(Q$points[p2018.data$Month==9,1], Q$points[p2018.data$Month==9,2], pch=21,bg=my.cols[5])
plot(en)
legend("topright", legend=c('May','June', 'July','Aug','Sept'), pch=21, pt.bg=c(my.cols),bty = "n") 
legend("bottomright",legend = c("k = 0.12"),bty = "n")
dev.off()

#plot w/ genera overlain - this is just for preliminary interpretation, 
#these will be difficult to read
png(filename = "./3_Visualization/NMDS_2018_1_2_genera.png",width = 9,height = 9,units = "in",res = 300)
fig <- ordiplot(Q, choices = c(1,2))
text(fig, "species", col="blue", cex=0.6)
legend("bottomright",legend = c("k = 0.12"),bty = "n")
dev.off()


#### 2019 ####

#subset 2019 data only
p2019.data <- phytos1 %>%
  filter(year(Date) == 2019)

#prep data for NMDS input
p2019 <- p2019.data %>%
  select(-Date,-EM1,-EM2,-EM3,-Month)
p2019 <- as.matrix(p2019)

#select out environmental variables that are constant in 2019
env2019 <- env[51:67,] %>%
  select(-Year, -EM1, -EM2)

##scree plot
#  Basic outline for making your own Scree plot
stress <- rep(NA,6)
k <- c(1:6)
for (i in 1:length(k)){
  Q <- metaMDS(p2019, distance='bray', k=i, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
  stress[i] <- Q$stress
}
plot(k, stress, type='b', xlab='Dimensions', ylab='Stress', pch=21, bg='grey75', cex=1.3)
legend("topright",legend = c("2019"),bty = 'n')

#based on scree plot the best value for k is 2

#run NMDS for all months
Q <- metaMDS(p2019, distance='bray', k=2, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
plot(Q, display=c('sites'), choices=c(1,2), type='p')

#get initial environmental relationships
en_0 = envfit(Q, env2019, permutations = 999, na.rm = TRUE)
en_0

# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
env_2019 <- env2019 %>%
  select(DOCmax_depth_m, Temp_C, perc_light_thermocline)

#run final correlations w/ significant environmental variables
en = envfit(Q, env_2019, permutations = 999, na.rm = TRUE)
en

#no EM plot because there was no EM in 2019

my.cols <- gg_color_hue(5)
#plot results
png(filename = "./3_Visualization/NMDS_2019_1_2_month.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[p2019.data$Month==5,1], Q$points[p2019.data$Month==5,2], pch=21,bg=my.cols[1])
points(Q$points[p2019.data$Month==6,1], Q$points[p2019.data$Month==6,2], pch=21,bg=my.cols[2])
points(Q$points[p2019.data$Month==7,1], Q$points[p2019.data$Month==7,2], pch=21,bg=my.cols[3])
points(Q$points[p2019.data$Month==8,1], Q$points[p2019.data$Month==8,2], pch=21,bg=my.cols[4])
points(Q$points[p2019.data$Month==9,1], Q$points[p2019.data$Month==9,2], pch=21,bg=my.cols[5])
plot(en)
legend("topright", legend=c('May','June', 'July','Aug','Sept'), pch=21, pt.bg=c(my.cols),bty = "n") 
legend("bottomright",legend = c("k = 0.05"),bty = "n")
dev.off()

#plot w/ genera overlain - this is just for preliminary interpretation, 
#these will be difficult to read
png(filename = "./3_Visualization/NMDS_2019_1_2_genera.png",width = 9,height = 9,units = "in",res = 300)
fig <- ordiplot(Q, choices = c(1,2))
text(fig, "species", col="blue", cex=0.6)
legend("bottomright",legend = c("k = 0.05"),bty = "n")
dev.off()

#### MC test for significance - COME BACK TO THIS ####
MC_stress <- rep(NA,50)
for (j in 1:50){
  df <- randomizeMatrix(hellinger_vars, null.model = c("frequency"))
  Q <- metaMDS(df, distance='bray', k=4, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
  MC_stress[j] <- Q$stress
}
mean(MC_stress)
min(MC_stress)
max(MC_stress)
Q <- readRDS(file = "alltimes_NMDS_results.rds")
Q$stress #[1] 0.1520395
test <- subset(MC_stress, MC_stress <= 0.1520395)
p = (1+length(test))/(1+length(MC_stress))

