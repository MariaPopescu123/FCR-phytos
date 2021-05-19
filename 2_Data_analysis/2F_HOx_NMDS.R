#Title: NMDS for FCR phytos
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
env <- my.cs.data[,c(2:6)]
env$Month <- month(my.cs.data$Date)
env <- data.frame(scale(env, center = TRUE, scale = TRUE)) 

# select out environmental response vars
dvars <- my.cs.data[,c(7:19,21:30)]
dvars <- data.frame(scale(dvars, center = TRUE, scale = TRUE)) 
dvars <- as.matrix(dvars)


####RUN THE NMDS FOR ALL YEARS####

# #calculate Hellinger distance - FOR COMMUNITY DATA
# WHY WOULD I DO THIS??
# hellinger_vars <- decostand(phytos[,-c(1:3)], method = "hellinger")

##scree plot
stress <- rep(NA,6)
k <- c(1:6)
for (i in 1:length(k)){
  Q <- metaMDS(dvars, distance='euclidean', k=i, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE, na.rm = TRUE)
  stress[i] <- Q$stress
}
plot(k, stress, type='b', xlab='Dimensions', ylab='Stress', pch=21, bg='grey75', cex=1.3)
legend("topright",legend = c("all years"),bty = 'n')

#based on results of scree plot, the best choice for k is 4

#run NMDS for all months
set.seed(1)
Q <- metaMDS(dvars, distance='euclidean', k=4, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE, na.rm = TRUE)
Q$stress
Q$points

#run initial correlations with environmental variables
set.seed(1)
en_0 = envfit(Q, env, permutations = 9999, na.rm = TRUE, choices = c(1,2))
en_0

set.seed(1)
en_1 = envfit(Q, env, permutations = 9999, na.rm = TRUE, choices = c(1,3))
en_1

set.seed(1)
en_2 = envfit(Q, env, permutations = 9999, na.rm = TRUE, choices = c(1,4))
en_2

set.seed(1)
en_3 = envfit(Q, env, permutations = 9999, na.rm = TRUE, choices = c(2,3))
en_3

set.seed(1)
en_4 = envfit(Q, env, permutations = 9999, na.rm = TRUE, choices = c(2,4))
en_4

set.seed(1)
en_5 = envfit(Q, env, permutations = 9999, na.rm = TRUE, choices = c(3,4))
en_5

#run ANOSIMs
dist <- vegdist(dvars, method = "euclidean", na.rm = TRUE)

ano_HOx = anosim(dist, env$HOx, permutations = 9999)
ano_HOx
# ANOSIM statistic R: 0.02591 
# Significance: 0.3219 

ano_Month = anosim(dist, env$Month, permutations = 9999)
ano_Month
# ANOSIM statistic R: 0.3714 
# Significance: 1e-04

ano_Year = anosim(dist, env$Year, permutations = 9999)
ano_Year
# ANOSIM statistic R: 0.1214 
# Significance: 3e-04

ano_EM1 = anosim(dist, env$EM1, permutations = 9999)
ano_EM1
# ANOSIM statistic R: -0.0188 
# Significance: 0.581  

ano_EM2 = anosim(dist, env$EM2, permutations = 9999)
ano_EM2
# ANOSIM statistic R: 0.1337 
# Significance: 1e-04 

ano_EM3 = anosim(dist, env$EM3, permutations = 9999)
ano_EM3
# ANOSIM statistic R: 0.1562 
# Significance: 0.0114  

# can't do indicator spp analysis w/ NA values
# but can try to kludge it w/ envfit
# ugh nm it deletes whole rows when things are missing
# so we end up running on only 1/3 of data
dvars.df <- data.frame(dvars)

# plot
par(mfrow = c(2,3))
plot(Q, display=c('sites'), choices=c(1,2), type='p')
plot(en_0)
plot(Q, display=c('sites'), choices=c(1,3), type='p')
plot(en_1)
plot(Q, display=c('sites'), choices=c(1,4), type='n')
plot(en_2)
plot(Q, display=c('sites'), choices=c(2,3), type='n')
plot(en_3)
plot(Q, display=c('sites'), choices=c(2,4), type='n')
plot(en_4)
plot(Q, display=c('sites'), choices=c(3,4), type='n')
plot(en_5)



# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
env_all_12 <- env %>%
  select(Year, HOx, EM2, schmidt.stability, thermo.depth, DO_mgL, pH, Kd, pz_depth_m, Max_biomass_ugL, Peak_depth_m, Month)
env_all_13 <- env %>%
  select(Year, Temp_C, thermo.depth, DO_mgL, pH, Kd, pz_depth_m, Max_biomass_ugL, Peak_depth_m, Month)
env_all <- env %>%
  select(Year, Temp_C, HOx, EM2, schmidt.stability, thermo.depth, DO_mgL, pH, Kd, pz_depth_m, Max_biomass_ugL, Peak_depth_m, Month)

#run final correlations w/ significant environmental variables
set.seed(1)
en12 = envfit(Q, env_all_12, permutations = 9999, na.rm = TRUE, choices = c(1,2))
en12
set.seed(1)
en13 = envfit(Q, env_all_13, permutations = 9999, na.rm = TRUE, choices = c(1,3))
en13

en_all = envfit(Q~EM2, env_all, permutations = 9999, na.rm = TRUE, choices = c(1:3))


#plot results with EM periods and relationships w/ environmental variables
ordirgl(Q, display = "sites", envfit = en_all)
plot(Q, type = "n",choices = c(1,3))
points(Q, display = "sites", cex = 1, pch = 16, col = "red",choices = c(1,3))
text(Q, display = "species", cex = 0.5, col = "blue",choices = c(1,3))
ordisurf(Q, env$Peak_depth_m, add = TRUE,choices = c(1,3))
ordisurf(Q,env$EM2, choices = c(1,2), display = "species")
ordisurf(Q,env$thermo.depth, choices = c(1,3), display = "species")
ordisurf(Q,env$Peak_depth_m, choices = c(1,3), display = "species")


#first and second axes by EM2
png(filename = "./3_Visualization/NMDS_all_1_2_EM2.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[phytos1$EM2==1,1], Q$points[phytos1$EM2==1,2], pch=21,bg="black", cex = 2)
points(Q$points[phytos1$EM2==0,1], Q$points[phytos1$EM2==0,2], pch=21,bg="white", cex = 2)
plot(en12)
legend("topright",legend = c("EM years","years w/ no EM"), pch = 21, pt.bg = c("black","white"), bty = "n", pt.cex = 2)
legend("bottomright",legend = c("k = 0.11"),bty = "n")
dev.off()

#first and second axes by year
my.cols <- gg_color_hue(4)
png(filename = "./3_Visualization/NMDS_all_1_2_year.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[phytos1$Year==2016,1], Q$points[phytos1$Year==2016,2], pch=21,bg=my.cols[1], cex = 2)
points(Q$points[phytos1$Year==2017,1], Q$points[phytos1$Year==2017,2], pch=21,bg=my.cols[2], cex = 2)
points(Q$points[phytos1$Year==2018,1], Q$points[phytos1$Year==2018,2], pch=21,bg=my.cols[3], cex = 2)
points(Q$points[phytos1$Year==2019,1], Q$points[phytos1$Year==2019,2], pch=21,bg=my.cols[4], cex = 2)
plot(en12)
legend("topright",legend = c("2016","2017","2018","2019"), pch = 21, pt.bg = my.cols, bty = "n", pt.cex = 2)
legend("bottomright",legend = c("k = 0.11"),bty = "n")
dev.off()

#first and second axes by peak depth
png(filename = "./3_Visualization/NMDS_all_1_2_peak_depth.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[fp$Peak_depth_indic==1,1], Q$points[fp$Peak_depth_indic==1,2], pch=21,bg=my.cols[1], cex = 2)
points(Q$points[fp$Peak_depth_indic==2,1], Q$points[fp$Peak_depth_indic==2,2], pch=21,bg=my.cols[2], cex = 2)
points(Q$points[fp$Peak_depth_indic==3,1], Q$points[fp$Peak_depth_indic==3,2], pch=21,bg=my.cols[3], cex = 2)
plot(en12)
legend("bottomleft",legend = c("peak depth < 3 m","peak depth = 3-6 m","peak depth > 6 m"), pch = 21, pt.bg = my.cols, bty = "n", pt.cex = 2)
legend("topleft",legend = c("k = 0.11"),bty = "n")
dev.off()

#first and second axes by month
my.cols <- gg_color_hue(5)
png(filename = "./3_Visualization/NMDS_all_1_2_month.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[phytos1$Month==5,1], Q$points[phytos1$Month==5,2], pch=21,bg=my.cols[1], cex = 2)
points(Q$points[phytos1$Month==6,1], Q$points[phytos1$Month==6,2], pch=21,bg=my.cols[2], cex = 2)
points(Q$points[phytos1$Month==7,1], Q$points[phytos1$Month==7,2], pch=21,bg=my.cols[3], cex = 2)
points(Q$points[phytos1$Month==8,1], Q$points[phytos1$Month==8,2], pch=21,bg=my.cols[4], cex = 2)
points(Q$points[phytos1$Month==9,1], Q$points[phytos1$Month==9,2], pch=21,bg=my.cols[5], cex = 2)
plot(en12)
legend("topleft", legend=c('May','June', 'July','Aug','Sept'), pch=21, pt.bg=c(my.cols),bty = "n", pt.cex = 2) 
legend("bottomright",legend = c("k = 0.11"),bty = "n")
dev.off()


#first and third axes by EM2
png(filename = "./3_Visualization/NMDS_all_1_3_EM2.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,3), type='n')
points(Q$points[phytos1$EM2==1,1], Q$points[phytos1$EM2==1,3], pch=21,bg="black", cex = 2)
points(Q$points[phytos1$EM2==0,1], Q$points[phytos1$EM2==0,3], pch=21,bg="white", cex = 2)
plot(en13)
legend("topright",legend = c("EM years","years w/ no EM"), pch = 21, pt.bg = c("black","white"), bty = "n", pt.cex = 2)
legend("bottomright",legend = c("k = 0.11"),bty = "n")
dev.off()

#first and third axes by year
my.cols <- gg_color_hue(4)
png(filename = "./3_Visualization/NMDS_all_1_3_year.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,3), type='n')
points(Q$points[phytos1$Year==2016,1], Q$points[phytos1$Year==2016,3], pch=21,bg=my.cols[1], cex = 2)
points(Q$points[phytos1$Year==2017,1], Q$points[phytos1$Year==2017,3], pch=21,bg=my.cols[2], cex = 2)
points(Q$points[phytos1$Year==2018,1], Q$points[phytos1$Year==2018,3], pch=21,bg=my.cols[3], cex = 2)
points(Q$points[phytos1$Year==2019,1], Q$points[phytos1$Year==2019,3], pch=21,bg=my.cols[4], cex = 2)
plot(en13)
legend("topright",legend = c("2016","2017","2018","2019"), pch = 21, pt.bg = my.cols, bty = "n", pt.cex = 2)
legend("bottomright",legend = c("k = 0.11"),bty = "n")
dev.off()

#first and third axes by peak depth
png(filename = "./3_Visualization/NMDS_all_1_3_peak_depth.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,3), type='n')
points(Q$points[fp$Peak_depth_indic==1,1], Q$points[fp$Peak_depth_indic==1,3], pch=21,bg=my.cols[1], cex = 2)
points(Q$points[fp$Peak_depth_indic==2,1], Q$points[fp$Peak_depth_indic==2,3], pch=21,bg=my.cols[2], cex = 2)
points(Q$points[fp$Peak_depth_indic==3,1], Q$points[fp$Peak_depth_indic==3,3], pch=21,bg=my.cols[3], cex = 2)
plot(en12)
legend("bottomleft",legend = c("peak depth < 3 m","peak depth = 3-6 m","peak depth > 6 m"), pch = 21, pt.bg = my.cols, bty = "n", pt.cex = 2)
legend("topleft",legend = c("k = 0.11"),bty = "n")
dev.off()

#first and third axes by month
png(filename = "./3_Visualization/NMDS_all_1_3_month.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,3), type='n')
points(Q$points[phytos1$Month==5,1], Q$points[phytos1$Month==5,3], pch=21,bg=my.cols[1], cex = 2)
points(Q$points[phytos1$Month==6,1], Q$points[phytos1$Month==6,3], pch=21,bg=my.cols[2], cex = 2)
points(Q$points[phytos1$Month==7,1], Q$points[phytos1$Month==7,3], pch=21,bg=my.cols[3], cex = 2)
points(Q$points[phytos1$Month==8,1], Q$points[phytos1$Month==8,3], pch=21,bg=my.cols[4], cex = 2)
points(Q$points[phytos1$Month==9,1], Q$points[phytos1$Month==9,3], pch=21,bg=my.cols[5], cex = 2)
plot(en13)
legend("topleft", legend=c('May','June', 'July','Aug','Sept'), pch=21, pt.bg=c(my.cols),bty = "n", pt.cex = 2) 
legend("bottomright",legend = c("k = 0.11"),bty = "n")
dev.off()

#plot w/ genera overlain - this is just for preliminary interpretation, 
#these will be difficult to read
png(filename = "./3_Visualization/NMDS_all_1_2_genera.png",width = 9,height = 9,units = "in",res = 300)
fig <- ordiplot(Q, choices = c(1,2))
points(Q$points[phytos1$EM1==1,1], Q$points[phytos1$EM1==1,2], pch=21,bg='black')
text(fig, "species", col="blue", cex=0.6)
plot(en12, col = "red")
legend("topright",legend = c("EM + 2 wks post","no EM"), pch = 21, pt.bg = c("black","white"), bty = "n")
dev.off()

png(filename = "./3_Visualization/NMDS_all_1_3_genera.png",width = 9,height = 9,units = "in",res = 300)
fig <- ordiplot(Q, choices = c(1,3))
points(Q$points[phytos1$EM1==1,1], Q$points[phytos1$EM1==1,3], pch=21,bg='black')
text(fig, "species", col="blue", cex=0.6)
plot(en13, col = "red")
legend("topright",legend = c("EM + 2 wks post","no EM"), pch = 21, pt.bg = c("black","white"), bty = "n")
dev.off()


#### NMDS FOR 2016 ONLY ####

#select 2016 data
p2016.data <- phytos1 %>%
  filter(year(Date) == 2016)

#prep data for NMDS input
p2016 <- p2016.data %>%
  select(-Year,-Date,-EM1,-EM2,-EM3,-Month,-MonthDay)
p2016 <- as.matrix(p2016)

#select out environmental variables that are constant in 2016
env2016 <- env[1:20,] %>%
  select(-Year, -HOx, -EM2)

#select fp data for 2016
fp2016 <- fp %>% filter(year(Date) == 2016)

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
set.seed(1)
Q <- metaMDS(p2016, distance='bray', k=2, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
plot(Q, display=c('sites'), choices=c(1,2), type='p')

#run ANOSIM
ano_EM1 = anosim(p2016, p2016.data$EM1, distance = "bray", permutations = 9999)
ano_EM1
ano_EM3 = anosim(p2016, p2016.data$EM3, distance = "bray", permutations = 9999)
ano_EM3
# ANOSIM statistic R: 0.3022 
# Significance: 0.0086 
ano_Month = anosim(p2016, p2016.data$Month, distance = "bray", permutations = 9999)
ano_Month
# ANOSIM statistic R: 0.496 
# Significance: 1e-04 
ano_FP = anosim(p2016, fp2016$Peak_depth_indic, distance = "bray", permutations = 9999)
ano_FP
#EM3 and Month significant; Month has highest dissimilarity at 0.49

#run month-to-month ANOSIMs
my.combn <- combn(0:3,2)

final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)

for (i in 1:ncol(my.combn)){
  
  #prep data for ANOSIM input
  p2016.temp <- p2016.data %>%
    filter(EM3 %in% my.combn[,i]) 
  p2016 <- p2016.temp %>%
    select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay,-Year)
  p2016 <- as.matrix(p2016)
  
  #anosim
  ano_EM3 = anosim(p2016, p2016.temp$EM3, distance = "bray", permutations = 999)
  
  #store output
  final[i,1] <- my.combn[1,i]
  final[i,2] <- my.combn[2,i]
  final[i,3] <- ano_EM3$statistic
  final[i,4] <- ano_EM3$signif
}
final <- data.frame(final)
colnames(final) <- c("EM_period_A","EM_period_B","R","p")
final.EM3.2016 <- final

#run month-to-month ANOSIMs
my.combn <- combn(5:9,2)

final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)

for (i in 1:ncol(my.combn)){
  
  #prep data for ANOSIM input
  p2016.temp <- p2016.data %>%
    filter(Month %in% my.combn[,i]) 
  p2016 <- p2016.temp %>%
    select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay,-Year)
  p2016 <- as.matrix(p2016)
  
  #anosim
  ano_month = anosim(p2016, p2016.temp$Month, distance = "bray", permutations = 999)
  
  #store output
  final[i,1] <- my.combn[1,i]
  final[i,2] <- my.combn[2,i]
  final[i,3] <- ano_month$statistic
  final[i,4] <- ano_month$signif
}
final <- data.frame(final)
colnames(final) <- c("Month_A","Month_B","R","p")
final.2016 <- final

#run indicator species analysis
inv_EM3 = multipatt(p2016, p2016.data$EM3, func = "r.g", control = how(nperm=9999))
summary(inv_EM3)
# Total number of species: 65
# Selected number of species: 7 
# Number of species associated to 1 group: 5 
# Number of species associated to 2 groups: 1 
# Number of species associated to 3 groups: 1 
# 
# List of species associated to each combination: 
#   
#   Group 0  #sps.  1 
# stat p.value  
# Dinobryon 0.724  0.0116 *
#   
#   Group 1  #sps.  1 
# stat p.value  
# Dolichospermum 0.503  0.0497 *
#   
#   Group 2  #sps.  1 
# stat p.value  
# Aphanocapsa 0.625  0.0292 *
#   
#   Group 3  #sps.  2 
# stat p.value  
# Trachelomonas 0.664  0.0193 *
#   Cyclotella    0.623  0.0242 *
#   
#   Group 0+2  #sps.  1 
# stat p.value  
# Rhodomonas 0.617  0.0371 *
#   
#   Group 0+1+2  #sps.  1 
# stat p.value  
# Ankistrodesmus 0.612  0.0392 *
#   ---

inv_Month = multipatt(p2016, p2016.data$Month, func = "r.g", control = how(nperm=9999))
summary(inv_Month)
# Total number of species: 65
# Selected number of species: 7 
# Number of species associated to 1 group: 3 
# Number of species associated to 2 groups: 4 
# Number of species associated to 3 groups: 0 
# Number of species associated to 4 groups: 0 
# 
# List of species associated to each combination: 
#   
#   Group 8  #sps.  3 
# stat p.value   
# Cryptomonas   0.764  0.0035 **
#   Trachelomonas 0.755  0.0161 * 
#   Woronichinia  0.552  0.0415 * 
#   
#   Group 5+6  #sps.  1 
# stat p.value    
# Ankistrodesmus 0.858   2e-04 ***
#   
#   Group 5+7  #sps.  1 
# stat p.value  
# Rhodomonas 0.63  0.0498 *
#   
#   Group 7+9  #sps.  1 
# stat p.value  
# Pseudanabaena 0.685  0.0214 *
#   
#   Group 8+9  #sps.  1 
# stat p.value  
# Cyclotella 0.627  0.0478 *
#   ---

#get initial environmental relationships
set.seed(1)
en_0 = envfit(Q, env2016, permutations = 9999, na.rm = TRUE)
en_0
p.values <- p.adjust(en_0$vectors$pvals, method = "holm")
p.values

# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
env_2016 <- env2016 %>%
  select(EM3, Temp_C, Kd, Max_biomass_ugL, Month)

#run final correlations w/ significant environmental variables
en = envfit(Q, env_2016, permutations = 999, na.rm = TRUE)
en

#plot results
my.cols <- gg_color_hue(4)
png(filename = "./3_Visualization/NMDS_2016_1_2_EM3.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[p2016.data$EM3==0,1], Q$points[p2016.data$EM3==0,2], pch=21,bg=my.cols[1], cex = 2)
points(Q$points[p2016.data$EM3==1,1], Q$points[p2016.data$EM3==1,2], pch=21,bg=my.cols[2], cex = 2)
points(Q$points[p2016.data$EM3==2,1], Q$points[p2016.data$EM3==2,2], pch=21,bg=my.cols[3], cex = 2)
points(Q$points[p2016.data$EM3==3,1], Q$points[p2016.data$EM3==3,2], pch=21,bg=my.cols[4], cex = 2)
plot(en)
legend("topright",legend = c("pre-mix","post-mix 1","post-mix 2",
                             "post-mix 3"), pch = 21, pt.bg = my.cols, bty = "n", pt.cex = 2)
legend("bottomright",legend = c("k = 0.13"),bty = "n")
dev.off()

#coded by month
my.cols <- gg_color_hue(5)
png(filename = "./3_Visualization/NMDS_2016_1_2_month.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[p2016.data$Month==5,1], Q$points[p2016.data$Month==5,2], pch=21,bg=my.cols[1], cex = 2)
points(Q$points[p2016.data$Month==6,1], Q$points[p2016.data$Month==6,2], pch=21,bg=my.cols[2], cex = 2)
points(Q$points[p2016.data$Month==7,1], Q$points[p2016.data$Month==7,2], pch=21,bg=my.cols[3], cex = 2)
points(Q$points[p2016.data$Month==8,1], Q$points[p2016.data$Month==8,2], pch=21,bg=my.cols[4], cex = 2)
points(Q$points[p2016.data$Month==9,1], Q$points[p2016.data$Month==9,2], pch=21,bg=my.cols[5], cex = 2)
plot(en)
legend("topright", legend=c('May','June', 'July','Aug','Sept'), pch=21, pt.bg=c(my.cols),bty = "n", pt.cex = 2) 
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
  select(-Year,-Date,-EM1,-EM2,-EM3,-Month,-MonthDay)
p2017 <- as.matrix(p2017)

#select out environmental variables that are constant in 2017
env2017 <- env[21:35,] %>%
  select(-Year, -HOx, -EM2, -pH)

#select fp data for 2017
fp2017 <- fp %>% filter(year(Date) == 2017)

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
set.seed(1)
Q <- metaMDS(p2017, distance='bray', k=2, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
plot(Q, display=c('sites'), choices=c(1,2), type='p')

#run ANOSIM
ano_EM1 = anosim(p2017, p2017.data$EM1, distance = "bray", permutations = 9999)
ano_EM1
ano_EM3 = anosim(p2017, p2017.data$EM3, distance = "bray", permutations = 9999)
ano_EM3
ano_Month = anosim(p2017, p2017.data$Month, distance = "bray", permutations = 9999)
ano_Month
# ANOSIM statistic R: 0.4738 
# Significance: 0.0018 
ano_FP = anosim(p2017, fp2017$Peak_depth_indic, distance = "bray", permutations = 9999)
ano_FP
#Month is significant at 0.47 

#run month-to-month ANOSIMs
my.combn <- combn(5:9,2)

final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)

for (i in 1:ncol(my.combn)){
  
  #prep data for ANOSIM input
  p2017.temp <- p2017.data %>%
    filter(Month %in% my.combn[,i]) 
  p2017 <- p2017.temp %>%
    select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay,-Year)
  p2017 <- as.matrix(p2017)
  
  #anosim
  ano_month = anosim(p2017, p2017.temp$Month, distance = "bray", permutations = 999)
  
  #store output
  final[i,1] <- my.combn[1,i]
  final[i,2] <- my.combn[2,i]
  final[i,3] <- ano_month$statistic
  final[i,4] <- ano_month$signif
}
final <- data.frame(final)
colnames(final) <- c("Month_A","Month_B","R","p")
final.2017 <- final

#run indicator species analysis
inv_Month = multipatt(p2017, p2017.data$Month, func = "r.g", control = how(nperm=9999))
summary(inv_Month)
# Total number of species: 65
# Selected number of species: 2 
# Number of species associated to 1 group: 1 
# Number of species associated to 2 groups: 1 
# Number of species associated to 3 groups: 0 
# Number of species associated to 4 groups: 0 
# 
# List of species associated to each combination: 
#   
#   Group 9  #sps.  1 
# stat p.value  
# Synedra 0.984  0.0103 *
#   
#   Group 8+9  #sps.  1 
# stat p.value   
# Cyclotella 0.859   0.005 **
#   ---

#get initial environmental relationships
en_0 = envfit(Q, env2017, permutations = 9999, na.rm = TRUE)
en_0
p.values <- p.adjust(en_0$vectors$pvals, method = "holm")
p.values

# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
env_2017 <- env2017 %>%
  select(Cmax_SRP_ugL, Cmax_DIN_ugL, Cmax_DOC_mgL, DO_mgL, Max_biomass_ugL, Month)

#run final correlations w/ significant environmental variables
en = envfit(Q, env_2017, permutations = 999, na.rm = TRUE)
en

#coded by month
my.cols <- gg_color_hue(5)
png(filename = "./3_Visualization/NMDS_2017_1_2_month.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[p2017.data$Month==5,1], Q$points[p2017.data$Month==5,2], pch=21,bg=my.cols[1], cex = 2)
points(Q$points[p2017.data$Month==6,1], Q$points[p2017.data$Month==6,2], pch=21,bg=my.cols[2], cex = 2)
points(Q$points[p2017.data$Month==7,1], Q$points[p2017.data$Month==7,2], pch=21,bg=my.cols[3], cex = 2)
points(Q$points[p2017.data$Month==8,1], Q$points[p2017.data$Month==8,2], pch=21,bg=my.cols[4], cex = 2)
points(Q$points[p2017.data$Month==9,1], Q$points[p2017.data$Month==9,2], pch=21,bg=my.cols[5], cex = 2)
plot(en)
legend("topright", legend=c('May','June', 'July','Aug','Sept'), pch=21, pt.bg=c(my.cols),bty = "n", pt.cex = 2) 
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
  select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay,-Year)
p2018 <- as.matrix(p2018)

#select out environmental variables that are constant in 2018
env2018 <- env[36:50,] %>%
  select(-Year, -EM1, -EM2, -EM3)

#select fp data for 2018
fp2018 <- fp %>% filter(year(Date) == 2018)

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
set.seed(1)
Q <- metaMDS(p2018, distance='bray', k=2, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
plot(Q, display=c('sites'), choices=c(1,2), type='p')

#run ANOSIM
ano_Month = anosim(p2018, p2018.data$Month, distance = "bray", permutations = 9999)
ano_Month
# ANOSIM statistic R: 0.238 
# Significance: 0.0432 
ano_FP = anosim(p2018, fp2018$Peak_depth_indic, distance = "bray", permutations = 9999)
ano_FP
#Month is significant at R = 0.24

#run month-to-month ANOSIMs
my.combn <- combn(5:9,2)

final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)

for (i in 1:ncol(my.combn)){
  
  #prep data for ANOSIM input
  p2018.temp <- p2018.data %>%
    filter(Month %in% my.combn[,i]) 
  p2018 <- p2018.temp %>%
    select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay,-Year)
  p2018 <- as.matrix(p2018)
  
  #anosim
  ano_month = anosim(p2018, p2018.temp$Month, distance = "bray", permutations = 999)
  
  #store output
  final[i,1] <- my.combn[1,i]
  final[i,2] <- my.combn[2,i]
  final[i,3] <- ano_month$statistic
  final[i,4] <- ano_month$signif
}
final <- data.frame(final)
colnames(final) <- c("Month_A","Month_B","R","p")
final.2018 <- final


#run indicator species analysis
inv_Month = multipatt(p2018, p2018.data$Month, func = "r.g", control = how(nperm=9999))
summary(inv_Month)
# Total number of species: 65
# Selected number of species: 2 
# Number of species associated to 1 group: 2 
# Number of species associated to 2 groups: 0 
# Number of species associated to 3 groups: 0 
# Number of species associated to 4 groups: 0 
# 
# List of species associated to each combination: 
#   
#   Group 5  #sps.  1 
# stat p.value  
# Chlamydomonas 0.688  0.0496 *
#   
#   Group 9  #sps.  1 
# stat p.value  
# Pseudanabaena 0.777  0.0462 *
#   ---

#get initial environmental relationships
set.seed(1)
en_0 = envfit(Q, env2018, permutations = 9999, na.rm = TRUE)
en_0
p.values <- p.adjust(en_0$vectors$pvals, method = "holm")
p.values

# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
env_2018 <- env2018 %>%
  select(HOx, Cmax_DOC_mgL, Temp_C, Kd, pz_depth_m, Peak_depth_m, Month)

#run final correlations w/ significant environmental variables
en = envfit(Q, env_2018, permutations = 9999, na.rm = TRUE)
en

my.cols <- gg_color_hue(5)
#plot results by month
png(filename = "./3_Visualization/NMDS_2018_1_2_month.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[p2018.data$Month==5,1], Q$points[p2018.data$Month==5,2], pch=21,bg=my.cols[1], cex = 2)
points(Q$points[p2018.data$Month==6,1], Q$points[p2018.data$Month==6,2], pch=21,bg=my.cols[2], cex = 2)
points(Q$points[p2018.data$Month==7,1], Q$points[p2018.data$Month==7,2], pch=21,bg=my.cols[3], cex = 2)
points(Q$points[p2018.data$Month==8,1], Q$points[p2018.data$Month==8,2], pch=21,bg=my.cols[4], cex = 2)
points(Q$points[p2018.data$Month==9,1], Q$points[p2018.data$Month==9,2], pch=21,bg=my.cols[5], cex = 2)
plot(en)
legend("topright", legend=c('May','June', 'July','Aug','Sept'), pch=21, pt.bg=c(my.cols),bty = "n", pt.cex = 2) 
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
  select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay,-Year)
p2019 <- as.matrix(p2019)

#select out environmental variables that are constant in 2019
env2019 <- env[51:67,] %>%
  select(-Year, -EM1, -EM2, -EM3, -DO_mgL, -pH, -perc_light_Cmax)

#select fp data for 2019
fp2019 <- fp %>% filter(year(Date) == 2019)

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
set.seed(1)
Q <- metaMDS(p2019, distance='bray', k=2, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
plot(Q, display=c('sites'), choices=c(1,2), type='p')

#run ANOSIM
ano_Month = anosim(p2019, p2019.data$Month, distance = "bray", permutations = 9999)
ano_Month
# ANOSIM statistic R: 0.3285 
# Significance: 0.0115 
ano_FP = anosim(p2019, fp2019$Peak_depth_indic, distance = "bray", permutations = 9999)
ano_FP
# ANOSIM statistic R: 0.2512 
# Significance: 0.0203 
#Month is significant at R = 0.33 and peak depth is significant at R = 0.25

#run month-to-month ANOSIMs
my.combn <- combn(5:9,2)

final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)

for (i in 1:ncol(my.combn)){
  
  #prep data for ANOSIM input
  p2019.temp <- p2019.data %>%
    filter(Month %in% my.combn[,i]) 
  p2019 <- p2019.temp %>%
    select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay,-Year)
  p2019 <- as.matrix(p2019)
  
  #anosim
  ano_month = anosim(p2019, p2019.temp$Month, distance = "bray", permutations = 999)
  
  #store output
  final[i,1] <- my.combn[1,i]
  final[i,2] <- my.combn[2,i]
  final[i,3] <- ano_month$statistic
  final[i,4] <- ano_month$signif
}
final <- data.frame(final)
colnames(final) <- c("Month_A","Month_B","R","p")
final.2019 <- final


#run indicator species analysis
inv_Month = multipatt(p2019, p2019.data$Month, func = "r.g", control = how(nperm=9999))
summary(inv_Month)
# Total number of species: 65
# Selected number of species: 5 
# Number of species associated to 1 group: 4 
# Number of species associated to 2 groups: 1 
# Number of species associated to 3 groups: 0 
# Number of species associated to 4 groups: 0 
# 
# List of species associated to each combination: 
#   
#   Group 7  #sps.  1 
# stat p.value   
# Spondylosium 0.853  0.0074 **
#   
#   Group 9  #sps.  3 
# stat p.value   
# Microcystis 0.902  0.0085 **
#   Euglena     0.729  0.0085 **
#   Cryptomonas 0.695  0.0325 * 
#   
#   Group 8+9  #sps.  1 
# stat p.value   
# Trachelomonas 0.839  0.0096 **
#   ---
  
#get initial environmental relationships
set.seed(1)
en_0 = envfit(Q, env2019, permutations = 9999, na.rm = TRUE)
en_0
p.values <- p.adjust(en_0$vectors$pvals, method = "holm")
p.values

# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
env_2019 <- env2019 %>%
  select(HOx, Temp_C, thermo.depth, Kd, pz_depth_m, Peak_depth_m, Month)

#run final correlations w/ significant environmental variables
#run final correlations w/ significant environmental variables
en = envfit(Q, env_2019, permutations = 9999, na.rm = TRUE)
en

#no EM plot because there was no EM in 2019

#plot results
my.cols <- gg_color_hue(5)
png(filename = "./3_Visualization/NMDS_2019_1_2_month.png",width = 9,height = 9,units = "in",res = 300)
plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
points(Q$points[p2019.data$Month==5,1], Q$points[p2019.data$Month==5,2], pch=21,bg=my.cols[1], cex = 2)
points(Q$points[p2019.data$Month==6,1], Q$points[p2019.data$Month==6,2], pch=21,bg=my.cols[2], cex = 2)
points(Q$points[p2019.data$Month==7,1], Q$points[p2019.data$Month==7,2], pch=21,bg=my.cols[3], cex = 2)
points(Q$points[p2019.data$Month==8,1], Q$points[p2019.data$Month==8,2], pch=21,bg=my.cols[4], cex = 2)
points(Q$points[p2019.data$Month==9,1], Q$points[p2019.data$Month==9,2], pch=21,bg=my.cols[5], cex = 2)
plot(en)
legend("topright", legend=c('May','June', 'July','Aug','Sept'), pch=21, pt.bg=c(my.cols),bty = "n", pt.cex = 2) 
legend("bottomright",legend = c("k = 0.05"),bty = "n")
dev.off()

#plot w/ genera overlain - this is just for preliminary interpretation, 
#these will be difficult to read
png(filename = "./3_Visualization/NMDS_2019_1_2_genera.png",width = 9,height = 9,units = "in",res = 300)
fig <- ordiplot(Q, choices = c(1,2))
text(fig, "species", col="blue", cex=0.6)
plot(en)
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

