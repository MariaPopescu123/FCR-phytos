#Title: NMDS for FCR phytos
#Author: Mary Lofton
#Date: 11DEC2017
#Updated: 01FEB2021

####SET-UP####
#load packages
pacman::p_load(tidyverse, lubridate, vegan, Hmisc, picante, indicspecies)
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

#read in fp data
my.fp.data <- read_csv("./2_Data_analysis/FP_megamatrix.csv")%>%
  mutate(Peak_depth_indic = ifelse(Peak_depth_m < 1,1,
                                   ifelse(Peak_depth_m >= 1 & Peak_depth_m <2,1,
                                          ifelse(Peak_depth_m >= 2 & Peak_depth_m <3,1,
                                                 ifelse(Peak_depth_m >= 3 & Peak_depth_m <4,2,
                                                        ifelse(Peak_depth_m >= 4 & Peak_depth_m <5,2,
                                                               ifelse(Peak_depth_m >= 5 & Peak_depth_m <6,2,
                                                                      ifelse(Peak_depth_m >= 6 & Peak_depth_m <7,3,
                                                                             ifelse(Peak_depth_m >= 7 & Peak_depth_m <8,3,
                                                                                    ifelse(Peak_depth_m >= 8 & Peak_depth_m <9,3,
                                                                                           ifelse(Peak_depth_m >= 9 & Peak_depth_m <10,3,NA)))))))))))


#read in community data
phytos <- read_csv("./00_Data_files/EDI_phytos/phytoplankton.csv") %>%
  select(Date, Genus, BV_um3mL) %>%
  spread(key = Genus, value = BV_um3mL)

# any non-detects filled with 0
phytos[is.na(phytos)]<- 0

# create data frame with different EM periods for plotting later on
em <- my.cs.data[,c(1,2,4,5,6)] %>%
  mutate(Month = month(Date))

phytos1 <- left_join(em, phytos, by = "Date")
fp <- left_join(phytos1,my.fp.data, by = "Date") %>%
  select(Date, Peak_depth_m, Peak_depth_indic)

# prepare phyto data for NMDS input
phytos2 <- phytos1 %>%
  select(-Date,-EM1,-EM2,-EM3,-Month,-Year)

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
set.seed(1)
Q <- metaMDS(phytos3, distance='bray', k=3, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
plot(Q, display=c('sites'), choices=c(1,2), type='p')
plot(Q, display=c('sites'), choices=c(1,3), type='p')

#run ANOSIM
ano_EM1 = anosim(phytos3, phytos1$EM1, distance = "bray", permutations = 9999)
ano_EM1
# ANOSIM statistic R: 0.09354 
# Significance: 0.0718 
ano_EM2 = anosim(phytos3, phytos1$EM2, distance = "bray", permutations = 9999)
ano_EM2
# ANOSIM statistic R: 0.1464 
# Significance: 2e-04 
ano_EM3 = anosim(phytos3, phytos1$EM3, distance = "bray", permutations = 9999)
ano_EM3
# ANOSIM statistic R: 0.0642 
# Significance: 0.1154 
ano_Month = anosim(phytos3, phytos1$Month, distance = "bray", permutations = 9999)
ano_Month
# ANOSIM statistic R: 0.1869 
# Significance: 1e-04 
ano_Year = anosim(phytos3, phytos1$Year, distance = "bray", permutations = 9999)
ano_Year
# ANOSIM statistic R: 0.1951 
# Significance: 1e-04 
ano_FP = anosim(phytos3, fp$Peak_depth_indic, distance = "bray", permutations = 9999)
ano_FP
# ANOSIM statistic R: 0.2408 
# Significance: 5e-04 

#EM1 and EM3 not significant; dissimilarity low even across significant groups;
#highest dissimilarity Year and Month

#run indicator species analysis
inv_EM2 = multipatt(phytos3, phytos1$EM2, func = "r.g", control = how(nperm=9999))
summary(inv_EM2)
# Group 0  #sps.  7 
# stat p.value    
# Oocystis     0.357  0.0005 ***
#   Rhodomonas   0.311  0.0074 ** 
#   Staurastrum  0.304  0.0083 ** 
#   Spondylosium 0.227  0.0361 *  
#   Monomastix   0.208  0.0225 *  
#   Staurodesmus 0.194  0.0212 *  
#   Selenastrum  0.146  0.0265 *
inv_Month = multipatt(phytos3, phytos1$Month, func = "r.g", control = how(nperm=9999))
summary(inv_Month)
# Total number of species: 65
# Selected number of species: 12 
# Number of species associated to 1 group: 10 
# Number of species associated to 2 groups: 1 
# Number of species associated to 3 groups: 1 
# Number of species associated to 4 groups: 0 
# 
# List of species associated to each combination: 
#   
#   Group 5  #sps.  1 
# stat p.value  
# Asterionella 0.412  0.0121 *
#   
#   Group 6  #sps.  1 
# stat p.value  
# Snowella 0.34  0.0468 *
#   
#   Group 7  #sps.  3 
# stat p.value  
# Nitzchia       0.392  0.0201 *
#   Spondylosium   0.390  0.0248 *
#   Dolichospermum 0.377  0.0244 *
#   
#   Group 8  #sps.  3 
# stat p.value   
# Cryptomonas   0.453  0.0040 **
#   Trachelomonas 0.390  0.0053 **
#   Woronichinia  0.361  0.0205 * 
#   
#   Group 9  #sps.  2 
# stat p.value    
# Microcystis   0.514  0.0008 ***
#   Pseudanabaena 0.421  0.0091 ** 
#   
#   Group 8+9  #sps.  1 
# stat p.value  
# Prorocentrum 0.371  0.0236 *
#   
#   Group 7+8+9  #sps.  1 
# stat p.value   
# Cyclotella 0.44  0.0064 **
#   ---
inv_Year = multipatt(phytos3, phytos1$Year, func = "r.g", control = how(nperm=9999))
summary(inv_Year)
# Total number of species: 65
# Selected number of species: 9 
# Number of species associated to 1 group: 8 
# Number of species associated to 2 groups: 1 
# Number of species associated to 3 groups: 0 
# 
# List of species associated to each combination: 
#   
#   Group 2017  #sps.  3 
# stat p.value  
# Chlorophyte sp. 4 0.393  0.0190 *
#   Nitzchia          0.352  0.0171 *
#   Micractinium      0.226  0.0185 *
#   
#   Group 2018  #sps.  4 
# stat p.value    
# Staurastrum  0.584  0.0001 ***
#   Oocystis     0.582  0.0001 ***
#   Staurodesmus 0.341  0.0049 ** 
#   Selenastrum  0.258  0.0078 ** 
#   
#   Group 2019  #sps.  1 
# stat p.value   
# Spondylosium 0.389  0.0047 **
#   
#   Group 2017+2018  #sps.  1 
# stat p.value  
# Parvodinium 0.338  0.0272 *
#   ---
inv_FP = multipatt(phytos3, fp$Peak_depth_indic, func = "r.g", control = how(nperm=9999))
summary(inv_FP)
# Total number of species: 65
# Selected number of species: 3 
# Number of species associated to 1 group: 3 
# Number of species associated to 2 groups: 0 
# 
# List of species associated to each combination: 
#   
#   Group 1  #sps.  2 
# stat p.value  
# Dolichospermum    0.469  0.0184 *
#   Chlorophyte sp. 1 0.380  0.0481 *
#   
#   Group 3  #sps.  1 
# stat p.value  
# Cryptomonas 0.411   0.027 *
#   ---


#run initial correlations with environmental variables
set.seed(1)
en_0 = envfit(Q, env, permutations = 999, na.rm = TRUE, choices = c(1,2))
en_0
set.seed(1)
en_1 = envfit(Q, env, permutations = 999, na.rm = TRUE, choices = c(1,3))
en_1

# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
env_all_12 <- env %>%
  select(Year, EM2, DINmax_ugL, DOCmax_depth_m,
         DOCmax_mgL, Kd)
env_all_13 <- env %>%
  select(Year, Temp_C,  Kd)

#run final correlations w/ significant environmental variables
set.seed(1)
en12 = envfit(Q, env_all_12, permutations = 999, na.rm = TRUE, choices = c(1,2))
en12
set.seed(1)
en13 = envfit(Q, env_all_13, permutations = 999, na.rm = TRUE, choices = c(1,3))
en13

#get colors for plotting
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#plot results with EM periods and relationships w/ environmental variables

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
  select(-Date,-EM1,-EM2,-EM3,-Month)
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
# ANOSIM statistic R: 0.2992 
# Significance: 0.011 
ano_Month = anosim(p2016, p2016.data$Month, distance = "bray", permutations = 9999)
ano_Month
ano_FP = anosim(p2016, fp2016$Peak_depth_indic, distance = "bray", permutations = 9999)
ano_FP
#EM3 and Month significant; Month has highest dissimilarity at 0.49

#run indicator species analysis
inv_EM3 = multipatt(p2016, p2016.data$EM3, func = "r.g", control = how(nperm=9999))
summary(inv_EM3)
# Total number of species: 66
# Selected number of species: 7 
# Number of species associated to 1 group: 5 
# Number of species associated to 2 groups: 1 
# Number of species associated to 3 groups: 1 
# 
# List of species associated to each combination: 
#   
#   Group 0  #sps.  1 
# stat p.value   
# Dinobryon 0.724  0.0099 **
#   
#   Group 1  #sps.  1 
# stat p.value  
# Dolichospermum 0.503  0.0471 *
#   
#   Group 2  #sps.  1 
# stat p.value  
# Aphanocapsa 0.625   0.029 *
#   
#   Group 3  #sps.  2 
# stat p.value  
# Trachelomonas 0.664  0.0198 *
#   Cyclotella    0.623  0.0243 *
#   
#   Group 0+2  #sps.  1 
# stat p.value  
# Rhodomonas 0.617  0.0331 *
#   
#   Group 0+1+2  #sps.  1 
# stat p.value  
# Ankistrodesmus 0.612  0.0395 *
#   ---
inv_Month = multipatt(p2016, p2016.data$Month, func = "r.g", control = how(nperm=9999))
summary(inv_Month)
# Total number of species: 66
# Selected number of species: 6 
# Number of species associated to 1 group: 3 
# Number of species associated to 2 groups: 3 
# Number of species associated to 3 groups: 0 
# Number of species associated to 4 groups: 0 
# 
# List of species associated to each combination: 
#   
#   Group 8  #sps.  3 
# stat p.value   
# Cryptomonas   0.764  0.0039 **
#   Trachelomonas 0.755  0.0148 * 
#   Woronichinia  0.552  0.0373 * 
#   
#   Group 5+6  #sps.  1 
# stat p.value    
# Ankistrodesmus 0.858   3e-04 ***
#   
#   Group 7+9  #sps.  1 
# stat p.value  
# Pseudanabaena 0.685  0.0187 *
#   
#   Group 8+9  #sps.  1 
# stat p.value  
# Cyclotella 0.627  0.0475 *
#   ---

#get initial environmental relationships
set.seed(1)
en_0 = envfit(Q, env2016, permutations = 999, na.rm = TRUE)
en_0

# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
env_2016 <- env2016 %>%
  select(Kd)

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
  select(-Date,-EM1,-EM2,-EM3,-Month)
p2017 <- as.matrix(p2017)

#select out environmental variables that are constant in 2017
env2017 <- env[21:35,] %>%
  select(-Year, -HOx, -EM2)

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
ano_FP = anosim(p2017, fp2017$Peak_depth_indic, distance = "bray", permutations = 9999)
ano_FP
#Month is significant at 0.47 

#run indicator species analysis
inv_Month = multipatt(p2017, p2017.data$Month, func = "r.g", control = how(nperm=9999))
summary(inv_Month)
# Total number of species: 66
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
# Synedra 0.984  0.0113 *
#   
#   Group 8+9  #sps.  1 
# stat p.value   
# Cyclotella 0.859  0.0046 **
#   ---

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
  select(-Date,-EM1,-EM2,-EM3,-Month)
p2018 <- as.matrix(p2018)

#select out environmental variables that are constant in 2018
env2018 <- env[36:50,] %>%
  select(-Year, -EM1, -EM2)

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
Q <- metaMDS(p2018, distance='bray', k=2, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
plot(Q, display=c('sites'), choices=c(1,2), type='p')

#run ANOSIM
ano_Month = anosim(p2018, p2018.data$Month, distance = "bray", permutations = 9999)
ano_Month
ano_FP = anosim(p2018, fp2018$Peak_depth_indic, distance = "bray", permutations = 9999)
ano_FP
#Month is significant at R = 0.24

#run indicator species analysis
inv_Month = multipatt(p2018, p2018.data$Month, func = "r.g", control = how(nperm=9999))
summary(inv_Month)
# Total number of species: 66
# Selected number of species: 3 
# Number of species associated to 1 group: 2 
# Number of species associated to 2 groups: 1 
# Number of species associated to 3 groups: 0 
# Number of species associated to 4 groups: 0 
# 
# List of species associated to each combination: 
#   
#   Group 5  #sps.  1 
# stat p.value  
# Chlamydomonas 0.688  0.0495 *
#   
#   Group 9  #sps.  1 
# stat p.value  
# Pseudanabaena 0.777  0.0472 *
#   
#   Group 7+9  #sps.  1 
# stat p.value  
# Cyclotella 0.719  0.0497 *
#   ---

#get initial environmental relationships
en_0 = envfit(Q, env2018, permutations = 999, na.rm = TRUE)
en_0

# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
# NONE!

#run final correlations w/ significant environmental variables
# NONE!

#no EM plot because there was no EM in 2018

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
  select(-Date,-EM1,-EM2,-EM3,-Month)
p2019 <- as.matrix(p2019)

#select out environmental variables that are constant in 2019
env2019 <- env[51:67,] %>%
  select(-Year, -EM1, -EM2)

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
Q <- metaMDS(p2019, distance='bray', k=2, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
plot(Q, display=c('sites'), choices=c(1,2), type='p')

#run ANOSIM
ano_Month = anosim(p2019, p2019.data$Month, distance = "bray", permutations = 9999)
ano_Month
ano_FP = anosim(p2019, fp2019$Peak_depth_indic, distance = "bray", permutations = 9999)
ano_FP
# ANOSIM statistic R: 0.2512 
# Significance: 0.0203 
#Month is significant at R = 0.33 and peak depth is significant at R = 0.35

#run indicator species analysis
inv_Month = multipatt(p2019, p2019.data$Month, func = "r.g", control = how(nperm=9999))
summary(inv_Month)
inv_FP = multipatt(p2019, fp2019$Peak_depth_indic, func = "r.g", control = how(nperm=9999))
summary(inv_FP)
# Total number of species: 66
# Selected number of species: 2 
# Number of species associated to 1 group: 2 
# Number of species associated to 2 groups: 0 
# 
# List of species associated to each combination: 
#   
#   Group 1  #sps.  1 
# stat p.value  
# Parvodinium 0.775  0.0298 *
#   
#   Group 3  #sps.  1 
# stat p.value  
# Ankistrodesmus 0.882  0.0214 *
#   ---

#get initial environmental relationships
en_0 = envfit(Q, env2019, permutations = 999, na.rm = TRUE)
en_0

# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
env_2019 <- env2019 %>%
  select(DOCmax_depth_m, Temp_C)

#run final correlations w/ significant environmental variables
en = envfit(Q, env_2019, permutations = 999, na.rm = TRUE)
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

