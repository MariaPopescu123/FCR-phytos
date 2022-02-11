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
env <- my.cs.data[,c(2,5:6,11:13,27,29:32,39,40:41,43)]
env$Month <- month(my.cs.data$Date)
env <- data.frame(scale(env, center = TRUE, scale = TRUE)) 

#read in fp data
my.fp.data <- read_csv("./2_Data_analysis/FP_megamatrix.csv")%>%
  mutate(Peak_depth_indic = ifelse(Peak_depth_m < 3.1,1,
                                   ifelse(Peak_depth_m >= 3.1 & Peak_depth_m <3.6,1,
                                          ifelse(Peak_depth_m >= 3.6 & Peak_depth_m <4.9,2,
                                                 ifelse(Peak_depth_m >= 4.9 ,2,NA)))))

# hist(my.fp.data$Peak_depth_m)
# mean(my.fp.data$Peak_depth_m, na.rm = TRUE)
# sd(my.fp.data$Peak_depth_m, na.rm = TRUE)
# quantile(my.fp.data$Peak_depth_m, probs = c(0.33, 0.5, 0.66), na.rm = TRUE)


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
  select(Date, Peak_depth_m, Peak_depth_indic)

# prepare phyto data for NMDS input
phytos2 <- phytos1 %>%
  select(-Date,-EM1,-EM2,-EM3,-Month,-Year,-MonthDay) 

phytos3 <- as.matrix(phytos2)

####RUN THE NMDS FOR ALL YEARS####

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
plot(Q, display=c('sites'), choices=c(1,2), type='p')
plot(Q, display=c('sites'), choices=c(1,3), type='p')
# Q$species
# Q$points

#run ANOSIMs on all year-month pairings
phytos1 <- phytos1 %>%
  mutate(MonthYear = as.numeric(factor(format(Date, "%y-%m"))))
my.combn.may <- combn(c(1,6,11,16),2)
my.combn.june <- combn(c(2,7,12,17),2)
my.combn.july <- combn(c(3,8,13,18),2)
my.combn.aug <- combn(c(4,9,14,19),2)
my.combn.sept <- combn(c(5,10,15,20),2)
my.combn <- cbind(my.combn.may,my.combn.june,my.combn.july,my.combn.aug,my.combn.sept)

final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)

for (i in 1:ncol(my.combn)){
  
  #prep data for ANOSIM input
  phytos1.temp <- phytos1 %>%
    filter(MonthYear %in% my.combn[,i]) 
  p1 <- phytos1.temp %>%
    select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay,-Year,-MonthYear)
  p1 <- as.matrix(p1)
  
  #anosim
  ano_month = anosim(p1, phytos1.temp$MonthYear, distance = "bray", permutations = 999)
  
  #store output
  final[i,1] <- my.combn[1,i]
  final[i,2] <- my.combn[2,i]
  final[i,3] <- ano_month$statistic
  final[i,4] <- ano_month$signif
}
final <- data.frame(final)
colnames(final) <- c("MonthYear_A","MonthYear_B","R","p")
final.m.y <- final
final.m.y$diff = ifelse(final.m.y$p < 0.05, "yes","no")

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
# ANOSIM statistic R: 0.281 
# Significance: 1e-04 

#EM1 and EM3 not significant; dissimilarity low even across significant groups;
#highest dissimilarity Year and Month

###comparing Years pairwise
#select 2016-2017 data
p2016_17.data <- phytos1 %>%
  filter(year(Date) %in% c(2016:2017))

#prep data for NMDS input
p2016_17 <- p2016_17.data %>%
  select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay)
p2016_17 <- as.matrix(p2016_17)

#run ANOSIM
ano_16_17 = anosim(p2016_17, p2016_17.data$Year, distance = "bray", permutations = 9999)
ano_16_17
# ANOSIM statistic R: 0.0393 
# Significance: 0.1718

#select 2016-2018 data
p2016_18.data <- phytos1 %>%
  filter(year(Date) %in% c(2016,2018))

#prep data for NMDS input
p2016_18 <- p2016_18.data %>%
  select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay)
p2016_18 <- as.matrix(p2016_18)

#run ANOSIM
ano_16_18 = anosim(p2016_18, p2016_18.data$Year, distance = "bray", permutations = 9999)
ano_16_18
# ANOSIM statistic R: 0.2313 
# Significance: 5e-04 

#select 2016-2019 data
p2016_19.data <- phytos1 %>%
  filter(year(Date) %in% c(2016,2019))

#prep data for NMDS input
p2016_19 <- p2016_19.data %>%
  select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay)
p2016_19 <- as.matrix(p2016_19)

#run ANOSIM
ano_16_19 = anosim(p2016_19, p2016_19.data$Year, distance = "bray", permutations = 9999)
ano_16_19
# ANOSIM statistic R: 0.3138 
# Significance: 1e-04 

#select 2017-2019 data
p2017_19.data <- phytos1 %>%
  filter(year(Date) %in% c(2017,2019))

#prep data for NMDS input
p2017_19 <- p2017_19.data %>%
  select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay)
p2017_19 <- as.matrix(p2017_19)

#run ANOSIM
ano_17_19 = anosim(p2017_19, p2017_19.data$Year, distance = "bray", permutations = 9999)
ano_17_19
# ANOSIM statistic R: 0.1858 
# Significance: 0.0112 

#select 2017-2018 data
p2017_18.data <- phytos1 %>%
  filter(year(Date) %in% c(2017,2018))

#prep data for NMDS input
p2017_18 <- p2017_18.data %>%
  select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay)
p2017_18 <- as.matrix(p2017_18)

#run ANOSIM
ano_17_18 = anosim(p2017_18, p2017_18.data$Year, distance = "bray", permutations = 9999)
ano_17_18
# ANOSIM statistic R: 0.1062 
# Significance: 0.0269

#select 2018-2019 data
p2018_19.data <- phytos1 %>%
  filter(year(Date) %in% c(2018,2019))

#prep data for NMDS input
p2018_19 <- p2018_19.data %>%
  select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay)
p2018_19 <- as.matrix(p2018_19)

#run ANOSIM
ano_18_19 = anosim(p2018_19, p2018_19.data$Year, distance = "bray", permutations = 9999)
ano_18_19
# ANOSIM statistic R: 0.2703 
# Significance: 0.0025 

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

inv_PD = multipatt(phytos3, fp$Peak_depth_indic, func = "r.g", control = how(nperm=9999))
summary(inv_PD)
# Multilevel pattern analysis
# 
#   
#   Association function: r.g
# Significance level (alpha): 0.05
# 
# Total number of species: 65
# Selected number of species: 10 
# Number of species associated to 1 group: 10 
# 
# List of species associated to each combination: 
#   
#   Group 1  #sps.  8 
#                     stat p.value    
#   Dolichospermum    0.405  0.0001 ***
#   Chlorophyte sp. 1 0.291  0.0019 ** 
#   Rhodomonas        0.276  0.0226 *  
#   Nitzchia          0.255  0.0056 ** 
#   Parvodinium       0.238  0.0373 *  
#   Spondylosium      0.225  0.0383 *  
#   Euglena           0.213  0.0471 *  
#   Selenastrum       0.149  0.0065 ** 
#   
#   Group 2  #sps.  2 
#                stat p.value  
#   Cryptomonas  0.243  0.0254 *
#   Elakatothrix 0.188  0.0122 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

#run initial correlations with environmental variables
set.seed(1)
en_0 = envfit(Q, env, permutations = 9999, na.rm = TRUE, choices = c(1,2))
en_0
# ***VECTORS
# 
# NMDS1    NMDS2     r2 Pr(>r)    
# Year              -0.53913 -0.84222 0.3222 0.0002 ***
# EM2                0.62498  0.78064 0.3325 0.0002 ***
# EM3                0.94967  0.31326 0.4016 0.0001 ***
# Grab_SRP_ugL       0.03239 -0.99948 0.0194 0.6159    
# Grab_DIN_ugL       0.89055 -0.45488 0.1079 0.0597 .  
# Grab_DOC_mgL      -0.18763 -0.98224 0.1595 0.0150 *  
# Temp_C_grab        0.54894 -0.83586 0.1834 0.0069 ** 
# schmidt.stability -0.39434 -0.91896 0.1410 0.0248 *  
# n2                 0.11560 -0.99330 0.1791 0.0072 ** 
# perc_light_grab   -0.74054 -0.67201 0.1106 0.0535 .  
# Max_biomass_ugL    0.74929 -0.66225 0.3431 0.0002 ***
# Peak_depth_m       0.68579  0.72780 0.2057 0.0032 ** 
# Peak_width_m       0.06416  0.99794 0.0720 0.1584    
# Month              0.83414 -0.55155 0.4523 0.0001 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 9999
# 
# 15 observations deleted due to missingness

##RESULTS INCLUDING LAKE NUMBER AND WEDDERBURN NUMBER

# ***VECTORS
# 
# NMDS1    NMDS2     r2 Pr(>r)    
# Year              -0.53913 -0.84222 0.3222 0.0002 ***
# EM2                0.62498  0.78064 0.3325 0.0002 ***
# EM3                0.94967  0.31326 0.4016 0.0001 ***
# Grab_SRP_ugL       0.03239 -0.99948 0.0194 0.6159    
# Grab_DIN_ugL       0.89055 -0.45488 0.1079 0.0597 .  
# Grab_DOC_mgL      -0.18763 -0.98224 0.1595 0.0150 *  
# Temp_C_grab        0.54894 -0.83586 0.1834 0.0069 ** 
# schmidt.stability -0.33322 -0.94285 0.1548 0.0153 *  
# n2                 0.11560 -0.99330 0.1791 0.0072 ** 
# lake.number        0.36277 -0.93188 0.3402 0.0001 ***
# wedderburn.number  0.95423 -0.29908 0.2849 0.0002 ***
# perc_light_grab   -0.74054 -0.67201 0.1106 0.0535 .  
# Max_biomass_ugL    0.74929 -0.66225 0.3431 0.0002 ***
# Peak_depth_m       0.68579  0.72780 0.2057 0.0032 ** 
# Peak_width_m       0.06416  0.99794 0.0720 0.1584    
# Month              0.83414 -0.55155 0.4523 0.0001 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 9999
# 
# 15 observations deleted due to missingness

set.seed(1)
en_1 = envfit(Q, env, permutations = 9999, na.rm = TRUE, choices = c(1,3))
en_1

# ***VECTORS
# 
# NMDS1    NMDS3     r2 Pr(>r)    
# Year              -0.22602 -0.97412 0.1581 0.0154 *  
# EM2                0.76190  0.64770 0.0915 0.0961 .  
# EM3                0.97155  0.23685 0.3608 0.0002 ***
# Grab_SRP_ugL       0.46447 -0.88559 0.0241 0.5545    
# Grab_DIN_ugL       0.93439  0.35625 0.0950 0.0858 .  
# Grab_DOC_mgL       0.36125 -0.93247 0.0549 0.2559    
# Temp_C_grab        0.66549 -0.74641 0.2838 0.0004 ***
# schmidt.stability  0.07271 -0.99735 0.0728 0.1647    
# n2                 0.99836 -0.05720 0.0159 0.6808    
# perc_light_grab   -0.23606 -0.97174 0.1800 0.0095 ** 
# Max_biomass_ugL    0.91520 -0.40300 0.2559 0.0006 ***
# Peak_depth_m       0.16264  0.98669 0.3226 0.0001 ***
# Peak_width_m      -0.38755 -0.92185 0.0053 0.8815    
# Month              0.98170 -0.19043 0.3287 0.0001 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 9999
# 
# 15 observations deleted due to missingness

##RESULTS INCLUDING LAKE NUMBER AND WEDDERBURN NUMBER

# ***VECTORS
# 
# NMDS1    NMDS3     r2 Pr(>r)    
# Year              -0.22602 -0.97412 0.1581 0.0154 *  
# EM2                0.76190  0.64770 0.0915 0.0961 .  
# EM3                0.97155  0.23685 0.3608 0.0002 ***
# Grab_SRP_ugL       0.46447 -0.88559 0.0241 0.5545    
# Grab_DIN_ugL       0.93439  0.35625 0.0950 0.0858 .  
# Grab_DOC_mgL       0.36125 -0.93247 0.0549 0.2559    
# Temp_C_grab        0.66549 -0.74641 0.2838 0.0004 ***
# schmidt.stability  0.15598 -0.98776 0.0682 0.1860    
# n2                 0.99836 -0.05720 0.0159 0.6808    
# lake.number        0.87427 -0.48545 0.1157 0.0497 *  
# wedderburn.number  0.94803  0.31818 0.2812 0.0006 ***
# perc_light_grab   -0.23606 -0.97174 0.1800 0.0095 ** 
# Max_biomass_ugL    0.91520 -0.40300 0.2559 0.0006 ***
# Peak_depth_m       0.16264  0.98669 0.3226 0.0001 ***
# Peak_width_m      -0.38755 -0.92185 0.0053 0.8815    
# Month              0.98170 -0.19043 0.3287 0.0001 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 9999
# 
# 15 observations deleted due to missingness


# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
env_all_12 <- env %>%
  select(Year, EM2, EM3, lake.number, wedderburn.number, Max_biomass_ugL, Month)
env_all_13 <- env %>%
  select(EM3, Temp_C_grab, wedderburn.number, Max_biomass_ugL, Peak_depth_m, Month)
env_all <- env %>%
  select(Year, EM2, EM3, Temp_C_grab, Max_biomass_ugL, Peak_depth_m, Month)

#run final correlations w/ significant environmental variables
set.seed(1)
en12 = envfit(Q, env_all_12, permutations = 9999, na.rm = TRUE, choices = c(1,2))
en12
check <- env_all_12 %>%
  filter(complete.cases(.))
set.seed(1)
en13 = envfit(Q, env_all_13, permutations = 9999, na.rm = TRUE, choices = c(1,3))
en13

en_all = envfit(Q, env_all, permutations = 9999, na.rm = TRUE, choices = c(1:2))


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
  select(-Year, -EM2)

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
# ANOSIM statistic R: 0.02853 
# Significance: 0.3237
#EM3 and Month significant; Month has highest dissimilarity at 0.49

#run mixing period ANOSIMs
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

inv_PD = multipatt(p2016, fp2016$Peak_depth_indic, func = "r.g", control = how(nperm=9999))
summary(inv_PD)
# Multilevel pattern analysis
# 
#   
#   Association function: r.g
# Significance level (alpha): 0.05
# 
# Total number of species: 65
# Selected number of species: 2 
# Number of species associated to 1 group: 2 
# 
# List of species associated to each combination: 
#   
#   Group 1  #sps.  2 
# stat p.value    
# Chlorophyte sp. 1 0.791  0.0006 ***
#   Aphanocapsa       0.484  0.0178 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

#get initial environmental relationships
set.seed(1)
en_0 = envfit(Q, env2016, permutations = 9999, na.rm = TRUE)
en_0
# ***VECTORS
# 
# NMDS1    NMDS2     r2 Pr(>r)    
# EM3                0.22855 -0.97353 0.7202 0.0002 ***
# Grab_SRP_ugL       0.59232 -0.80570 0.0954 0.4920    
# Grab_DIN_ugL       0.81772 -0.57562 0.1048 0.4867    
# Grab_DOC_mgL      -0.04159 -0.99913 0.1412 0.3561    
# Temp_C_grab       -0.09686 -0.99530 0.5530 0.0038 ** 
# schmidt.stability -0.37250 -0.92803 0.2261 0.1587    
# n2                 0.13502 -0.99084 0.5813 0.0035 ** 
# perc_light_grab   -0.96004 -0.27986 0.1312 0.3775    
# Max_biomass_ugL    0.63913 -0.76910 0.6576 0.0006 ***
# Peak_depth_m       0.61060  0.79194 0.4217 0.0239 *  
# Peak_width_m      -0.98983 -0.14223 0.2325 0.1623    
# Month              0.12325 -0.99238 0.5974 0.0024 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 9999
# 
# 3 observations deleted due to missingness

##RESULTS INCLUDING LAKE NUMBER AND WEDDERBURN NUMBER
# ***VECTORS
# 
# NMDS1    NMDS2     r2 Pr(>r)    
# EM3                0.22855 -0.97353 0.7202 0.0002 ***
# Grab_SRP_ugL       0.59232 -0.80570 0.0954 0.4920    
# Grab_DIN_ugL       0.81772 -0.57562 0.1048 0.4867    
# Grab_DOC_mgL      -0.04159 -0.99913 0.1412 0.3561    
# Temp_C_grab       -0.09686 -0.99530 0.5530 0.0038 ** 
# schmidt.stability -0.37250 -0.92803 0.2261 0.1587    
# n2                 0.13502 -0.99084 0.5813 0.0035 ** 
# lake.number        0.53514 -0.84476 0.4099 0.0255 *  
# wedderburn.number  0.99699  0.07749 0.5560 0.0031 ** 
# perc_light_grab   -0.96004 -0.27986 0.1312 0.3775    
# Max_biomass_ugL    0.63913 -0.76910 0.6576 0.0006 ***
# Peak_depth_m       0.61060  0.79194 0.4217 0.0239 *  
# Peak_width_m      -0.98983 -0.14223 0.2325 0.1623    
# Month              0.12325 -0.99238 0.5974 0.0024 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 9999
# 
# 3 observations deleted due to missingness

# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
env_2016 <- env2016 %>%
  select(EM3, Temp_C, Kd, Max_biomass_ugL, Month)

#run final correlations w/ significant environmental variables
en = envfit(Q, env_2016, permutations = 999, na.rm = TRUE)
en


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
  select(-Year, -EM2)

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
# ANOSIM statistic R: -0.04584 
# Significance: 0.5818
ano_Month = anosim(p2017, p2017.data$Month, distance = "bray", permutations = 9999)
ano_Month
# ANOSIM statistic R: 0.4738 
# Significance: 0.0018 
ano_FP = anosim(p2017, fp2017$Peak_depth_indic, distance = "bray", permutations = 9999)
ano_FP
# ANOSIM statistic R: 0.4611 
# Significance: 0.01
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

inv_PD = multipatt(p2017, fp2017$Peak_depth_indic, func = "r.g", control = how(nperm=9999))
summary(inv_PD)
# Association function: r.g
# Significance level (alpha): 0.05
# 
# Total number of species: 65
# Selected number of species: 5 
# Number of species associated to 1 group: 5 
# 
# List of species associated to each combination: 
#   
#   Group 1  #sps.  3 
# stat p.value    
# Nitzchia          0.672   3e-04 ***
#   Dolichospermum    0.495   4e-03 ** 
#   Chlorophyte sp. 1 0.487   5e-04 ***
#   
#   Group 2  #sps.  2 
# stat p.value  
# Cyclotella      0.592  0.0329 *
#   Dictyosphaerium 0.470  0.0465 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

#get initial environmental relationships
set.seed(1)
en_0 = envfit(Q, env2017, permutations = 9999, na.rm = TRUE)
en_0

# ***VECTORS
# 
# NMDS1    NMDS2     r2 Pr(>r)    
# EM3               -0.08278 -0.99657 0.4206 0.0605 .  
# Grab_SRP_ugL       0.94454 -0.32839 0.7486 0.0015 ** 
# Grab_DIN_ugL      -0.44983 -0.89311 0.4266 0.0562 .  
# Grab_DOC_mgL       0.96241 -0.27161 0.8032 0.0005 ***
# Temp_C_grab        0.76203 -0.64754 0.4986 0.0293 *  
# schmidt.stability  0.86627  0.49957 0.1780 0.3754    
# n2                -0.22652 -0.97401 0.1539 0.4409    
# perc_light_grab    0.99744 -0.07150 0.3772 0.0858 .  
# Max_biomass_ugL   -0.01728 -0.99985 0.5384 0.0244 *  
# Peak_depth_m      -0.99775 -0.06711 0.3303 0.1336    
# Peak_width_m       0.18602  0.98255 0.0774 0.6693    
# Month             -0.19940 -0.97992 0.6001 0.0079 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 9999
# 
# 2 observations deleted due to missingness

##RESULTS INCLUDING LAKE NUMBER AND WEDDERBURN NUMBER
# ***VECTORS
# 
# NMDS1    NMDS2     r2 Pr(>r)    
# EM3               -0.08278 -0.99657 0.4206 0.0605 .  
# Grab_SRP_ugL       0.94454 -0.32839 0.7486 0.0015 ** 
# Grab_DIN_ugL      -0.44983 -0.89311 0.4266 0.0562 .  
# Grab_DOC_mgL       0.96241 -0.27161 0.8032 0.0005 ***
# Temp_C_grab        0.76203 -0.64754 0.4986 0.0293 *  
# schmidt.stability  0.97966  0.20065 0.1253 0.5168    
# n2                -0.22652 -0.97401 0.1539 0.4409    
# lake.number       -0.00031 -1.00000 0.1973 0.3241    
# wedderburn.number -0.58681 -0.80973 0.4915 0.0334 *  
# perc_light_grab    0.99744 -0.07150 0.3772 0.0858 .  
# Max_biomass_ugL   -0.01728 -0.99985 0.5384 0.0244 *  
# Peak_depth_m      -0.99775 -0.06711 0.3303 0.1336    
# Peak_width_m       0.18602  0.98255 0.0774 0.6693    
# Month             -0.19940 -0.97992 0.6001 0.0079 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 9999
# 
# 2 observations deleted due to missingness


# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
env_2017 <- env2017 %>%
  select(Cmax_SRP_ugL, Cmax_DIN_ugL, Cmax_DOC_mgL, DO_mgL, Max_biomass_ugL, Month)

#run final correlations w/ significant environmental variables
en = envfit(Q, env_2017, permutations = 999, na.rm = TRUE)
en



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
  select(-Year,  -EM2, -EM3)

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
# ANOSIM statistic R: -0.016 
# Significance: 0.4661 
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

inv_PD = multipatt(p2018, fp2018$Peak_depth_indic, func = "r.g", control = how(nperm=9999))
summary(inv_PD)
# Association function: r.g
# Significance level (alpha): 0.05
# 
# Total number of species: 65
# Selected number of species: 1 
# Number of species associated to 1 group: 1 
# 
# List of species associated to each combination: 
#   
#   Group 1  #sps.  1 
# stat p.value  
# Cyclotella 0.535   0.022 *

#get initial environmental relationships
set.seed(1)
en_0 = envfit(Q, env2018, permutations = 9999, na.rm = TRUE)
en_0
# ***VECTORS
# 
# NMDS1    NMDS2     r2 Pr(>r)   
# Grab_SRP_ugL      -0.83084 -0.55652 0.4210 0.0808 . 
# Grab_DIN_ugL      -0.92336 -0.38393 0.1958 0.3734   
# Grab_DOC_mgL      -0.95544  0.29518 0.6785 0.0064 **
# Temp_C_grab       -0.98938  0.14533 0.7681 0.0023 **
# schmidt.stability -0.96841  0.24936 0.2778 0.2304   
# n2                -0.51400  0.85779 0.1676 0.4505   
# perc_light_grab   -0.99623 -0.08678 0.5672 0.0238 * 
# Max_biomass_ugL   -0.17364  0.98481 0.4750 0.0570 . 
# Peak_depth_m       0.99254  0.12189 0.6239 0.0156 * 
# Peak_width_m       0.94085  0.33883 0.3092 0.1948   
# Month             -0.37116  0.92857 0.7318 0.0027 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 9999
# 
# 3 observations deleted due to missingness

##RESULTS INCLUDING LAKE NUMBER AND WEDDERBURN NUMBER
# ***VECTORS
# 
# NMDS1    NMDS2     r2 Pr(>r)   
# Grab_SRP_ugL      -0.83084 -0.55652 0.4210 0.0808 . 
# Grab_DIN_ugL      -0.92336 -0.38393 0.1958 0.3734   
# Grab_DOC_mgL      -0.95544  0.29518 0.6785 0.0064 **
# Temp_C_grab       -0.98938  0.14533 0.7681 0.0023 **
# schmidt.stability -0.96841  0.24936 0.2778 0.2304   
# n2                -0.51400  0.85779 0.1676 0.4505   
# lake.number       -0.33209  0.94325 0.4723 0.0588 . 
# wedderburn.number -0.15658  0.98766 0.2583 0.2403   
# perc_light_grab   -0.99623 -0.08678 0.5672 0.0238 * 
# Max_biomass_ugL   -0.17364  0.98481 0.4750 0.0570 . 
# Peak_depth_m       0.99254  0.12189 0.6239 0.0156 * 
# Peak_width_m       0.94085  0.33883 0.3092 0.1948   
# Month             -0.37116  0.92857 0.7318 0.0027 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 9999
# 
# 3 observations deleted due to missingness


# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
env_2018 <- env2018 %>%
  select(HOx, Cmax_DOC_mgL, Temp_C, Kd, pz_depth_m, Peak_depth_m, Month)

#run final correlations w/ significant environmental variables
en = envfit(Q, env_2018, permutations = 9999, na.rm = TRUE)
en


#### 2019 ####

#subset 2019 data only
p2019.data <- phytos1 %>%
  filter(year(Date) == 2019) %>%
  mutate(Storm = ifelse(Date < "2019-06-08",1,2)) 

#prep data for NMDS input
p2019 <- p2019.data %>%
  select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay,-Year,-Storm)
p2019 <- as.matrix(p2019)

#select out environmental variables that are constant in 2019
env2019 <- env[51:67,] %>%
  select(-Year,  -EM2, -EM3, -perc_light_grab)

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
# ANOSIM statistic R: 0.7793 
# Significance: 0.0006 
ano_Storm = anosim(p2019, p2019.data$Storm, distance = "bray", permutations = 9999)
ano_Storm
# ANOSIM statistic R: 0.9838 
# Significance: 3e-04 

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

inv_PD = multipatt(p2019, fp2019$Peak_depth_indic, func = "r.g", control = how(nperm=9999))
summary(inv_PD)
# Total number of species: 65
# Selected number of species: 2 
# Number of species associated to 1 group: 2 
# 
# List of species associated to each combination: 
#   
#   Group 1  #sps.  1 
# stat p.value   
# Dolichospermum 0.658  0.0079 **
#   
#   Group 2  #sps.  1 
# stat p.value  
# Asterionella 0.543  0.0269 *

inv_Storm = multipatt(p2019, p2019.data$Storm, func = "r.g", control = how(nperm=9999))
summary(inv_Storm)
# Association function: r.g
# Significance level (alpha): 0.05
# 
# Total number of species: 65
# Selected number of species: 3 
# Number of species associated to 1 group: 3 
# 
# List of species associated to each combination: 
#   
#   Group 2  #sps.  3 
#   stat p.value   
#   Dolichospermum 0.697  0.0022 **
#   Rhodomonas     0.513  0.0281 * 
#   Cyclotella     0.464  0.0357 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
  
#get initial environmental relationships
set.seed(1)
en_0 = envfit(Q, env2019, permutations = 9999, na.rm = TRUE)
en_0
# ***VECTORS
# 
# NMDS1    NMDS2     r2 Pr(>r)   
# Grab_SRP_ugL      -0.77705  0.62944 0.3931 0.0663 . 
# Grab_DIN_ugL       0.19360 -0.98108 0.0209 0.9024   
# Grab_DOC_mgL      -0.48923  0.87216 0.3148 0.1271   
# Temp_C_grab       -0.16502  0.98629 0.7123 0.0011 **
# schmidt.stability  0.48799  0.87285 0.2324 0.2315   
# n2                 0.71135  0.70284 0.0884 0.5976   
# Max_biomass_ugL    0.13446  0.99092 0.3163 0.1158   
# Peak_depth_m       0.39175 -0.92007 0.5647 0.0086 **
# Peak_width_m       0.92366  0.38322 0.0250 0.8653   
# Month             -0.90255  0.43059 0.5465 0.0134 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 9999
# 
# 3 observations deleted due to missingness

##RESULTS INCLUDING LAKE NUMBER AND WEDDERBURN NUMBER
# ***VECTORS
# 
# NMDS1    NMDS2     r2 Pr(>r)   
# Grab_SRP_ugL      -0.77705  0.62944 0.3931 0.0663 . 
# Grab_DIN_ugL       0.19360 -0.98108 0.0209 0.9024   
# Grab_DOC_mgL      -0.48923  0.87216 0.3148 0.1271   
# Temp_C_grab       -0.16502  0.98629 0.7123 0.0011 **
# schmidt.stability  0.48799  0.87285 0.2324 0.2315   
# n2                 0.71135  0.70284 0.0884 0.5976   
# lake.number       -0.55691  0.83057 0.4478 0.0363 * 
# wedderburn.number -0.42855  0.90352 0.4180 0.0546 . 
# Max_biomass_ugL    0.13446  0.99092 0.3163 0.1158   
# Peak_depth_m       0.39175 -0.92007 0.5647 0.0086 **
# Peak_width_m       0.92366  0.38322 0.0250 0.8653   
# Month             -0.90255  0.43059 0.5465 0.0134 * 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Permutation: free
# Number of permutations: 9999
# 
# 3 observations deleted due to missingness


# selecting environmental drivers that were significant in NMDS for
# all years so plot is less busy
env_2019 <- env2019 %>%
  select(HOx, Temp_C, thermo.depth, Kd, pz_depth_m, Peak_depth_m, Month)

#run final correlations w/ significant environmental variables
#run final correlations w/ significant environmental variables
en = envfit(Q, env_2019, permutations = 9999, na.rm = TRUE)
en

