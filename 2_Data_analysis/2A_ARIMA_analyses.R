#2A_ARIMA_analyses
#Author: Mary Lofton
#Date: 13OCT20

#need to check normality of all vars and correlation of drivers and fix as needed
pacman::p_load(PerformanceAnalytics, tidyverse, lubridate, forecast, MuMIn)
rm(list=ls())


#get data
mydata <- read_csv("./2_Data_analysis/FP_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_Depth_m)
mydata <- mydata[,c(2,1,3:21)]

#look at correlations - cutoff of 0.5
driver_cor <- cor(mydata[,2:17],method = "spearman",use = "complete.obs")
#thermo.depth and Year
#SRPmax_ugL and Cmax_SRP_ugL
#DINmax_ugL and DOC_max_mgL
#DOCmax_mgL and Cmax_DOC_mgL
#schmidt.stability and thermo.depth

#look at correlations of collinear drivers w/ responses to select drivers
mydata1 <- mydata %>%
  select(Peak_width_m,Max_biomass_ugL,Peak_magnitude_ugL,Peak_depth_m,
         thermo.depth,Year,SRPmax_ugL,Cmax_SRP_ugL,DINmax_ugL,DOCmax_mgL,
         Cmax_DOC_mgL,schmidt.stability)
response_cor <- cor(mydata1,method = "spearman",use = "complete.obs")
#get rid of Year
#get rid of SRPmax_ugL
#get rid of DOCmax_mgL
#get rid of schmidt.stability

mydata2 <- mydata %>%
  select(-Year,-SRPmax_ugL,-DOCmax_mgL,-schmidt.stability)

#check for skewness in drivers and responses and whether logging improves it
#ideally we want skewness to approach 0

for (i in 2:17){
  print(colnames(mydata2)[i])
  print(skewness(mydata2[,i], na.rm = TRUE))
  print(skewness(log(mydata2[,i]), na.rm = TRUE))
}
#should log:
#Cmax_DIN_ugL
#Cmax_DOC_mgL
#DINmax_ugL
#DOCmax_depth_m
#Max_biomass_ugL
#Peak_magnitude_ugL

mydata3 <- mydata2 %>%
  mutate(Cmax_DIN_ugL = log(Cmax_DIN_ugL),
         Cmax_DOC_mgL = log(Cmax_DOC_mgL),
         DINmax_ugL = log(DINmax_ugL),
         DOCmax_depth_m = log(DOCmax_depth_m),
         Max_biomass_ugL = log(Max_biomass_ugL),
         Peak_magnitude_ugL = log(Peak_magnitude_ugL))

#scale predictor variables to allow comparison of coefficients in model
mydata4 <- mydata3 %>%
  mutate_at(vars(-Date,-Max_biomass_ugL,-Peak_width_m,-Peak_depth_m,-Peak_magnitude_ugL),scale)

#plot histograms to be sure what's going on here
for (i in 2:17){
  var <- mydata4[,i]
  hist(as.matrix(var),main = colnames(mydata4)[i])
}

#check what order of ARIMA model is correct for response vars

#Peak_depth_m
fit1 <- auto.arima(mydata4$Peak_depth_m)#364.15
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,2]))#365.7
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,3]))#366.19
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,4]))#328.5
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,5]))#367.28
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,6]))#363.71
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,7]))#327.45
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,8]))#363.32
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,9]))#361.29
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,10]))#334.6
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,11]))#267.99
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,12]))#328.01
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,13]))#289.58
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,c(11,13)]))#210.81
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,c(7,11,13)]))#191.5
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,c(7,11,12,13)]))#180.17
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,c(4,7,11,12,13)]))#166.7
#best fit model?? Temp_C,Kd,SRPmax_depth_m,thermo.depth,Cmax_SRP_ugL,DOCmax_depth_m
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,c(4,7,10,11,12,13)]))#154.97
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,c(4,7,9,10,11,12,13)]))#155.89
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,c(4,7,8,9,10,11,12,13)]))#158.2
fit1 <- auto.arima(mydata4$Peak_depth_m,xreg = as.matrix(mydata4[,2:13]))#164.15
fit1

#Peak_width_m
fit2 <- auto.arima(mydata4$Peak_width_m)#396.2
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,2]))#398.32
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,3]))#396.02
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,4]))#376.99
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,5]))#377.75
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,6]))#382.94
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,7]))#321.9
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,8]))#374.97
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,9]))#374.6
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,10]))#358.9
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,11]))#372.74
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,12]))#351.36
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,13]))#322.32
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,c(7,13)]))#266.47
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,c(7,12,13)]))#255.8
#best fit model?? SRPmax_depth_m, Kd, thermo.depth, DOCmax_depth_m
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,c(7,10,12,13)]))#227.36
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,c(7,10,11,12,13)]))#228.1
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,c(7,9,10,11,12,13)]))#225.98
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,c(7,8,9,10,11,12,13)]))#226.81
fit2 <- auto.arima(mydata4$Peak_width_m,xreg = as.matrix(mydata4[,2:13]))#228.84
fit2

#Peak_magnitude_ugL
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL)#249.96
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,2]))#251.11
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,3]))#252.09
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,4]))#240.52
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,5]))#241.17
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,6]))#238.88
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,7]))#213.62
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,8]))#238.66
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,9]))#227.21
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,10]))#227.41
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,11]))#231.96
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,12]))#222.34
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,13]))#204.66
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,c(7,13)]))#176.3
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,c(7,12,13)]))#166.76
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,c(7,9,12,13)]))#160.05
#best fit model?? Kd,SRPmax_depth_m,thermo.depth,DINmax_ugL,DOCmax_depth_m
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,c(7,9,10,12,13)]))#150.97
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,c(7,9,10,11,12,13)]))#153.44
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,c(7,8,9,10,11,12,13)]))#153.27
fit3 <- auto.arima(mydata4$Peak_magnitude_ugL,xreg = as.matrix(mydata4[,2:13]))#156.97
