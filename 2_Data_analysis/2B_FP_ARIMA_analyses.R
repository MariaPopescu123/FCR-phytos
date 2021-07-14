#2A_FP_ARIMA_analyses
#Author: Mary Lofton
#Date: 13OCT20

####SET-UP####

#need to check normality of all vars and correlation of drivers and fix as needed
pacman::p_load(PerformanceAnalytics, tidyverse, lubridate, forecast, utils, igraph,RColorBrewer)
rm(list=ls())


#get data
my.fp.data <- read_csv("./2_Data_analysis/FP_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_DO_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20" | MonthDay == "01-01") %>%
  select(-MonthDay)
my.fp.data <- my.fp.data[,c(2,1,3:47)]

mydata <- my.fp.data
colnames(mydata)
#peak depth drivers = c(14,16,18,31,37)
#peak width drivers = c(23,24,25,29,30,34)
#max biomass drivers = c(7,8,9,28,29,30,36)

#look at correlations - cutoff of 0.7

#for peak depth
driver_cor <- cor(mydata[,c(40,14,16,18,29,30,31,34)],method = "spearman",use = "pairwise.complete.obs")
driver_cor[lower.tri(driver_cor)] = ""
write.csv(driver_cor, file = "./2_Data_analysis/PD_driver_cor.csv", row.names = FALSE)
#schmidt.stability and n2

mydata1 <- mydata %>%
  select(Peak_depth_m,schmidt.stability,n2,pz_depth_m,DINmax_depth_m)
response_cor <- cor(mydata1,method = "spearman",use = "pairwise.complete.obs")
response_cor[lower.tri(response_cor)] = ""
#use schmidt.stability 


#for each year
d <- mydata %>% filter(Year %in% c(2016:2017))
driver_cor <- cor(d[,c(40,14,16,18,29,30,31,34)],method = "spearman",use = "pairwise.complete.obs")
driver_cor[lower.tri(driver_cor)] = ""
write.csv(driver_cor, file = "./2_Data_analysis/PD_driver_cor_man.csv", row.names = FALSE)
#Kd

mydata1 <- mydata %>%
  filter(Year %in% c(2016:2017)) %>%
  select(Peak_depth_m,Kd,pz_depth_m,DINmax_depth_m)
response_cor <- cor(mydata1,method = "spearman",use = "pairwise.complete.obs")
response_cor[lower.tri(response_cor)] = ""

d <- mydata %>% filter(Year %in% c(2018:2019))
driver_cor <- cor(d[,c(40,14,16,18,29,30,31,34)],method = "spearman",use = "pairwise.complete.obs")
driver_cor[lower.tri(driver_cor)] = ""
write.csv(driver_cor, file = "./2_Data_analysis/PD_driver_ref.csv", row.names = FALSE)
#schmidt.stability and n2


#look at correlations of collinear drivers w/ responses to select drivers
mydata1 <- mydata %>%
  filter(Year %in% c(2018:2019)) %>%
  select(Peak_depth_m,schmidt.stability,n2,pz_depth_m,thermo.depth,Kd)
response_cor <- cor(mydata1,method = "spearman",use = "pairwise.complete.obs")
response_cor[lower.tri(response_cor)] = ""
#use schmidt.stability
#Kd and thermo.depth better than pz_depth_m






####for peak width
driver_cor <- cor(mydata[,c(42,23,24,25,29,30,34)],method = "spearman",use = "pairwise.complete.obs")
driver_cor[lower.tri(driver_cor)] = ""
write.csv(driver_cor, file = "./2_Data_analysis/PW_driver_cor.csv", row.names = FALSE)
#schmidt.stability and n2

#look at correlations of collinear drivers w/ responses to select drivers
mydata1 <- mydata %>%
  select(Peak_width_m,schmidt.stability,n2)
response_cor <- cor(mydata1,method = "spearman",use = "pairwise.complete.obs")
response_cor[lower.tri(response_cor)] = ""
#n2 better 

d <- mydata %>% filter(Year %in% c(2016:2017))
driver_cor <- cor(d[,c(42,23,24,25,29,30,34)],method = "spearman",use = "pairwise.complete.obs")
driver_cor[lower.tri(driver_cor)] = ""
write.csv(driver_cor, file = "./2_Data_analysis/PW_driver_cor_man.csv", row.names = FALSE)


d <- mydata %>% filter(Year %in% c(2018:2019))
driver_cor <- cor(d[,c(42,23,24,25,29,30,34)],method = "spearman",use = "pairwise.complete.obs")
driver_cor[lower.tri(driver_cor)] = ""
write.csv(driver_cor, file = "./2_Data_analysis/PW_driver_cor_ref.csv", row.names = FALSE)


#look at correlations of collinear drivers w/ responses to select drivers
mydata1 <- mydata %>%
  filter(Year %in% c(2018:2019)) %>%
  select(Peak_width_m,schmidt.stability,n2)
response_cor <- cor(mydata1,method = "spearman",use = "pairwise.complete.obs")
response_cor[lower.tri(response_cor)] = ""
#schmidt.stability better 

#going w/ n2 b/c better for whole time series




####for max biomass
driver_cor <- cor(mydata[,c(39,7,8,9,28,29,30,36)],method = "spearman",use = "pairwise.complete.obs")
driver_cor[lower.tri(driver_cor)] = ""
write.csv(driver_cor, file = "./2_Data_analysis/MB_driver_cor.csv", row.names = FALSE)
#schmidt.stability and n2

mydata1 <- mydata %>%
  select(Max_biomass_ugL,schmidt.stability,n2)
response_cor <- cor(mydata1,method = "spearman",use = "pairwise.complete.obs")
response_cor[lower.tri(response_cor)] = ""
#use n2


d <- mydata %>% filter(Year %in% c(2016:2017))
driver_cor <- cor(d[,c(39,7,8,9,28,29,30,36)],method = "spearman",use = "pairwise.complete.obs")
driver_cor[lower.tri(driver_cor)] = ""
write.csv(driver_cor, file = "./2_Data_analysis/MB_driver_cor_man.csv", row.names = FALSE)

d <- mydata %>% filter(Year %in% c(2018:2019))
driver_cor <- cor(d[,c(39,7,8,9,28,29,30,36)],method = "spearman",use = "pairwise.complete.obs")
driver_cor[lower.tri(driver_cor)] = ""
write.csv(driver_cor, file = "./2_Data_analysis/MB_driver_cor_ref.csv", row.names = FALSE)

#Temp_C_Cmax and Cmax_DOC_mgL
#schmidt.stability and n2

#look at correlations of collinear drivers w/ responses to select drivers
mydata1 <- mydata %>%
  filter(Year %in% c(2018:2019)) %>%
  select(Max_biomass_ugL,Temp_C_Cmax,Cmax_DOC_mgL,schmidt.stability,n2)
response_cor <- cor(mydata1,method = "spearman",use = "pairwise.complete.obs")
response_cor[lower.tri(response_cor)] = ""
#use Temp_C_Cmax
#schmidt.stability and n2 about the same so going w/ n2 based on whole timeseries




colnames(mydata)
#subset to drivers and responses
#peak depth: get rid of n2
#peak width: get rid of schmidt.stability
#max biomass: get rid of Cmax_DOC_mgL and schmidt.stability
#peak depth drivers = c(14,16,18,29,30,31,34)
#peak width drivers = c(23,24,25,29,30,34)
#max biomass drivers = c(7,8,28,29,30,36)
mydata2 <- mydata[,c(1,40,42,39,14,16,18,29,30,31,34,23,24,25,7,8,28,36)]

#check for skewness in drivers and responses and whether logging improves it
#ideally we want skewness to approach 0

for (i in 2:18){
  print(colnames(mydata2)[i])
  var <- mydata2[,i]
  hist(as.matrix(var), main = colnames(mydata2)[i])
  print(skewness(mydata2[,i], na.rm = TRUE))
  print(skewness(log(mydata2[,i]+0.0001), na.rm = TRUE))
  var <- log(mydata2[,i])
  hist(as.matrix(var), main = c("Log",colnames(mydata2)[i]))
}
#should log:

#Max_biomass_ugL 
#SRPmax_depth_m 
#DOCmax_depth_m 
#thermo.depth
#pz_depth_m
#pz_DIN_CV
#Cmax_DIN_ugL
#perc_light_Cmax



mydata3 <- mydata2 %>%
  mutate(Cmax_DIN_ugL = log(Cmax_DIN_ugL),
         SRPmax_depth_m = log(SRPmax_depth_m),
         DOCmax_depth_m = log(DOCmax_depth_m),
         Max_biomass_ugL = log(Max_biomass_ugL),
         pz_DIN_CV = log(pz_DIN_CV),
         perc_light_Cmax = log(perc_light_Cmax),
         thermo.depth = log(thermo.depth)
         )

#scale predictor/response variables to allow comparison of coefficients in model
mydata4 <- mydata3 %>%
  mutate_at(vars(-Date),scale) 

# ##check skewness for EM years
d <- mydata[,c(1,2,40,42,39,14,16,18,29,30,31,34,23,24,25,7,8,28,36)] %>%
  filter(Year %in% c(2016:2017)) %>%
  select(-Year)
for (i in 2:18){
  print(colnames(d)[i])
  var <- d[,i]
  hist(as.matrix(var), main = colnames(d)[i])
  print(skewness(d[,i], na.rm = TRUE))
  print(skewness(log(d[,i]+0.0001), na.rm = TRUE))
  var <- log(d[,i])
  hist(as.matrix(var), main = c("Log",colnames(d)[i]))
}
#should log:
#Max_biomass_ugL
#DOCmax_depth_m
#pz_depth_m
#pz_DIN_CV
#pz_DOC_CV
#Cmax_DIN_ugL
#perc_light_Cmax


d1 <- d %>%
  mutate(Cmax_DIN_ugL = log(Cmax_DIN_ugL),
         Max_biomass_ugL = log(Max_biomass_ugL),
         pz_DIN_CV = log(pz_DIN_CV),
         pz_DOC_CV = log(pz_DOC_CV),
         DOCmax_depth_m = log(DOCmax_depth_m),
         perc_light_Cmax = log(perc_light_Cmax)
  )

#scale predictor/response variables to allow comparison of coefficients in model
dEM <- d1 %>%
  mutate_at(vars(-Date),scale) 


##non-EM years
d <- mydata[,c(1,2,40,42,39,14,16,18,29,30,31,34,23,24,25,7,8,28,36)] %>%
  filter(Year %in% c(2018:2019)) %>%
  select(-Year)
for (i in 2:18){
  print(colnames(d)[i])
  var <- d[,i]
  hist(as.matrix(var), main = colnames(d)[i])
  print(skewness(d[,i], na.rm = TRUE))
  print(skewness(log(d[,i]+0.0001), na.rm = TRUE))
  var <- log(d[,i])
  hist(as.matrix(var), main = c("Log",colnames(d)[i]))
}
#should log:
#Max_biomass_ugL
#DINmax_depth_m
#DOCmax_depth_m
#thermo.depth
#Cmax_DIN_ugL
#perc_light_Cmax

d1 <- d %>%
  mutate(Max_biomass_ugL = log(Max_biomass_ugL),
         DOCmax_depth_m = log(DOCmax_depth_m),
         DINmax_depth_m = log(DINmax_depth_m),
         thermo.depth = log(thermo.depth),
         Cmax_DIN_ugL = log(Cmax_DIN_ugL),
         perc_light_Cmax = log(perc_light_Cmax))

#scale predictor/response variables to allow comparison of coefficients in model
dnoEM <- d1 %>%
  mutate_at(vars(-Date),scale) 

####ARIMAS FOR ALL YEARS####

#find best-fit ARIMA model for peak depth
#subset to drivers and responses
#peak depth drivers = c(10,14,31)
#peak width drivers = c(19,20,24,28)
#max biomass drivers = c(8,22)
colnames(mydata4)


cols <- c(5:8,10:11)
sub.final <- NULL
final <- NULL

y <- mydata4[,2]

  for (i in 1:length(cols)){
    my.combn <- combn(cols,i)
    sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
    for (j in 1:ncol(my.combn)){
      
      skip_to_next <- FALSE
      
      tryCatch(fit <- auto.arima(y,xreg = as.matrix(mydata4[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
      
      if(skip_to_next) { 
        sub.sub.final[j,4] <- NA
        sub.sub.final[j,3] <- j
        sub.sub.final[j,2] <- i
        sub.sub.final[j,1] <- "Peak_depth_m"
        next }

      sub.sub.final[j,4] <- fit$aicc
      sub.sub.final[j,3] <- j
      sub.sub.final[j,2] <- i
      sub.sub.final[j,1] <- "Peak_depth_m"
    }
    
    sub.final <- rbind(sub.final,sub.sub.final)
    print(paste("I have finished with all combinations of length",i,"for Peak_depth_m",sep = " "))
  }
  
  final <- rbind(final, sub.final)


#run null models for comparison
null <- matrix(NA, nrow = 1, ncol = 4)

  fit <- auto.arima(y, max.p = 1, max.P = 1)
  null[1,4] <- fit$aicc
  null[1,3] <- NA
  null[1,2] <- NA
  null[1,1] <- "Peak_depth_m"


final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  slice(which.min(AICc))

best.vars <- colnames(mydata4)[combn(cols,5)[,1]]
best.vars.cols <- combn(cols,5)[,1]

best.fit <- auto.arima(y,xreg = as.matrix(mydata4[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
accuracy(best.fit)
hist(unlist(mydata4[,2]))
plot(mydata4[,2],fitted(best.fit),xlim = c(-3,3),ylim = c(-3,3))
abline(a = 0, b = 1)
median((unlist(mydata4[,2])-unlist(fitted(best.fit))), na.rm = TRUE)

good <- final %>%
  filter(AICc > as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){
  good.vars.1 <- colnames(mydata4)[combn(cols,good[i,2])[,good[i,3]]]
  
  good.vars.1
  
  good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]
  
  
  good.fit.1 <- auto.arima(y,xreg = as.matrix(mydata4[,good.vars.cols.1]),max.p = 1, max.P = 1)
  print(good.fit.1)
  print(accuracy(good.fit.1))
  
  
}

#find best-fit ARIMA model for peak width
#peak width drivers = c(16:21,23:24,28,36)
colnames(mydata4)

cols <- c(9,11:14)
sub.final <- NULL
final <- NULL

y <- mydata4[,3]

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(mydata4[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      sub.sub.final[j,4] <- NA
      sub.sub.final[j,3] <- j
      sub.sub.final[j,2] <- i
      sub.sub.final[j,1] <- "Peak_width_m"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "Peak_width_m"
  }
  
  sub.final <- rbind(sub.final,sub.sub.final)
  print(paste("I have finished with all combinations of length",i,"for Peak_width_m",sep = " "))
}

final <- rbind(final, sub.final)


#run null models for comparison
null <- matrix(NA, nrow = 1, ncol = 4)

fit <- auto.arima(y, max.p = 1, max.P = 1)
null[1,4] <- fit$aicc
null[1,3] <- NA
null[1,2] <- NA
null[1,1] <- "Peak_width_m"


final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  slice(which.min(AICc))
best.vars <- colnames(mydata4)[combn(cols,3)[,3]]
best.vars.cols <- combn(cols,3)[,3]

best.fit <- auto.arima(y,xreg = as.matrix(mydata4[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
accuracy(best.fit)
hist(unlist(mydata4[,3]))
plot(mydata4[,3],fitted(best.fit),xlim = c(-3,3),ylim = c(-3,3))
abline(a = 0, b = 1)

median((unlist(mydata4[,3])-unlist(fitted(best.fit))), na.rm = TRUE)


good <- final %>%
  filter(AICc > as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){
good.vars.1 <- colnames(mydata4)[combn(cols,good[i,2])[,good[i,3]]]

good.vars.1

good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]


good.fit.1 <- auto.arima(y,xreg = as.matrix(mydata4[,good.vars.cols.1]),max.p = 1, max.P = 1)
print(good.fit.1)
print(accuracy(good.fit.1))



}

#find best-fit ARIMA model for max_biomass
colnames(mydata4)

cols <- c(9,15:18)
sub.final <- NULL
final <- NULL

y <- mydata4[,4]

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(mydata4[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      sub.sub.final[j,4] <- NA
      sub.sub.final[j,3] <- j
      sub.sub.final[j,2] <- i
      sub.sub.final[j,1] <- "Max_biomass_ugL"
      next }
    
    sub.sub.final[j,4] <- fit$aicc
    sub.sub.final[j,3] <- j
    sub.sub.final[j,2] <- i
    sub.sub.final[j,1] <- "Max_biomass_ugL"
  }
  
  sub.final <- rbind(sub.final,sub.sub.final)
  print(paste("I have finished with all combinations of length",i,"for Max_biomass_ugL",sep = " "))
}

final <- rbind(final, sub.final)


#run null models for comparison
null <- matrix(NA, nrow = 1, ncol = 4)

fit <- auto.arima(y, max.p = 1, max.P = 1)
null[1,4] <- fit$aicc
null[1,3] <- NA
null[1,2] <- NA
null[1,1] <- "Max_biomass_ugL"


final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  slice(which.min(AICc))
best.vars <- colnames(mydata4)[combn(cols,3)[,9]]
best.vars
best.vars.cols <- combn(cols,3)[,9]

best.fit <- auto.arima(y,xreg = as.matrix(mydata4[,best.vars.cols]),max.p = 1, max.P = 1)
best.fit
hist(resid(best.fit))
accuracy(best.fit)
hist(unlist(mydata4[,4]))
plot(mydata4[,4],fitted(best.fit),xlim = c(-3,3),ylim = c(-3,3))
abline(a = 0, b = 1)


good <- final %>%
  filter(AICc > as.numeric(best$AICc[1]) & AICc <= (as.numeric(best$AICc[1]) + 2)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(good)){
  good.vars.1 <- colnames(mydata4)[combn(cols,good[i,2])[,good[i,3]]]
  
  good.vars.1
  
  good.vars.cols.1 <- combn(cols,good[i,2])[,good[i,3]]
  
  
  good.fit.1 <- auto.arima(y,xreg = as.matrix(mydata4[,good.vars.cols.1]),max.p = 1, max.P = 1)
  print(good.fit.1)
  print(accuracy(good.fit.1))
  
  
}

####LOOKING AT MANIPULATION VS. REFERENCE YEARS####

#peak depth

#find best-fit ARIMA model for peak_depth_m

sub.final <- NULL
final <- NULL

years <- c("EM","noEM")

for(k in 1:length(years)){
  
  if(years[k] == "EM"){
    y <- dEM[,2]
    x <- dEM
    cols <- c(5:8,10:11)
    }
  if(years[k] == "noEM"){
    y <- dnoEM[,2]
    x <- dnoEM
    cols <- c(5:8,10:11)
    }
  

for (i in 1:length(cols)){
  my.combn <- combn(cols,i)
  sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 5)
  
  for (j in 1:ncol(my.combn)){
    
    skip_to_next <- FALSE
    
    tryCatch(fit <- auto.arima(y,xreg = as.matrix(x[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      sub.sub.final[j,5] <- NA
      sub.sub.final[j,4] <- j
      sub.sub.final[j,3] <- i
      sub.sub.final[j,2] <- years[k]
      sub.sub.final[j,1] <- "Peak_depth_m"
      next }
    
    sub.sub.final[j,5] <- fit$aicc
    sub.sub.final[j,4] <- j
    sub.sub.final[j,3] <- i
    sub.sub.final[j,2] <- years[k]
    sub.sub.final[j,1] <- "Peak_depth_m"
  }
  
  sub.final <- rbind(sub.final,sub.sub.final)
  print(paste("I have finished with all combinations of length",i,"for Peak_depth_m in",years[k],sep = " "))
}

final <- rbind(final, sub.final)

#run null models for comparison
null <- matrix(NA, nrow = 1, ncol = 5)

fit <- auto.arima(y, max.p = 1, max.P = 1)
null[1,5] <- fit$aicc
null[1,4] <- NA
null[1,3] <- NA
null[1,2] <- years[k]
null[1,1] <- "Peak_depth_m"


final <- rbind(final, null)

}


final <- data.frame(final)
colnames(final) <- c("Response.variable","Year","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  group_by(Year) %>%
  slice(which.min(AICc)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))

for (i in 1:nrow(best)){
  
  if(years[i] == "EM"){
    y <- dEM[,2]
    x <- dEM
    cols <- c(5:8,10:11)
  }
  if(years[i] == "noEM"){
    y <- dnoEM[,2]
    x <- dnoEM
    cols <- c(5:8,10:11)
  }
  
  best.vars <- colnames(x)[combn(cols,unlist(best[i,3]))[,unlist(best[i,4])]]
  best.vars.cols <- combn(cols,unlist(best[i,3]))[,unlist(best[i,4])]
  
  best.fit <- auto.arima(y,xreg = as.matrix(x[,best.vars.cols]),max.p = 1, max.P = 1)
  print(best[i,2])
  print(best.fit)
  print(accuracy(best.fit))

  good <- final %>%
    filter(Year == years[i]) %>%
    filter(AICc > as.numeric(best$AICc[i]) & AICc <= (as.numeric(best$AICc[i]) + 2)) %>%
    mutate(Num.covars = as.numeric(Num.covars),
           Covar.cols = as.numeric(Covar.cols))
  print(good)
  
  if(nrow(good)>0){
    for (j in 1:nrow(good)){
      good.vars.1 <- colnames(x)[combn(cols,unlist(good[j,3]))[,unlist(good[j,4])]]
      
      good.vars.1
      
      good.vars.cols.1 <- combn(cols,unlist(good[j,3]))[,unlist(good[j,4])]
      
      good.fit.1 <- auto.arima(y,xreg = as.matrix(x[,good.vars.cols.1]),max.p = 1, max.P = 1)
      print(good.fit.1)
      print(accuracy(good.fit.1))
      
      
    }}
  
  
}





#peak width

#find best-fit ARIMA model for peak_width_m

sub.final <- NULL
final <- NULL

years <- c("EM","noEM")

for(k in 1:length(years)){
  
  if(years[k] == "EM"){
    y <- dEM[,3]
    x <- dEM
    cols <- c(9,11:14)
  }
  if(years[k] == "noEM"){
    y <- dnoEM[,3]
    x <- dnoEM
    cols <- c(9,11:14)
  }
  
  
  
  for (i in 1:length(cols)){
    my.combn <- combn(cols,i)
    sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 5)
    
    for (j in 1:ncol(my.combn)){
      
      skip_to_next <- FALSE
      
      tryCatch(fit <- auto.arima(y,xreg = as.matrix(x[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
      
      if(skip_to_next) { 
        sub.sub.final[j,5] <- NA
        sub.sub.final[j,4] <- j
        sub.sub.final[j,3] <- i
        sub.sub.final[j,2] <- years[k]
        sub.sub.final[j,1] <- "Peak_width_m"
        next }
      
      sub.sub.final[j,5] <- fit$aicc
      sub.sub.final[j,4] <- j
      sub.sub.final[j,3] <- i
      sub.sub.final[j,2] <- years[k]
      sub.sub.final[j,1] <- "Peak_width_m"
    }
    
    sub.final <- rbind(sub.final,sub.sub.final)
    print(paste("I have finished with all combinations of length",i,"for Peak_width_m in",years[k],sep = " "))
  }
  
  final <- rbind(final, sub.final)
  
  #run null models for comparison
  null <- matrix(NA, nrow = 1, ncol = 5)
  
  fit <- auto.arima(y, max.p = 1, max.P = 1)
  null[1,5] <- fit$aicc
  null[1,4] <- NA
  null[1,3] <- NA
  null[1,2] <- years[k]
  null[1,1] <- "Peak_width_m"
  
  
  final <- rbind(final, null)
  
}


final <- data.frame(final)
colnames(final) <- c("Response.variable","Year","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  group_by(Year) %>%
  slice(which.min(AICc)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))


for (i in 1:nrow(best)){
  
  if(years[i] == "EM"){
    y <- dEM[,3]
    x <- dEM
    cols <- c(9,11:14)
  }
  if(years[i] == "noEM"){
    y <- dnoEM[,3]
    x <- dnoEM
    cols <- c(9,11:14)
  }
  
  best.vars <- colnames(x)[combn(cols,unlist(best[i,3]))[,unlist(best[i,4])]]
  best.vars.cols <- combn(cols,unlist(best[i,3]))[,unlist(best[i,4])]
  
  best.fit <- auto.arima(y,xreg = as.matrix(x[,best.vars.cols]),max.p = 1, max.P = 1)
  print(best[i,2])
  print(best.fit)
  print(accuracy(best.fit))

  good <- final %>%
    filter(Year == years[i]) %>%
    filter(AICc > as.numeric(best$AICc[i]) & AICc <= (as.numeric(best$AICc[i]) + 2)) %>%
    mutate(Num.covars = as.numeric(Num.covars),
           Covar.cols = as.numeric(Covar.cols))
  print(good)
  
  if(nrow(good)>0){
    for (j in 1:nrow(good)){
      good.vars.1 <- colnames(x)[combn(cols,unlist(good[j,3]))[,unlist(good[j,4])]]
      
      good.vars.1
      
      good.vars.cols.1 <- combn(cols,unlist(good[j,3]))[,unlist(good[j,4])]
      
      good.fit.1 <- auto.arima(y,xreg = as.matrix(x[,good.vars.cols.1]),max.p = 1, max.P = 1)
      print(good.fit.1)
      print(accuracy(good.fit.1))
      
      
      
    }}
  
}




#max biomass
#drop Temp_C from 2018
#drop Cmax_SRP_ugL and Cmax_DOC_mgL from 2019

#find best-fit ARIMA model for peak_depth_m

sub.final <- NULL
final <- NULL

years <- c("EM","noEM")

for(k in 1:length(years)){
  
  if(years[k] == "EM"){
    y <- dEM[,4]
    x <- dEM
    cols <- c(9,15:18)
  }
  if(years[k] == "noEM"){
    y <- dnoEM[,4]
    x <- dnoEM
    cols <- c(9,15:18)
  }
  
  
  for (i in 1:length(cols)){
    my.combn <- combn(cols,i)
    sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 5)
    
    for (j in 1:ncol(my.combn)){
      
      skip_to_next <- FALSE
      
      tryCatch(fit <- auto.arima(y,xreg = as.matrix(x[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
      
      if(skip_to_next) { 
        sub.sub.final[j,5] <- NA
        sub.sub.final[j,4] <- j
        sub.sub.final[j,3] <- i
        sub.sub.final[j,2] <- years[k]
        sub.sub.final[j,1] <- "Max_biomass_ugL"
        next }
      
      sub.sub.final[j,5] <- fit$aicc
      sub.sub.final[j,4] <- j
      sub.sub.final[j,3] <- i
      sub.sub.final[j,2] <- years[k]
      sub.sub.final[j,1] <- "Max_biomass_ugL"
    }
    
    sub.final <- rbind(sub.final,sub.sub.final)
    print(paste("I have finished with all combinations of length",i,"for Max_biomass_ugL in",years[k],sep = " "))
  }
  
  final <- rbind(final, sub.final)
  
  #run null models for comparison
  null <- matrix(NA, nrow = 1, ncol = 5)
  
  fit <- auto.arima(y, max.p = 1, max.P = 1)
  null[1,5] <- fit$aicc
  null[1,4] <- NA
  null[1,3] <- NA
  null[1,2] <- years[k]
  null[1,1] <- "Max_biomass_ugL"
  
  
  final <- rbind(final, null)
  
}


final <- data.frame(final)
colnames(final) <- c("Response.variable","Year","Num.covars","Covar.cols","AICc")
final <- distinct(final)

best <- final %>%
  group_by(Year) %>%
  slice(which.min(AICc)) %>%
  mutate(Num.covars = as.numeric(Num.covars),
         Covar.cols = as.numeric(Covar.cols))


for (i in 1:nrow(best)){
  
  if(years[i] == "EM"){
    y <- dEM[,4]
    x <- dEM
    cols <- c(9,15:18)
  }
  if(years[i] == "noEM"){
    y <- dnoEM[,4]
    x <- dnoEM
    cols <- c(9,15:18)
  }
  
  best.vars <- colnames(x)[combn(cols,unlist(best[i,3]))[,unlist(best[i,4])]]
  best.vars.cols <- combn(cols,unlist(best[i,3]))[,unlist(best[i,4])]
  
  best.fit <- auto.arima(y,xreg = as.matrix(x[,best.vars.cols]),max.p = 1, max.P = 1)
  print(best[i,2])
  print(best.fit)
  print(accuracy(best.fit))
  
  good <- final %>%
    filter(Year == years[i]) %>%
    filter(AICc > as.numeric(best$AICc[i]) & AICc <= (as.numeric(best$AICc[i]) + 2)) %>%
    mutate(Num.covars = as.numeric(Num.covars),
           Covar.cols = as.numeric(Covar.cols))
  print(good)
  
  if(nrow(good)>0){
    for (j in 1:nrow(good)){
      good.vars.1 <- colnames(x)[combn(cols,unlist(good[j,3]))[,unlist(good[j,4])]]
      
      good.vars.1
      
      good.vars.cols.1 <- combn(cols,unlist(good[j,3]))[,unlist(good[j,4])]
      
      good.fit.1 <- auto.arima(y,xreg = as.matrix(x[,good.vars.cols.1]),max.p = 1, max.P = 1)
      print(good.fit.1)
      print(accuracy(good.fit.1))
      
      
    }}
  
}
