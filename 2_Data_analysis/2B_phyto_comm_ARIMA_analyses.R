#2B_phyto_comm_ARIMA_analyses
#Author: Mary Lofton
#Date: 21JAN21

#need to check normality of all vars and correlation of drivers and fix as needed
pacman::p_load(PerformanceAnalytics, tidyverse, lubridate, forecast, utils, igraph,RColorBrewer)
rm(list=ls())


#get data
mydata <- read_csv("./2_Data_analysis/CS_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") %>%
  select(-MonthDay)

AR_NAs <- read_csv("./2_Data_analysis/CS_megamatrix.csv")%>%
  select(-Chem_Depth_m,-Temp_Depth_m)
AR_NAs <- AR_NAs[c(24,47,68),]

mydata <- bind_rows(mydata, AR_NAs)%>%
  arrange(Date)

#look at correlations - cutoff of 0.5
driver_cor <- cor(mydata[,2:18],method = "spearman",use = "complete.obs")
driver_cor[lower.tri(driver_cor)] = ""
#Cmax_DOC_mgL and Year
#SRPmax_ugL and Cmax_SRP_ugL
#DOCmax_mgL and Year
#DOCmax_mgL and Cmax_DOC_mgL
#DOCmax_mgL and DINmax_ugL
#schmidt.stability and Year
#thermo.depth and Year
#thermo.depth and schmidt.stability
#perc_light_thermocline and thermo.depth

#look at correlations of collinear drivers w/ responses to select drivers
mydata1 <- mydata %>%
  select(shannon:BV_TOTAL,
         Cmax_DOC_mgL, Year, SRPmax_ugL, Cmax_SRP_ugL, DOCmax_mgL, DINmax_ugL,
         schmidt.stability, thermo.depth, perc_light_thermocline)
response_cor <- cor(mydata1,method = "spearman",use = "complete.obs")
response_cor <- data.frame(response_cor)
response_cor <- response_cor[c(10:18),c(1:9)]

#Cmax_DOC_mgL and Year
mean(abs(as.numeric(response_cor[1,])))#Cmax_DOC_mgL = 0.20
mean(abs(as.numeric(response_cor[2,])))#Year = 0.25

#SRPmax_ugL and Cmax_SRP_ugL
mean(abs(as.numeric(response_cor[3,])))#SRPmax_ugL = 0.15
mean(abs(as.numeric(response_cor[4,])))#Cmax_SRP_ugL = 0.09

#DOCmax_mgL and Year
mean(abs(as.numeric(response_cor[5,])))#DOCmax_mgL = 0.27

#DOCmax_mgL and Cmax_DOC_mgL

#DOCmax_mgL and DINmax_ugL
mean(abs(as.numeric(response_cor[6,])))#DINmax_ugL = 0.28

#schmidt.stability and Year
mean(abs(as.numeric(response_cor[7,])))#schmidt.stability = 0.16

#thermo.depth and Year
mean(abs(as.numeric(response_cor[8,])))#thermo.depth = 0.06

#thermo.depth and schmidt.stability

#perc_light_thermocline and thermo.depth
mean(abs(as.numeric(response_cor[9,])))#perc_light_thermocline = 0.19

#OPTION 1
#drop Cmax_DOC_mgL = 0.20
#drop Cmax_SRP_ugL = 0.09
#drop DOCmax_mgL = 0.27
#drop thermo.depth = 0.06
#drop schmidt.stability = 0.16

#OR

#OPTION 2
#drop Year = 0.25
#drop Cmax_SRP_ugL = 0.09
#drop DOCmax_mgL = 0.27
#drop thermo.depth = 0.06

#GOING WITH OPTION 1; I think it better maximizes keeping in strong predictors

mydata2 <- mydata %>%
  select(-Cmax_DOC_mgL, -Cmax_SRP_ugL, -DOCmax_mgL, -thermo.depth,
         -schmidt.stability)

#check for skewness in drivers and responses and whether logging improves it
#ideally we want skewness to approach 0

for (i in 2:22){
  print(colnames(mydata2)[i])
  var <- mydata2[,i]
  hist(as.matrix(var), main = colnames(mydata2)[i])
  print(skewness(mydata2[,i], na.rm = TRUE))
  print(skewness(log(mydata2[,i]+0.0001), na.rm = TRUE))
  var <- log(mydata2[,i])
  hist(as.matrix(var), main = c("Log",colnames(mydata2)[i]))
}
#should log:
#Cmax_DIN_ugL
#SRPmax_depth_m
#SRPmax_ugL
#DINmax_ugL
#DOCmax_depth_m
#Kd
#perc_light_thermocline
#BV_Cyanobacteria
#rel_abund_Bacillaria
#BV_TOTAL

mydata3 <- mydata2 %>%
  mutate(Cmax_DIN_ugL = log(Cmax_DIN_ugL),
         SRPmax_depth_m = log(SRPmax_depth_m),
         SRPmax_ugL = log(SRPmax_ugL),
         DINmax_ugL = log(DINmax_ugL),
         DOCmax_depth_m = log(DOCmax_depth_m),
         Kd = log(Kd),
         perc_light_thermocline = log(perc_light_thermocline),
         BV_Cyanobacteria = log(BV_Cyanobacteria+0.001),
         rel_abund_Bacillaria = log(rel_abund_Bacillaria+0.001),
         BV_TOTAL = log(BV_TOTAL))

#scale predictor/response variables to allow comparison of coefficients in model
mydata4 <- mydata3 %>%
  mutate_at(vars(-Date),scale) 

#plot histograms to be sure what's going on here
for (i in 2:22){
  var <- mydata4[,i]
  hist(as.matrix(var),main = colnames(mydata4)[i])
}

#find best-fit ARIMA model for full timeseries 

cols <- c(2:13)
sub.final <- NULL
final <- NULL

for (k in 14:22){
  y <- mydata4[,k]
  
  for (i in 1:12){
    my.combn <- combn(cols,i)
    sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
    
    for (j in 1:ncol(my.combn)){
      
      skip_to_next <- FALSE
      
      tryCatch(fit <- auto.arima(y,xreg = as.matrix(mydata4[,my.combn[,j]]),max.p = 1, max.P = 1), error = function(e) { skip_to_next <<- TRUE})
      
      if(skip_to_next) { 
        sub.sub.final[j,4] <- NA
        sub.sub.final[j,3] <- j
        sub.sub.final[j,2] <- i
        sub.sub.final[j,1] <- colnames(mydata4)[k]
        next }
      
      sub.sub.final[j,4] <- fit$aicc
      sub.sub.final[j,3] <- j
      sub.sub.final[j,2] <- i
      sub.sub.final[j,1] <- colnames(mydata4)[k]
    }
    
    sub.final <- rbind(sub.final,sub.sub.final)
    print(paste("I have finished with all combinations of length",i,"for",colnames(mydata4)[k],sep = " "))
  }
  
  final <- rbind(final, sub.final)
}

#run null models for comparison
null <- matrix(NA, nrow = 9, ncol = 4)

for (i in 14:22){
  y <- mydata4[,i]
  fit <- auto.arima(y, max.p = 1, max.P = 1)
  null[i-13,4] <- fit$aicc
  null[i-13,3] <- NA
  null[i-13,2] <- NA
  null[i-13,1] <- colnames(mydata4)[i]
}

final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)
write.csv(final, "./2_Data_analysis/CS_arima_results.csv",row.names = FALSE)

check <- final %>%
  filter(is.na(AICc))

#pull out best fit model and models <2 AICc of best fit
#also null and global models
#into table with covariates and coef values
cols <- c(2:13)
response.vars <- colnames(mydata4)[14:22]
top.models <- NULL
final <- read.csv("./2_Data_analysis/CS_arima_results.csv") %>%
  filter(!is.na(Num.covars))

for (i in 1:9){
  
  focal.var <- subset(final, final$Response.variable == response.vars[i])
  
  best.model <- subset(focal.var, focal.var$AICc == min(focal.var$AICc, na.rm = TRUE))
  best.fit <- auto.arima(mydata4[,response.vars[i]],xreg = as.matrix(mydata4[,combn(cols,unlist(best.model[1,2]))[,best.model[1,3]]]), max.p = 1, max.P = 1)
  
  model.info <- matrix(NA, nrow = length(names(best.fit$coef)), ncol = 9)
  model.info[,1] <- response.vars[i]
  model.info[,2] <- 1
  model.info[,3] <- best.fit$aicc
  model.info[,4] <- 0
  model.info[,5] <- best.fit$sigma2
  model.info[,6] <- best.fit$loglik
  model.info[,7] <- names(best.fit$coef)
  model.info[,8] <- unname(best.fit$coef)
  model.info[,9] <- unname(sqrt(diag(vcov(best.fit))))
  
  good.models <- subset(focal.var, abs(as.numeric(focal.var$AICc) - as.numeric(min(focal.var$AICc)))<2) %>%
    distinct() %>%
    arrange(AICc)
  
  if(length(good.models[,1])>1){
    
    for (k in 2:length(good.models[,1])){
      
      current.model <- subset(focal.var, focal.var$AICc == good.models$AICc[k])
      current.fit <- auto.arima(mydata4[,response.vars[i]],xreg = as.matrix(mydata4[,combn(cols,current.model[1,2])[,current.model[1,3]]]), max.p = 1, max.P = 1)
      
      current.model.info <- matrix(NA, nrow = length(names(current.fit$coef)), ncol = 9)
      current.model.info[,1] <- response.vars[i]
      current.model.info[,2] <- k
      current.model.info[,3] <- current.fit$aicc
      current.model.info[,4] <- abs(best.fit$aicc-current.fit$aicc)
      current.model.info[,5] <- current.fit$sigma2
      current.model.info[,6] <- current.fit$loglik
      current.model.info[,7] <- names(current.fit$coef)
      current.model.info[,8] <- unname(current.fit$coef)
      current.model.info[,9] <- unname(sqrt(diag(vcov(current.fit))))
      
      model.info <- rbind(model.info, current.model.info)
      
    }}
  
  top.models <- rbind(top.models, model.info)
}

#null models
#run null models for comparison
null.models <- NULL

for (i in 1:9){
  focal.var <- mydata4[,response.vars[i]]
  null.fit <- auto.arima(focal.var, max.p = 1, max.P = 1)
  
  null.model.info <- matrix(NA, nrow = length(names(null.fit$coef)), ncol = 9)
  null.model.info[,1] <- response.vars[i]
  null.model.info[,2] <- "null"
  null.model.info[,3] <- null.fit$aicc
  null.model.info[,4] <- abs(best.fit$aicc-null.fit$aicc)
  null.model.info[,5] <- null.fit$sigma2
  null.model.info[,6] <- null.fit$loglik
  null.model.info[,7] <- names(null.fit$coef)
  null.model.info[,8] <- unname(null.fit$coef)
  null.model.info[,9] <- unname(sqrt(diag(vcov(null.fit))))
  
  null.models <- rbind(null.models, null.model.info)
  
}

top.models <- rbind(top.models,null.models)

#global models
#run global models for comparison
global.models <- NULL

for (i in 1:9){
  focal.var <- mydata4[,response.vars[i]]
  global.fit <- auto.arima(focal.var, xreg = as.matrix(mydata4[,cols]),max.p = 1, max.P = 1)
  
  global.model.info <- matrix(NA, nrow = length(names(global.fit$coef)), ncol = 9)
  global.model.info[,1] <- response.vars[i]
  global.model.info[,2] <- "global"
  global.model.info[,3] <- global.fit$aicc
  global.model.info[,4] <- abs(best.fit$aicc-global.fit$aicc)
  global.model.info[,5] <- global.fit$sigma2
  global.model.info[,6] <- global.fit$loglik
  global.model.info[,7] <- names(global.fit$coef)
  global.model.info[,8] <- unname(global.fit$coef)
  global.model.info[,9] <- unname(sqrt(diag(vcov(global.fit))))
  
  global.models <- rbind(global.models, global.model.info)
  
}

top.models <- rbind(top.models,global.models)

top.models <- data.frame(top.models)
top.models[,3] <- round(as.numeric(top.models[,3]),2)
top.models[,4] <- round(as.numeric(top.models[,4]),2)
top.models[,5] <- round(as.numeric(top.models[,5]),2)
top.models[,6] <- round(as.numeric(top.models[,6]),2)
top.models[,8] <- round(as.numeric(top.models[,8]),2)
top.models[,9] <- round(as.numeric(top.models[,9]),2)
colnames(top.models) <- c("Response.var","Rank","AICc","Del.AICc","Sigma.2","Log.likelihood","Covar",
                          "Covar.coefs","Covar.coef.SE")
write.csv(top.models,"./2_Data_analysis/CS_top_models.csv",row.names = FALSE)
