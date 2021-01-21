#2B_phyto_comm_ARIMA_analyses
#Author: Mary Lofton
#Date: 13OCT20

#need to check normality of all vars and correlation of drivers and fix as needed
pacman::p_load(PerformanceAnalytics, tidyverse, lubridate, forecast, utils, igraph,RColorBrewer)
rm(list=ls())


#get data
mydata <- read_csv("./2_Data_analysis/CS_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20")

#look at correlations - cutoff of 0.5
driver_cor <- cor(mydata[,2:18],method = "spearman",use = "complete.obs")
driver_cor[lower.tri(driver_cor)] = ""
#thermo.depth and Year - yes
#SRPmax_ugL and Cmax_SRP_ugL - yes
#DINmax_ugL and DOC_max_mgL - yes
#DOCmax_mgL and Cmax_DOC_mgL - yes
#schmidt.stability and thermo.depth - yes
#thermo.depth and perc_light_thermocline - yes

#look at correlations of collinear drivers w/ responses to select drivers
mydata1 <- mydata %>%
  select(shannon:BV_TOTAL,
         thermo.depth,Year,SRPmax_ugL,Cmax_SRP_ugL,DINmax_ugL,DOCmax_mgL,
         Cmax_DOC_mgL,schmidt.stability, perc_light_thermocline)
response_cor <- cor(mydata1,method = "spearman",use = "complete.obs")
response_cor <- data.frame(response_cor)
response_cor <- response_cor[c(10:18),c(1:9)]
#thermo.depth and Year
mean(abs(as.numeric(response_cor[1,])))#thermo.depth = 0.07
mean(abs(as.numeric(response_cor[2,])))#year = 0.28
mean(abs(as.numeric(response_cor[8,])))#schmidt.stability = 0.10
mean(abs(as.numeric(response_cor[9,])))#perc_light_thermocline = 0.20
#SRPmax_ugL and Cmax_SRP_ugL
mean(abs(as.numeric(response_cor[3,])))#SRPmax_ugL = 0.13
mean(abs(as.numeric(response_cor[4,])))#Cmax_SRP_ugL = 0.09
#DINmax_ugL and DOCmax_mgL
mean(abs(as.numeric(response_cor[5,])))#DINmax_ugL = 0.30
mean(abs(as.numeric(response_cor[6,])))#DOCmax_ugL = 0.28
#DOCmax_mgL and Cmax_DOC_mgL
mean(abs(as.numeric(response_cor[7,])))#0.19

#get rid of thermo.depth
#get rid of Cmax_SRP_ugL
#get rid of DOCmax_mgL

mydata2 <- mydata %>%
  select(-thermo.depth, -Cmax_SRP_ugL, -DOCmax_mgL)

#check for skewness in drivers and responses and whether logging improves it
#ideally we want skewness to approach 0

for (i in 2:24){
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
#Cmax_DOC_mgL
#SRPmax_depth_m
#SRPmax_ugL
#DINmax_ugL
#DOCmax_depth_m
#perc_light_thermocline
#BV_Bacillaria
#BV_Cyanobacteria
#rel_abund_Bacillaria
#rel_abund_Cryptophytes
#BV_TOTAL

mydata3 <- mydata2 %>%
  mutate(Cmax_DIN_ugL = log(Cmax_DIN_ugL),
         Cmax_DOC_mgL = log(Cmax_DOC_mgL),
         SRPmax_depth_m = log(SRPmax_depth_m),
         SRPmax_ugL = log(SRPmax_ugL),
         DINmax_ugL = log(DINmax_ugL),
         DOCmax_depth_m = log(DOCmax_depth_m),
         perc_light_thermocline = log(perc_light_thermocline),
         BV_Bacillaria = log(BV_Bacillaria+0.001),
         BV_Cyanobacteria = log(BV_Cyanobacteria+0.001),
         rel_abund_Bacillaria = log(rel_abund_Bacillaria+0.001),
         rel_abund_Cryptophytes = log(rel_abund_Cryptophytes+0.001),
         BV_TOTAL = log(BV_TOTAL+0.001))


#scale predictor/response variables to allow comparison of coefficients in model
mydata4 <- mydata3 %>%
  mutate_at(vars(-Date),scale) 

#plot histograms to be sure what's going on here
for (i in 2:24){
  var <- mydata4[,i]
  hist(as.matrix(var),main = colnames(mydata4)[i])
}

#find best-fit ARIMA model for full timeseries 

cols <- c(2:12)
sub.final <- NULL
final <- NULL

for (k in 13:16){
  y <- mydata4[,k]

  for (i in 1:11){
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
null <- matrix(NA, nrow = 4, ncol = 4)

for (i in 13:16){
  y <- mydata4[,i]
  fit <- auto.arima(y, max.p = 1, max.P = 1)
  null[i-14,4] <- fit$aicc
  null[i-14,3] <- NA
  null[i-14,2] <- NA
  null[i-14,1] <- colnames(mydata4)[i]
}

final <- rbind(final, null)
final <- data.frame(final)
colnames(final) <- c("Response.variable","Num.covars","Covar.cols","AICc")
final <- distinct(final)
write.csv(final, "./2_Data_analysis/FP_arima_results.csv",row.names = FALSE)

check <- final %>%
  filter(is.na(AICc))

#pull out best fit model and models <2 AICc of best fit
#also null and global models
#into table with covariates and coef values
cols <- c(2:12)
response.vars <- colnames(mydata4)[13:16]
top.models <- NULL
final <- read.csv("./2_Data_analysis/FP_arima_results.csv") %>%
  filter(!is.na(Num.covars))

for (i in 1:4){
  
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

for (i in 1:4){
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

for (i in 1:4){
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
write.csv(top.models,"./2_Data_analysis/FP_top_models.csv",row.names = FALSE)