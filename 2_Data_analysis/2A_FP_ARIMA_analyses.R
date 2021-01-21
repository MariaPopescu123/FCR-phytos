#2A_FP_ARIMA_analyses
#Author: Mary Lofton
#Date: 13OCT20

#need to check normality of all vars and correlation of drivers and fix as needed
pacman::p_load(PerformanceAnalytics, tidyverse, lubridate, forecast, utils, igraph,RColorBrewer)
rm(list=ls())


#get data
mydata <- read_csv("./2_Data_analysis/FP_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") %>%
  select(-MonthDay)
mydata <- mydata[,c(2,1,3:22)]

AR_NAs <- read_csv("./2_Data_analysis/FP_megamatrix.csv")%>%
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

#transform correlation matrix into adjacency matrix for network plot
driver_cor <- data.frame(driver_cor)

dc <- driver_cor %>%
  mutate(source = row.names(driver_cor)) %>%
  gather(Year:Kd, key = "target", value = "importance") %>%
  filter(!importance %in% c(1,0,"")) %>%
  mutate(importance = round(abs(as.numeric(importance)*10),0)) %>%
  filter(!importance == 0) %>%
  filter(importance >=5) %>%
  mutate(importance = (importance-5)*3)

names <- colnames(mydata)[2:17]

#set up nodes
nodes <- data.frame(
  name=names,
  carac=c("year",rep("mgmt",2),rep("nutrients",9),rep("physics",4)))

# build the network graph object
network <- graph_from_data_frame(d=dc, vertices=nodes, directed=F)
# Make a palette of 3 colors
coul  <- brewer.pal(4, "Set1") 
# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(network)$carac))]
#plot
resfactor = 3
png(filename = "./3_Visualization/driver_cor.png", 
    height = 600*resfactor, width = 600*resfactor, res = 72*resfactor)
plot(network, vertex.color=my_color,edge.width=E(network)$importance,
     layout = layout.circle,
     vertex.size = 35,
     vertex.label.color="black",
     vertex.label.cex = 0.9)
dev.off()

data.2016 <- mydata[c(1:23),] %>%
  select(-HOx, -Year) 

#look at correlations in 2016 - cutoff of 0.5
driver_cor.2016 <- cor(data.2016[,2:15],method = "spearman",use = "complete.obs")
driver_cor.2016[lower.tri(driver_cor.2016)] = ""

#DOCmax_depth_m and EM
#schmidt.stability and EM
#DINmax_depth_m and Cmax_SRP_ugL
#DINmax_ugL and Cmax_DOC_mgL
#DOCmax_mgL and Cmax_DOC_mgL
#Temp_C and Cmax_DOC_mgL
#schmidt.stability and SRPmax_ugL
#DINmax_ugL and Kd
#DINmax_ugL and DOCmax_depth_m
#DINmax_ugL and DOCmax_mgL
#DINmax_ugL and Temp_C
#thermo.depth and DOCmax_depth_m
#DOCmax_depth_m and Temp_C
#DOCmax_ugL and Temp_C
#Temp_C and thermo.depth
#DOCmax_ugL and Kd

#transform correlation matrix into adjacency matrix for network plot
driver_cor.2016 <- data.frame(driver_cor.2016)

dc_2016 <- driver_cor.2016 %>%
  mutate(source = row.names(driver_cor.2016)) %>%
  gather(EM:Kd, key = "target", value = "importance") %>%
  filter(!importance %in% c(1,0,"")) %>%
  mutate(importance = round(abs(as.numeric(importance)*10),0)) %>%
  filter(!importance == 0) %>%
  filter(importance >=5) %>%
  mutate(importance = (importance-5)*3)

names <- colnames(data.2016)[2:15]

#set up nodes
nodes <- data.frame(
  name=names,
  carac=c("mgmt",rep("nutrients",9),rep("physics",4)))

# build the network graph object
network.2016 <- graph_from_data_frame(d=dc_2016, vertices=nodes, directed=F)
# Make a palette of 3 colors
coul  <- brewer.pal(3, "Set1") 
# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(network.2017)$carac))]
#plot
resfactor = 3
png(filename = "./3_Visualization/driver_cor_2016.png", 
    height = 600*resfactor, width = 600*resfactor, res = 72*resfactor)
plot(network.2016, vertex.color=my_color,edge.width=E(network.2016)$importance,
     layout = layout.circle,
     vertex.size = 35,
     vertex.label.color="black",
     vertex.label.cex = 0.9)
dev.off()

data.2017 <- mydata[c(25:46),] %>%
  select(-HOx, -Year)

#look at correlations in 2017 - cutoff of 0.5
driver_cor.2017 <- cor(data.2017[,c(2:15)],method = "spearman",use = "complete.obs")
driver_cor.2017[lower.tri(driver_cor.2017)] = ""
#EM and Cmax_SRP_ugL
#Cmax_SRP_ugL and Temp_C
#Cmax_DOC_mgL and Temp_C
#Cmax_DOC_mgL and schmidt.stability
#SRPmax_depth_m and DOCmax_depth_m
#SRPmax_depth_m and Temp_C
#SRPmax_ugL and DOC_max_mgL
#DINmax_depth_m and schmidt.stability
#DINmax_depth_m and thermo.depth
#DINmax_depth_m and Kd
#DINmax_ugL and DOCmax_depth_m
#Temp_C and schmidt.stability
#thermo.depth and Kd

#transform correlation matrix into adjacency matrix for network plot
driver_cor.2017 <- data.frame(driver_cor.2017)

dc_2017 <- driver_cor.2017 %>%
  mutate(source = row.names(driver_cor.2017)) %>%
  gather(EM:Kd, key = "target", value = "importance") %>%
  filter(!importance %in% c(1,0,"")) %>%
  mutate(importance = round(abs(as.numeric(importance)*10),0)) %>%
  filter(!importance == 0) %>%
  filter(importance >=5) %>%
  mutate(importance = (importance-5)*3)

names <- colnames(data.2017)[2:15]

#set up nodes
nodes <- data.frame(
  name=names,
  carac=c("mgmt",rep("nutrients",9),rep("physics",4)))

# build the network graph object
network.2017 <- graph_from_data_frame(d=dc_2017, vertices=nodes, directed=F)
# Make a palette of 3 colors
coul  <- brewer.pal(3, "Set1") 
# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(network.2017)$carac))]
#plot
resfactor = 3
png(filename = "./3_Visualization/driver_cor_2017.png", 
    height = 600*resfactor, width = 600*resfactor, res = 72*resfactor)
plot(network.2017, vertex.color=my_color,edge.width=E(network.2017)$importance,
     layout = layout.circle,
     vertex.size = 35,
     vertex.label.color="black",
     vertex.label.cex = 0.9)
dev.off()

data.2018 <- mydata[c(48:67),] %>%
  select(-EM,-Year)

#look at correlations in 2018 - cutoff of 0.5
driver_cor.2018 <- cor(data.2018[,c(2:15)],method = "spearman",use = "complete.obs")
driver_cor.2018[lower.tri(driver_cor.2018)] = ""
#HOx and Cmax_DIN_ugL
#HOx and Cmax_DOC_mgL
#HOx and SRPmax_depth_m
#HOx and DINmax_depth_m
#HOx and DINmax_ugL
#HOx and DOCmax_depth_m
#HOx and DOCmax_mgL
#HOx and thermo.depth
#Cmax_SRP_ugL and SRPmax_ugL
#Cmax_SRP_ugL and Kd
#Cmax_DIN_ugL and DINmax_ugL
#Cmax_DIN_ugL and Temp_C
#Cmax_DIN_ugL and schmidt.stability
#Cmax_DOC_mgL and DINmax_depth_m
#Cmax_DOC_mgL and DINmax_ugL
#Cmax_DOC_mgL and DOCmax_depth_m
#Cmax_DOC_mgL and DOCmax_mgL
#Cmax_DOC_mgL and Temp_C
#Cmax_DOC_mgL and thermo.depth
#SRP_max_depth_m and DINmax_ugL
#SRP_max_depth_m and schmidt.stability
#SRPmax_ugL and Kd
#DINmax_depth_m and DINmax_ugL
#DINmax_depth_m and DOCmax_depth_m
#DINmax_depth_m and DOCmax_mgL
#DINmax_depth_m and thermo.depth
#DINmax_ugL and DOCmax_mgL
#DINmax_ugL and schmidt.stability
#DINmax_ugL and thermo.depth
#DOCmax_depth_m and DOCmax_mgL
#DOCmax_depth_m and thermo.depth
#DOCmax_mgL and thermo.depth
#schmidt.stability and thermo.depth

#transform correlation matrix into adjacency matrix for network plot
driver_cor.2018 <- data.frame(driver_cor.2018)

dc_2018 <- driver_cor.2018 %>%
  mutate(source = row.names(driver_cor.2018)) %>%
  gather(HOx:Kd, key = "target", value = "importance") %>%
  filter(!importance %in% c(1,0,"")) %>%
  mutate(importance = round(abs(as.numeric(importance)*10),0)) %>%
  filter(!importance == 0) %>%
  filter(importance >=5) %>%
  mutate(importance = (importance-5)*3)

names <- colnames(data.2018)[2:15]

#set up nodes
nodes <- data.frame(
  name=names,
  carac=c("mgmt",rep("nutrients",9),rep("physics",4)))

# build the network graph object
network.2018 <- graph_from_data_frame(d=dc_2018, vertices=nodes, directed=F)
# Make a palette of 3 colors
coul  <- brewer.pal(3, "Set1") 
# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(network.2018)$carac))]
#plot
resfactor = 3
png(filename = "./3_Visualization/driver_cor_2018.png", 
    height = 600*resfactor, width = 600*resfactor, res = 72*resfactor)
plot(network.2018, vertex.color=my_color,edge.width=E(network.2018)$importance,
     layout = layout.circle,
     vertex.size = 35,
     vertex.label.color="black",
     vertex.label.cex = 0.9)
dev.off()

data.2019 <- mydata[c(75:104),] %>%
  select(-EM,-Year)

#look at correlations in 2019 - cutoff of 0.5
driver_cor.2019 <- cor(data.2019[,c(2:15)],method = "spearman",use = "complete.obs")
driver_cor.2019[lower.tri(driver_cor.2019)] = ""
#HOx and SRPmax_ugL
#Cmax_SRP_ugL and Cmax_DOC_mgL
#Cmax_SRP_ugL and SRPmax_depth_m
#Cmax_SRP_ugL and SRPmax_ugL
#Cmax_SRP_ugL and DOCmax_mgL
#Cmax_SRP_ugL and Temp_C
#Cmax_DIN_ugL and SRPmax_depth_m
#Cmax_DIN_ugL and thermo.depth
#Cmax_DOC_mgL and SRPmax_depth_m
#Cmax_DOC_mgL and DOCmax_mgL
#Cmax_DOC_mgL and thermo.depth
#SRPmax_depth_m and DINmax_ugL
#SRPmax_depth_m and thermo.depth
#SRPmax_ugL and Temp_C
#DINmax_ugL and DOCmax_mgL
#DINmax_ugL and thermo.depth
#Temp_C and schmidt.stability
#Temp_C and Kd
#schmidt.stability and Kd

#transform correlation matrix into adjacency matrix for network plot
driver_cor.2019 <- data.frame(driver_cor.2019)

dc_2019 <- driver_cor.2019 %>%
  mutate(source = row.names(driver_cor.2019)) %>%
  gather(HOx:Kd, key = "target", value = "importance") %>%
  filter(!importance %in% c(1,0,"")) %>%
  mutate(importance = round(abs(as.numeric(importance)*10),0)) %>%
  filter(!importance == 0) %>%
  filter(importance >=5) %>%
  mutate(importance = (importance-5)*3)

names <- colnames(data.2019)[2:15]

#set up nodes
nodes <- data.frame(
  name=names,
  carac=c("mgmt",rep("nutrients",9),rep("physics",4)))

# build the network graph object
network.2019 <- graph_from_data_frame(d=dc_2019, vertices=nodes, directed=F)
# Make a palette of 3 colors
coul  <- brewer.pal(3, "Set1") 
# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(network.2019)$carac))]
#plot
resfactor = 3
png(filename = "./3_Visualization/driver_cor_2019.png", 
    height = 600*resfactor, width = 600*resfactor, res = 72*resfactor)
plot(network.2019, vertex.color=my_color,edge.width=E(network.2019)$importance,
     layout = layout.circle,
     vertex.size = 35,
     vertex.label.color="black",
     vertex.label.cex = 0.9)
dev.off()



#look at correlations of collinear drivers w/ responses to select drivers
mydata1 <- mydata %>%
  select(Max_biomass_ugL:Peak_width_m,
         Cmax_DOC_mgL, Year, SRPmax_ugL, Cmax_SRP_ugL, DOCmax_mgL, DINmax_ugL,
         schmidt.stability, thermo.depth, perc_light_thermocline)
response_cor <- cor(mydata1,method = "spearman",use = "complete.obs")
response_cor <- data.frame(response_cor)
response_cor <- response_cor[c(5:13),c(1:4)]

#Cmax_DOC_mgL and Year
mean(abs(as.numeric(response_cor[1,])))#Cmax_DOC_mgL = 0.25
mean(abs(as.numeric(response_cor[2,])))#Year = 0.15

#SRPmax_ugL and Cmax_SRP_ugL
mean(abs(as.numeric(response_cor[3,])))#SRPmax_ugL = 0.05
mean(abs(as.numeric(response_cor[4,])))#Cmax_SRP_ugL = 0.15

#DOCmax_mgL and Year
mean(abs(as.numeric(response_cor[5,])))#DOCmax_mgL = 0.33

#DOCmax_mgL and Cmax_DOC_mgL

#DOCmax_mgL and DINmax_ugL
mean(abs(as.numeric(response_cor[6,])))#DINmax_ugL = 0.29

#schmidt.stability and Year
mean(abs(as.numeric(response_cor[7,])))#schmidt.stability = 0.19

#thermo.depth and Year
mean(abs(as.numeric(response_cor[8,])))#thermo.depth = 0.27

#thermo.depth and schmidt.stability

#perc_light_thermocline and thermo.depth
mean(abs(as.numeric(response_cor[9,])))#perc_light_thermocline = 0.20

#get rid of Year, schmidt.stability, and perc_light_thermocline
#get rid of SRPmax_ugL
#get rid of DINmax_ugL and Cmax_DOC_mgL

mydata2 <- mydata %>%
  select(-Year, -schmidt.stability, -perc_light_thermocline,
         -SRPmax_ugL, -DINmax_ugL, -Cmax_DOC_mgL)

#check for skewness in drivers and responses and whether logging improves it
#ideally we want skewness to approach 0

for (i in 2:16){
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
#DOCmax_depth_m
#DOCmax_mgL
#Max_biomass_ugL
#Peak_magnitude_ugL

mydata3 <- mydata2 %>%
  mutate(Cmax_DIN_ugL = log(Cmax_DIN_ugL),
         SRPmax_depth_m = log(SRPmax_depth_m),
         DOCmax_depth_m = log(DOCmax_depth_m),
         DOCmax_mgL = log(DOCmax_mgL),
         Max_biomass_ugL = log(Max_biomass_ugL),
         Peak_magnitude_ugL = log(Peak_magnitude_ugL))

#scale predictor/response variables to allow comparison of coefficients in model
mydata4 <- mydata3 %>%
  mutate_at(vars(-Date),scale) 

#plot histograms to be sure what's going on here
for (i in 2:16){
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
