#2A_ARIMA_analyses
#Author: Mary Lofton
#Date: 13OCT20

#need to check normality of all vars and correlation of drivers and fix as needed
pacman::p_load(PerformanceAnalytics, tidyverse, lubridate, forecast, utils, igraph,RColorBrewer)
rm(list=ls())


#get data
mydata <- read_csv("./2_Data_analysis/megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_Depth_m)
mydata <- mydata[,c(1,4,2:3,5:25,29:30,38:39,44)]

#look at correlations - cutoff of 0.5
driver_cor <- cor(mydata[,2:17],method = "spearman",use = "complete.obs")
driver_cor[lower.tri(driver_cor)] = ""
#thermo.depth and Year
#SRPmax_ugL and Cmax_SRP_ugL
#DINmax_ugL and DOC_max_mgL
#DOCmax_mgL and Cmax_DOC_mgL
#schmidt.stability and thermo.depth

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
  select(Max_biomass_ugL:BV_TOTAL,
         thermo.depth,Year,SRPmax_ugL,Cmax_SRP_ugL,DINmax_ugL,DOCmax_mgL,
         Cmax_DOC_mgL,schmidt.stability)
response_cor <- cor(mydata1,method = "spearman",use = "complete.obs")
response_cor <- data.frame(response_cor)
response_cor <- response_cor[,c(1:13)]
#thermo.depth and Year
mean(abs(as.numeric(response_cor[14,])))#thermo.depth = 0.13
mean(abs(as.numeric(response_cor[15,])))#year = 0.18
#SRPmax_ugL and Cmax_SRP_ugL
mean(abs(as.numeric(response_cor[16,])))#SRPmax_ugL = 0.10
mean(abs(as.numeric(response_cor[17,])))#Cmax_SRP_ugL = 0.19
#DINmax_ugL and DOC_max_mgL
mean(abs(as.numeric(response_cor[18,])))#DINmax_ugL = 0.28
mean(abs(as.numeric(response_cor[19,])))#DOCmax_ugL = 0.26
#DOCmax_mgL and Cmax_DOC_mgL
mean(abs(as.numeric(response_cor[20,])))#0.28
#schmidt.stability and thermo.depth
mean(abs(as.numeric(response_cor[21,])))#0.10
#get rid of thermo.depth
#get rid of SRPmax_ugL
#get rid of DOCmax_mgL

mydata2 <- mydata %>%
  select(-thermo.depth,-SRPmax_ugL,-DOCmax_mgL)

#check for skewness in drivers and responses and whether logging improves it
#ideally we want skewness to approach 0

for (i in 2:27){
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
#DINmax_ugL
#DOCmax_depth_m
#Max_biomass_ugL
#Peak_magnitude_ugL
#BV_Cyanobacteria
#rel_abund_Cryptophytes
#BV_TOTAL

mydata3 <- mydata2 %>%
  mutate(Cmax_DIN_ugL = log(Cmax_DIN_ugL),
         Cmax_DOC_mgL = log(Cmax_DOC_mgL),
         SRPmax_depth_m = log(SRPmax_depth_m),
         DINmax_ugL = log(DINmax_ugL),
         DOCmax_depth_m = log(DOCmax_depth_m),
         Max_biomass_ugL = log(Max_biomass_ugL),
         Peak_magnitude_ugL = log(Peak_magnitude_ugL),
         BV_Cyanobacteria = log(BV_Cyanobacteria+0.001),
         rel_abund_Cryptophytes = log(rel_abund_Cryptophytes+0.001),
         BV_TOTAL = log(BV_TOTAL))

#scale predictor/response variables to allow comparison of coefficients in model
mydata4 <- mydata3 %>%
  mutate_at(vars(-Date),scale) %>%
  select(-braycurtis, -jaccard)

#plot histograms to be sure what's going on here
for (i in 2:17){
  var <- mydata4[,i]
  hist(as.matrix(var),main = colnames(mydata4)[i])
}

#find best-fit ARIMA model for full timeseries 

cols <- c(2:14)
sub.final <- NULL
final <- NULL

for (k in 15:25){
  y <- mydata4[,k]

  for (i in 1:13){
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
null <- matrix(NA, nrow = 11, ncol = 4)

for (i in 15:25){
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
write.csv(final, "./2_Data_analysis/arima_results_all.csv",row.names = FALSE)

check <- final %>%
  filter(is.na(AICc))

#pull out best fit model and models <2 AICc of best fit
#also null and global models
#into table with covariates and coef values
cols <- c(2:14)
response.vars <- colnames(mydata4)[15:25]
top.models <- NULL
final <- read.csv("./2_Data_analysis/arima_results_all.csv") %>%
  filter(!is.na(Num.covars))

for (i in 1:11){
  
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

for (i in 1:11){
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

for (i in 1:11){
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
write.csv(top.models,"./2_Data_analysis/top_models_all.csv",row.names = FALSE)


# #find best-fit ARIMA model for 2016
# #need to check collinearity of drivers for each year separately?
# #then use most conservative base model for all years
# #document my decision tree
# data.2016 <- mydata[c(1:23),] %>%
#   select(-HOx, -Year) 
# 
# #look at correlations - cutoff of 0.5
# driver_cor <- cor(data.2016[,2:15],method = "spearman",use = "complete.obs")
# #DOC_max_depth_m and EM
# #schmidt.stability and EM
# #DINmax_depth_m and Cmax_SRP_ugL
# #DINmax_ugL and Cmax_DOC_mgL
# #DOCmax_mgL and Cmax_DOC_mgL
# #Temp_C and Cmax_DOC_mgL
# #schmidt.stability and SRPmax_ugL
# #DINmax_ugL and Kd
# #DINmax_ugL and DOC_max_depth_m
# #DINmax_ugL and DOCmax_mgL
# #DINmax_ugL and Temp_C
# #thermo.depth and DOCmax_depth_m
# #DOCmax_depth_m and Temp_C
# #DOCmax_ugL and Temp_C
# #Temp_C and thermo.depth
# #DOCmax_ugL and Kd
# 
# #look at correlations of collinear drivers w/ responses to select drivers
# mydata1 <- mydata %>%
#   select(Peak_width_m,Max_biomass_ugL,Peak_magnitude_ugL,Peak_depth_m,
#          thermo.depth,Year,SRPmax_ugL,Cmax_SRP_ugL,DINmax_ugL,DOCmax_mgL,
#          Cmax_DOC_mgL,schmidt.stability)
# response_cor <- cor(mydata1,method = "spearman",use = "complete.obs")
# #get rid of Year
# #get rid of SRPmax_ugL
# #get rid of DOCmax_mgL
# #get rid of schmidt.stability
# 
# mydata2 <- mydata %>%
#   select(-Year,-SRPmax_ugL,-DOCmax_mgL,-schmidt.stability)
# 
# #check for skewness in drivers and responses and whether logging improves it
# #ideally we want skewness to approach 0
# 
# for (i in 2:17){
#   print(colnames(mydata2)[i])
#   var <- mydata2[,i]
#   hist(as.matrix(var), main = colnames(mydata2)[i])
#   print(skewness(mydata2[,i], na.rm = TRUE))
#   print(skewness(log(mydata2[,i]), na.rm = TRUE))
#   var <- log(mydata2[,i])
#   hist(as.matrix(var), main = c("Log",colnames(mydata2)[i]))
# }
# #should log:
# #Cmax_DIN_ugL
# #Cmax_DOC_mgL
# #DINmax_ugL
# #DOCmax_depth_m
# #Max_biomass_ugL
# #Peak_magnitude_ugL
# 
# mydata3 <- mydata2 %>%
#   mutate(Cmax_DIN_ugL = log(Cmax_DIN_ugL),
#          Cmax_DOC_mgL = log(Cmax_DOC_mgL),
#          DINmax_ugL = log(DINmax_ugL),
#          DOCmax_depth_m = log(DOCmax_depth_m),
#          Max_biomass_ugL = log(Max_biomass_ugL),
#          Peak_magnitude_ugL = log(Peak_magnitude_ugL))
# 
# #scale predictor variables to allow comparison of coefficients in model
# mydata4 <- mydata3 %>%
#   mutate_at(vars(-Date,-Max_biomass_ugL,-Peak_width_m,-Peak_depth_m,-Peak_magnitude_ugL),scale)
# 
# #plot histograms to be sure what's going on here
# for (i in 2:17){
#   var <- mydata4[,i]
#   hist(as.matrix(var),main = colnames(mydata4)[i])
# }
# 
# cols <- c(2:12)
# sub.final.2016 <- NULL
# final.2016 <- NULL
# 
# for (k in 13:16){
#   y <- data.2016[,k]
#   
#   for (i in 1:11){
#     my.combn <- combn(cols,i)
#     sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
#     
#     for (j in 1:ncol(my.combn)){
#       
#       fit <- auto.arima(y,xreg = as.matrix(data.2016[,my.combn[,j]]))
#       sub.sub.final[j,4] <- fit$aicc
#       sub.sub.final[j,3] <- j
#       sub.sub.final[j,2] <- i
#       sub.sub.final[j,1] <- colnames(data.2016)[k]
#     }
#     
#     sub.final.2016 <- rbind(sub.final.2016,sub.sub.final)
#     print(paste("I have finished with all combinations of length",i,"for",colnames(data.2016)[k],sep = " "))
#   }
#   
#   final.2016 <- rbind(final.2016, sub.final.2016)
# }
# 
# #run null models for comparison
# null.2016 <- matrix(NA, nrow = 4, ncol = 4)
# 
# for (i in 13:16){
#   y <- data.2016[,i]
#   fit <- auto.arima(y)
#   null.2016[i-12,4] <- fit$aicc
#   null.2016[i-12,3] <- NA
#   null.2016[i-12,2] <- NA
#   null.2016[i-12,1] <- colnames(data.2016)[i]
# }
# 
# final.2016 <- rbind(final.2016, null.2016)
# final.2016 <- data.frame(final.2016)
# colnames(final.2016) <- c("Response.variable","Num.covars","Covar.cols","AICc")
# final.2016 <- distinct(final.2016)
# write.csv(final.2016, "./2_Data_analysis/arima_results_2016.csv",row.names = FALSE)
# 
# #pull out best fit model and models <2 AICc of best fit
# #also null and global models
# #into table with covariates and coef values
# response.vars <- c("Max_biomass_ugL","Peak_depth_m","Peak_magnitude_ugL","Peak_width_m")
# top.models <- NULL
# 
# for (i in 1:4){
#   
#   focal.var <- subset(final.2016, final.2016$Response.variable == response.vars[i])
#   focal.var[,2] <- as.numeric(focal.var[,2])
#   focal.var[,3] <- as.numeric(focal.var[,3])
#   focal.var[,4] <- as.numeric(focal.var[,4])
#   
#   best.model <- subset(focal.var, focal.var$AICc == min(focal.var$AICc))
#   best.fit <- auto.arima(data.2016[,response.vars[i]],xreg = as.matrix(data.2016[,combn(cols,best.model[1,2])[,best.model[1,3]]]))
#   
#   model.info <- matrix(NA, nrow = length(names(best.fit$coef)), ncol = 9)
#   model.info[,1] <- response.vars[i]
#   model.info[,2] <- 1
#   model.info[,3] <- best.fit$aicc
#   model.info[,4] <- 0
#   model.info[,5] <- best.fit$sigma2
#   model.info[,6] <- best.fit$loglik
#   model.info[,7] <- names(best.fit$coef)
#   model.info[,8] <- unname(best.fit$coef)
#   model.info[,9] <- unname(sqrt(diag(vcov(best.fit))))
#   
#   good.models <- subset(focal.var, abs(as.numeric(focal.var$AICc) - as.numeric(min(focal.var$AICc)))<2) %>%
#     distinct() %>%
#     arrange(AICc)
#   
#   if(length(good.models[,1])>1){
#   
#   for (k in 2:length(good.models[,1])){
#     
#     current.model <- subset(focal.var, focal.var$AICc == good.models$AICc[k])
#     current.fit <- auto.arima(data.2016[,response.vars[i]],xreg = as.matrix(data.2016[,combn(cols,current.model[1,2])[,current.model[1,3]]]))
#     
#     current.model.info <- matrix(NA, nrow = length(names(current.fit$coef)), ncol = 9)
#     current.model.info[,1] <- response.vars[i]
#     current.model.info[,2] <- k
#     current.model.info[,3] <- current.fit$aicc
#     current.model.info[,4] <- abs(best.fit$aicc-current.fit$aicc)
#     current.model.info[,5] <- current.fit$sigma2
#     current.model.info[,6] <- current.fit$loglik
#     current.model.info[,7] <- names(current.fit$coef)
#     current.model.info[,8] <- unname(current.fit$coef)
#     current.model.info[,9] <- unname(sqrt(diag(vcov(current.fit))))
#     
#     model.info <- rbind(model.info, current.model.info)
#     
#   }}
#   
#   top.models <- rbind(top.models, model.info)
# }
# 
# #null models
# #run null models for comparison
# null.models <- NULL
# 
# for (i in 1:4){
#   focal.var <- data.2016[,response.vars[i]]
#   null.fit <- auto.arima(focal.var)
#   
#   null.model.info <- matrix(NA, nrow = length(names(null.fit$coef)), ncol = 9)
#   null.model.info[,1] <- response.vars[i]
#   null.model.info[,2] <- "null"
#   null.model.info[,3] <- null.fit$aicc
#   null.model.info[,4] <- abs(best.fit$aicc-null.fit$aicc)
#   null.model.info[,5] <- null.fit$sigma2
#   null.model.info[,6] <- null.fit$loglik
#   null.model.info[,7] <- names(null.fit$coef)
#   null.model.info[,8] <- unname(null.fit$coef)
#   null.model.info[,9] <- unname(sqrt(diag(vcov(null.fit))))
#   
#   null.models <- rbind(null.models, null.model.info)
#   
# }
# 
# top.models <- rbind(top.models,null.models)
# 
# #global models
# #run global models for comparison
# global.models <- NULL
# 
# for (i in 1:4){
#   focal.var <- data.2016[,response.vars[i]]
#   global.fit <- auto.arima(focal.var, xreg = as.matrix(data.2016[,cols]))
#   
#   global.model.info <- matrix(NA, nrow = length(names(global.fit$coef)), ncol = 9)
#   global.model.info[,1] <- response.vars[i]
#   global.model.info[,2] <- "global"
#   global.model.info[,3] <- global.fit$aicc
#   global.model.info[,4] <- abs(best.fit$aicc-global.fit$aicc)
#   global.model.info[,5] <- global.fit$sigma2
#   global.model.info[,6] <- global.fit$loglik
#   global.model.info[,7] <- names(global.fit$coef)
#   global.model.info[,8] <- unname(global.fit$coef)
#   global.model.info[,9] <- unname(sqrt(diag(vcov(global.fit))))
#   
#   global.models <- rbind(global.models, global.model.info)
#   
# }
# 
# top.models <- rbind(top.models,global.models)
# 
# top.models <- data.frame(top.models)
# top.models[,3] <- round(as.numeric(top.models[,3]),2)
# top.models[,4] <- round(as.numeric(top.models[,4]),2)
# top.models[,5] <- round(as.numeric(top.models[,5]),2)
# top.models[,6] <- round(as.numeric(top.models[,6]),2)
# top.models[,8] <- round(as.numeric(top.models[,8]),2)
# top.models[,9] <- round(as.numeric(top.models[,9]),2)
# colnames(top.models) <- c("Response.var","Rank","AICc","Del.AICc","Sigma.2","Log.likelihood","Covar",
#                           "Covar.coefs","Covar.coef.SE")
# write.csv(top.models,"./2_Data_analysis/top_models_2016.csv",row.names = FALSE)
# 
# 
# 
# #find best-fit ARIMA model for 2017
# cols <- c(2:12)
# data.2017 <- mydata2[c(25:46),] %>%
#   select(-HOx) %>%
#   mutate_at(vars(-Date,-Max_biomass_ugL,-Peak_width_m,-Peak_depth_m,-Peak_magnitude_ugL),scale) 
# for (i in 2:12){
#   var <- data.2017[,i]
#   hist(as.matrix(var),main = colnames(data.2017)[i])
# }
# sub.final.2017 <- NULL
# final.2017 <- NULL
# 
# for (k in 13:16){
#   y <- data.2017[,k]
# 
#   for (i in 1:11){
#     my.combn <- combn(cols,i)
#     sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
#     
#     for (j in 1:ncol(my.combn)){
#       fit <- auto.arima(y,xreg = as.matrix(data.2017[,my.combn[,j]]))
#       sub.sub.final[j,4] <- fit$aicc
#       sub.sub.final[j,3] <- j
#       sub.sub.final[j,2] <- i
#       sub.sub.final[j,1] <- colnames(data.2017)[k]
#     }
#     
#     sub.final.2017 <- rbind(sub.final.2017,sub.sub.final)
#     print(paste("I have finished with all combinations of length",i,"for",colnames(data.2016)[k],sep = " "))
#   }
#   
#   final.2017 <- rbind(final.2017, sub.final.2017)
# }
# 
# #run null models for comparison
# null.2017 <- matrix(NA, nrow = 4, ncol = 4)
# 
# for (i in 13:16){
#   y <- data.2017[,i]
#   fit <- auto.arima(y)
#   null.2017[i-12,4] <- fit$aicc
#   null.2017[i-12,3] <- NA
#   null.2017[i-12,2] <- NA
#   null.2017[i-12,1] <- colnames(data.2017)[i]
# }
# 
# final.2017 <- rbind(final.2017, null.2017)
# final.2017 <- data.frame(final.2017)
# colnames(final.2017) <- c("Response.variable","Num.covars","Covar.cols","AICc")
# final.2017 <- distinct(final.2017)
# write.csv(final.2017, "./2_Data_analysis/arima_results_2017.csv",row.names = FALSE)
# 
# #pull out best fit model and models <2 AICc of best fit
# #also null and global models
# #into table with covariates and coef values
# response.vars <- c("Max_biomass_ugL","Peak_depth_m","Peak_magnitude_ugL","Peak_width_m")
# top.models <- NULL
# 
# for (i in 1:4){
#   
#   focal.var <- subset(final.2017, final.2017$Response.variable == response.vars[i])
#   focal.var[,2] <- as.numeric(focal.var[,2])
#   focal.var[,3] <- as.numeric(focal.var[,3])
#   focal.var[,4] <- as.numeric(focal.var[,4])
#   
#   best.model <- subset(focal.var, focal.var$AICc == min(focal.var$AICc))
#   best.fit <- auto.arima(data.2017[,response.vars[i]],xreg = as.matrix(data.2017[,combn(cols,best.model[1,2])[,best.model[1,3]]]))
#   
#   model.info <- matrix(NA, nrow = length(names(best.fit$coef)), ncol = 9)
#   model.info[,1] <- response.vars[i]
#   model.info[,2] <- 1
#   model.info[,3] <- best.fit$aicc
#   model.info[,4] <- 0
#   model.info[,5] <- best.fit$sigma2
#   model.info[,6] <- best.fit$loglik
#   model.info[,7] <- names(best.fit$coef)
#   model.info[,8] <- unname(best.fit$coef)
#   model.info[,9] <- unname(sqrt(diag(vcov(best.fit))))
#   
#   good.models <- subset(focal.var, abs(as.numeric(focal.var$AICc) - as.numeric(min(focal.var$AICc)))<2) %>%
#     distinct() %>%
#     arrange(AICc)
#   
#   if(length(good.models[,1])>1){
#   
#   for (k in 2:length(good.models[,1])){
#     
#     current.model <- subset(focal.var, focal.var$AICc == good.models$AICc[k])
#     current.fit <- auto.arima(data.2017[,response.vars[i]],xreg = as.matrix(data.2017[,combn(cols,current.model[1,2])[,current.model[1,3]]]))
#     
#     current.model.info <- matrix(NA, nrow = length(names(current.fit$coef)), ncol = 9)
#     current.model.info[,1] <- response.vars[i]
#     current.model.info[,2] <- k
#     current.model.info[,3] <- current.fit$aicc
#     current.model.info[,4] <- abs(best.fit$aicc-current.fit$aicc)
#     current.model.info[,5] <- current.fit$sigma2
#     current.model.info[,6] <- current.fit$loglik
#     current.model.info[,7] <- names(current.fit$coef)
#     current.model.info[,8] <- unname(current.fit$coef)
#     current.model.info[,9] <- unname(sqrt(diag(vcov(current.fit))))
#     
#     model.info <- rbind(model.info, current.model.info)
#     
#   }}
#   
#   top.models <- rbind(top.models, model.info)
# }
# 
# #null models
# #run null models for comparison
# null.models <- NULL
# 
# for (i in 1:4){
#   focal.var <- data.2017[,response.vars[i]]
#   null.fit <- auto.arima(focal.var)
#   
#   null.model.info <- matrix(NA, nrow = length(names(null.fit$coef)), ncol = 9)
#   null.model.info[,1] <- response.vars[i]
#   null.model.info[,2] <- "null"
#   null.model.info[,3] <- null.fit$aicc
#   null.model.info[,4] <- abs(best.fit$aicc-null.fit$aicc)
#   null.model.info[,5] <- null.fit$sigma2
#   null.model.info[,6] <- null.fit$loglik
#   null.model.info[,7] <- names(null.fit$coef)
#   null.model.info[,8] <- unname(null.fit$coef)
#   null.model.info[,9] <- unname(sqrt(diag(vcov(null.fit))))
#   
#   null.models <- rbind(null.models, null.model.info)
#   
# }
# 
# top.models <- rbind(top.models,null.models)
# 
# #global models
# #run global models for comparison
# global.models <- NULL
# 
# for (i in 1:4){
#   focal.var <- data.2017[,response.vars[i]]
#   global.fit <- auto.arima(focal.var, xreg = as.matrix(data.2017[,cols]))
#   
#   global.model.info <- matrix(NA, nrow = length(names(global.fit$coef)), ncol = 9)
#   global.model.info[,1] <- response.vars[i]
#   global.model.info[,2] <- "global"
#   global.model.info[,3] <- global.fit$aicc
#   global.model.info[,4] <- abs(best.fit$aicc-global.fit$aicc)
#   global.model.info[,5] <- global.fit$sigma2
#   global.model.info[,6] <- global.fit$loglik
#   global.model.info[,7] <- names(global.fit$coef)
#   global.model.info[,8] <- unname(global.fit$coef)
#   global.model.info[,9] <- unname(sqrt(diag(vcov(global.fit))))
#   
#   global.models <- rbind(global.models, global.model.info)
#   
# }
# 
# top.models <- rbind(top.models,global.models)
# 
# top.models <- data.frame(top.models)
# top.models[,3] <- round(as.numeric(top.models[,3]),2)
# top.models[,4] <- round(as.numeric(top.models[,4]),2)
# top.models[,5] <- round(as.numeric(top.models[,5]),2)
# top.models[,6] <- round(as.numeric(top.models[,6]),2)
# top.models[,8] <- round(as.numeric(top.models[,8]),2)
# top.models[,9] <- round(as.numeric(top.models[,9]),2)
# colnames(top.models) <- c("Response.var","Rank","AICc","Del.AICc","Sigma.2","Log.likelihood","Covar",
#                           "Covar.coefs","Covar.coef.SE")
# write.csv(top.models,"./2_Data_analysis/top_models_2017.csv",row.names = FALSE)
# 
# 
# 
# #find best-fit ARIMA model for 2018
# 
# cols <- c(2:12)
# data.2018 <- mydata4[c(48:67),] %>%
#   select(-EM) %>%
#   mutate_at(vars(-Date,-Max_biomass_ugL,-Peak_width_m,-Peak_depth_m,-Peak_magnitude_ugL),scale) 
# for (i in 2:12){
#   var <- data.2018[,i]
#   hist(as.matrix(var),main = colnames(data.2018)[i])
# }
# sub.final.2018 <- NULL
# final.2018 <- NULL
# 
# for (k in 13:16){
#   y <- data.2018[,k]
#   
#   for (i in 1:11){
#     my.combn <- combn(cols,i)
#     sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
#     
#     for (j in 1:ncol(my.combn)){
#       fit <- auto.arima(y,xreg = as.matrix(data.2018[,my.combn[,j]]))
#       sub.sub.final[j,4] <- fit$aicc
#       sub.sub.final[j,3] <- j
#       sub.sub.final[j,2] <- i
#       sub.sub.final[j,1] <- colnames(data.2018)[k]
#     }
#     
#     sub.final.2018 <- rbind(sub.final.2018,sub.sub.final)
#     print(paste("I have finished with all combinations of length",i,"for",colnames(data.2018)[k],sep = " "))
#   }
#   
#   final.2018 <- rbind(final.2018, sub.final.2018)
# }
# 
# #run null models for comparison
# null.2018 <- matrix(NA, nrow = 4, ncol = 4)
# 
# for (i in 13:16){
#   y <- data.2018[,i]
#   fit <- auto.arima(y)
#   null.2018[i-12,4] <- fit$aicc
#   null.2018[i-12,3] <- NA
#   null.2018[i-12,2] <- NA
#   null.2018[i-12,1] <- colnames(data.2018)[i]
# }
# 
# final.2018 <- rbind(final.2018, null.2018)
# final.2018 <- data.frame(final.2018)
# colnames(final.2018) <- c("Response.variable","Num.covars","Covar.cols","AICc")
# final.2018 <- distinct(final.2018)
# write.csv(final.2018, "./2_Data_analysis/arima_results_2018.csv",row.names = FALSE)
# 
# #pull out best fit model and models <2 AICc of best fit
# #also null and global models
# #into table with covariates and coef values
# response.vars <- c("Max_biomass_ugL","Peak_depth_m","Peak_magnitude_ugL","Peak_width_m")
# top.models <- NULL
# 
# for (i in 1:4){
#   
#   focal.var <- subset(final.2018, final.2018$Response.variable == response.vars[i])
#   focal.var[,2] <- as.numeric(focal.var[,2])
#   focal.var[,3] <- as.numeric(focal.var[,3])
#   focal.var[,4] <- as.numeric(focal.var[,4])
#   
#   best.model <- subset(focal.var, focal.var$AICc == min(focal.var$AICc))
#   best.fit <- auto.arima(data.2018[,response.vars[i]],xreg = as.matrix(data.2018[,combn(cols,best.model[1,2])[,best.model[1,3]]]))
#   
#   model.info <- matrix(NA, nrow = length(names(best.fit$coef)), ncol = 9)
#   model.info[,1] <- response.vars[i]
#   model.info[,2] <- 1
#   model.info[,3] <- best.fit$aicc
#   model.info[,4] <- 0
#   model.info[,5] <- best.fit$sigma2
#   model.info[,6] <- best.fit$loglik
#   model.info[,7] <- names(best.fit$coef)
#   model.info[,8] <- unname(best.fit$coef)
#   model.info[,9] <- unname(sqrt(diag(vcov(best.fit))))
#   
#   good.models <- subset(focal.var, abs(as.numeric(focal.var$AICc) - as.numeric(min(focal.var$AICc)))<2) %>%
#     distinct() %>%
#     arrange(AICc)
#   
#   if(length(good.models[,1])>1){
#   
#   for (k in 2:length(good.models[,1])){
#     
#     current.model <- subset(focal.var, focal.var$AICc == good.models$AICc[k])
#     current.fit <- auto.arima(data.2018[,response.vars[i]],xreg = as.matrix(data.2018[,combn(cols,current.model[1,2])[,current.model[1,3]]]))
#     
#     current.model.info <- matrix(NA, nrow = length(names(current.fit$coef)), ncol = 9)
#     current.model.info[,1] <- response.vars[i]
#     current.model.info[,2] <- k
#     current.model.info[,3] <- current.fit$aicc
#     current.model.info[,4] <- abs(best.fit$aicc-current.fit$aicc)
#     current.model.info[,5] <- current.fit$sigma2
#     current.model.info[,6] <- current.fit$loglik
#     current.model.info[,7] <- names(current.fit$coef)
#     current.model.info[,8] <- unname(current.fit$coef)
#     current.model.info[,9] <- unname(sqrt(diag(vcov(current.fit))))
#     
#     model.info <- rbind(model.info, current.model.info)
#     
#   }}
#   
#   top.models <- rbind(top.models, model.info)
# }
# 
# #null models
# #run null models for comparison
# null.models <- NULL
# 
# for (i in 1:4){
#   focal.var <- data.2018[,response.vars[i]]
#   null.fit <- auto.arima(focal.var)
#   
#   null.model.info <- matrix(NA, nrow = length(names(null.fit$coef)), ncol = 9)
#   null.model.info[,1] <- response.vars[i]
#   null.model.info[,2] <- "null"
#   null.model.info[,3] <- null.fit$aicc
#   null.model.info[,4] <- abs(best.fit$aicc-null.fit$aicc)
#   null.model.info[,5] <- null.fit$sigma2
#   null.model.info[,6] <- null.fit$loglik
#   null.model.info[,7] <- names(null.fit$coef)
#   null.model.info[,8] <- unname(null.fit$coef)
#   null.model.info[,9] <- unname(sqrt(diag(vcov(null.fit))))
#   
#   null.models <- rbind(null.models, null.model.info)
#   
# }
# 
# top.models <- rbind(top.models,null.models)
# 
# #global models
# #run global models for comparison
# global.models <- NULL
# 
# for (i in 1:4){
#   focal.var <- data.2018[,response.vars[i]]
#   global.fit <- auto.arima(focal.var, xreg = as.matrix(data.2018[,cols]))
#   
#   global.model.info <- matrix(NA, nrow = length(names(global.fit$coef)), ncol = 9)
#   global.model.info[,1] <- response.vars[i]
#   global.model.info[,2] <- "global"
#   global.model.info[,3] <- global.fit$aicc
#   global.model.info[,4] <- abs(best.fit$aicc-global.fit$aicc)
#   global.model.info[,5] <- global.fit$sigma2
#   global.model.info[,6] <- global.fit$loglik
#   global.model.info[,7] <- names(global.fit$coef)
#   global.model.info[,8] <- unname(global.fit$coef)
#   global.model.info[,9] <- unname(sqrt(diag(vcov(global.fit))))
#   
#   global.models <- rbind(global.models, global.model.info)
#   
# }
# 
# top.models <- rbind(top.models,global.models)
# 
# top.models <- data.frame(top.models)
# top.models[,3] <- round(as.numeric(top.models[,3]),2)
# top.models[,4] <- round(as.numeric(top.models[,4]),2)
# top.models[,5] <- round(as.numeric(top.models[,5]),2)
# top.models[,6] <- round(as.numeric(top.models[,6]),2)
# top.models[,8] <- round(as.numeric(top.models[,8]),2)
# top.models[,9] <- round(as.numeric(top.models[,9]),2)
# colnames(top.models) <- c("Response.var","Rank","AICc","Del.AICc","Sigma.2","Log.likelihood","Covar",
#                           "Covar.coefs","Covar.coef.SE")
# write.csv(top.models,"./2_Data_analysis/top_models_2018.csv",row.names = FALSE)
# 
# 
# 
# 
# #find best-fit ARIMA model for 2019
# cols <- c(2:12)
# data.2019 <- mydata4[c(75:104),] %>%
#   select(-EM) %>%
#   mutate_at(vars(-Date,-Max_biomass_ugL,-Peak_width_m,-Peak_depth_m,-Peak_magnitude_ugL),scale) 
# for (i in 2:12){
#   var <- data.2019[,i]
#   hist(as.matrix(var),main = colnames(data.2019)[i])
# }
# sub.final.2019 <- NULL
# final.2019 <- NULL
# 
# for (k in 13:16){
#   y <- data.2019[,k]
#   
#   for (i in 1:11){
#     my.combn <- combn(cols,i)
#     sub.sub.final <- matrix(NA, nrow = ncol(my.combn), ncol = 4)
#     
#     for (j in 1:ncol(my.combn)){
#       
#       skip_to_next <- FALSE
#       
#       tryCatch(fit <- auto.arima(y,xreg = as.matrix(data.2019[,my.combn[,j]])), error = function(e) { skip_to_next <<- TRUE})
#       
#       if(skip_to_next) { next } 
#       
#       sub.sub.final[j,4] <- fit$aicc
#       sub.sub.final[j,3] <- j
#       sub.sub.final[j,2] <- i
#       sub.sub.final[j,1] <- colnames(data.2019)[k]
#     }
#     
#     sub.final.2019 <- rbind(sub.final.2019,sub.sub.final)
#     print(paste("I have finished with all combinations of length",i,"for",colnames(data.2019)[k],sep = " "))
#   }
#   
#   final.2019 <- rbind(final.2019, sub.final.2019)
# }
# 
# #run null models for comparison
# null.2019 <- matrix(NA, nrow = 4, ncol = 4)
# 
# for (i in 13:16){
#   y <- data.2019[,i]
#   fit <- auto.arima(y)
#   null.2019[i-12,4] <- fit$aicc
#   null.2019[i-12,3] <- NA
#   null.2019[i-12,2] <- NA
#   null.2019[i-12,1] <- colnames(data.2019)[i]
# }
# 
# final.2019 <- rbind(final.2019, null.2019)
# final.2019 <- data.frame(final.2019)
# colnames(final.2019) <- c("Response.variable","Num.covars","Covar.cols","AICc")
# final.2019 <- distinct(final.2019)
# write.csv(final.2019, "./2_Data_analysis/arima_results_2019.csv",row.names = FALSE)
# 
# #pull out best fit model and models <2 AICc of best fit
# #also null and global models
# #into table with covariates and coef values
# response.vars <- c("Max_biomass_ugL","Peak_depth_m","Peak_magnitude_ugL","Peak_width_m")
# top.models <- NULL
# 
# for (i in 1:4){
#   
#   focal.var <- subset(final.2019, final.2019$Response.variable == response.vars[i])
#   focal.var[,2] <- as.numeric(focal.var[,2])
#   focal.var[,3] <- as.numeric(focal.var[,3])
#   focal.var[,4] <- as.numeric(focal.var[,4])
#   
#   best.model <- subset(focal.var, focal.var$AICc == min(focal.var$AICc))
#   best.fit <- auto.arima(data.2019[,response.vars[i]],xreg = as.matrix(data.2019[,combn(cols,best.model[1,2])[,best.model[1,3]]]))
#   
#   model.info <- matrix(NA, nrow = length(names(best.fit$coef)), ncol = 9)
#   model.info[,1] <- response.vars[i]
#   model.info[,2] <- 1
#   model.info[,3] <- best.fit$aicc
#   model.info[,4] <- 0
#   model.info[,5] <- best.fit$sigma2
#   model.info[,6] <- best.fit$loglik
#   model.info[,7] <- names(best.fit$coef)
#   model.info[,8] <- unname(best.fit$coef)
#   model.info[,9] <- unname(sqrt(diag(vcov(best.fit))))
#   
#   good.models <- subset(focal.var, abs(as.numeric(focal.var$AICc) - as.numeric(min(focal.var$AICc)))<2) %>%
#     distinct() %>%
#     arrange(AICc)
#   
#   if(length(good.models[,1])>1){
#   
#   for (k in 2:length(good.models[,1])){
#     
#     current.model <- subset(focal.var, focal.var$AICc == good.models$AICc[k])
#     current.fit <- auto.arima(data.2019[,response.vars[i]],xreg = as.matrix(data.2019[,combn(cols,current.model[1,2])[,current.model[1,3]]]))
#     
#     current.model.info <- matrix(NA, nrow = length(names(current.fit$coef)), ncol = 9)
#     current.model.info[,1] <- response.vars[i]
#     current.model.info[,2] <- k
#     current.model.info[,3] <- current.fit$aicc
#     current.model.info[,4] <- abs(best.fit$aicc-current.fit$aicc)
#     current.model.info[,5] <- current.fit$sigma2
#     current.model.info[,6] <- current.fit$loglik
#     current.model.info[,7] <- names(current.fit$coef)
#     current.model.info[,8] <- unname(current.fit$coef)
#     current.model.info[,9] <- unname(sqrt(diag(vcov(current.fit))))
#     
#     model.info <- rbind(model.info, current.model.info)
#     
#   }}
#   
#   top.models <- rbind(top.models, model.info)
# }
# 
# #null models
# #run null models for comparison
# null.models <- NULL
# 
# for (i in 1:4){
#   focal.var <- data.2019[,response.vars[i]]
#   null.fit <- auto.arima(focal.var)
#   
#   null.model.info <- matrix(NA, nrow = length(names(null.fit$coef)), ncol = 9)
#   null.model.info[,1] <- response.vars[i]
#   null.model.info[,2] <- "null"
#   null.model.info[,3] <- null.fit$aicc
#   null.model.info[,4] <- abs(best.fit$aicc-null.fit$aicc)
#   null.model.info[,5] <- null.fit$sigma2
#   null.model.info[,6] <- null.fit$loglik
#   null.model.info[,7] <- names(null.fit$coef)
#   null.model.info[,8] <- unname(null.fit$coef)
#   null.model.info[,9] <- unname(sqrt(diag(vcov(null.fit))))
#   
#   null.models <- rbind(null.models, null.model.info)
#   
# }
# 
# top.models <- rbind(top.models,null.models)
# 
# #global models
# #run global models for comparison
# global.models <- NULL
# 
# for (i in 1:4){
#   focal.var <- data.2019[,response.vars[i]]
#   global.fit <- auto.arima(focal.var, xreg = as.matrix(data.2019[,cols]))
#   
#   global.model.info <- matrix(NA, nrow = length(names(global.fit$coef)), ncol = 9)
#   global.model.info[,1] <- response.vars[i]
#   global.model.info[,2] <- "global"
#   global.model.info[,3] <- global.fit$aicc
#   global.model.info[,4] <- abs(best.fit$aicc-global.fit$aicc)
#   global.model.info[,5] <- global.fit$sigma2
#   global.model.info[,6] <- global.fit$loglik
#   global.model.info[,7] <- names(global.fit$coef)
#   global.model.info[,8] <- unname(global.fit$coef)
#   global.model.info[,9] <- unname(sqrt(diag(vcov(global.fit))))
#   
#   global.models <- rbind(global.models, global.model.info)
#   
# }
# 
# top.models <- rbind(top.models,global.models)
# 
# top.models <- data.frame(top.models)
# top.models[,3] <- round(as.numeric(top.models[,3]),2)
# top.models[,4] <- round(as.numeric(top.models[,4]),2)
# top.models[,5] <- round(as.numeric(top.models[,5]),2)
# top.models[,6] <- round(as.numeric(top.models[,6]),2)
# top.models[,8] <- round(as.numeric(top.models[,8]),2)
# top.models[,9] <- round(as.numeric(top.models[,9]),2)
# colnames(top.models) <- c("Response.var","Rank","AICc","Del.AICc","Sigma.2","Log.likelihood","Covar",
#                           "Covar.coefs","Covar.coef.SE")
# write.csv(top.models,"./2_Data_analysis/top_models_2019.csv",row.names = FALSE)
# 
# 
# driver_cor <- cor(data.2016[,my.combn[,j]],method = "spearman",use = "complete.obs")
