#Peak_depth_correlations
#Author: Mary Lofton
#Date: 22FEB21

#load packages
#install.packages('pacman')
pacman::p_load(tidyverse, lubridate, data.table)
rm(list=ls())
options(scipen = 999)

#read in FP data and drivers
fp <- read_csv("./2_Data_analysis/FP_megamatrix.csv")
colnames(fp)
fp <- fp[,c(1:6,11,13,15,19:20,23,26,31:33,35,28,27)] %>%
  mutate(WRT_day = log(WRT_day))

output <- matrix(NA, nrow = length(c(7:17)), ncol = 4)

for (i in 7:17){
  
  x <- unlist(fp[,i])
  y <- fp$Peak_depth_m
  
  output[i-6,1] <- colnames(fp)[i]
  output[i-6,2] <- round(cor(x,y, method = "spearman", use = "complete.obs"),2)
  spearman_p <- cor.test(x,y, method = "spearman", use = "complete.obs", exact = F)
  output[i-6,3] <- spearman_p$p.value
  
}

output[,4] <- p.adjust(output[,3], method = "holm")

output <- data.frame(output)
colnames(output) <- c("Driver","Spearmans_rho","Spearmans_p","adj_Spearmans_p")

#For all years: schmidt.stability, pz_depth_m, pz_Temp_C, pz_pH
drivers <- fp %>%
  select(pz_depth_m, pz_Temp_C, Peak_depth_m)

ggplot(data = drivers, aes(x = pz_depth_m, y = Peak_depth_m))+
  geom_point(size = 2)+
  theme_classic()

ggplot(data = drivers, aes(x = pz_Temp_C, y = Peak_depth_m))+
  geom_point(size = 2)+
  theme_classic()

##for 2019
d2019 <- fp %>% filter(Year == 2019)
output <- matrix(NA, nrow = length(c(7:17)), ncol = 4)

for (i in 7:17){
  
  x <- unlist(d2019[,i])
  y <- d2019$Peak_depth_m
  
  output[i-6,1] <- colnames(d2019)[i]
  output[i-6,2] <- round(cor(x,y, method = "spearman", use = "complete.obs"),2)
  spearman_p <- cor.test(x,y, method = "spearman", use = "complete.obs", exact = F)
  output[i-6,3] <- round(spearman_p$p.value,2)
  
}

output[,4] <- p.adjust(output[,3], method = "holm")

output <- data.frame(output)
colnames(output) <- c("Driver","Spearmans_rho","Spearmans_p","adj_Spearmans_p")

###max biomass
#read in FP data and drivers
fp <- read_csv("./2_Data_analysis/FP_megamatrix.csv")
colnames(fp)
fp <- fp[,c(1:6,8:10,18:23,25,28,35,27)] %>%
  mutate(WRT_day = log(WRT_day))

output <- matrix(NA, nrow = length(c(7:18)), ncol = 4)

for (i in 7:18){
  
  x <- unlist(fp[,i])
  y <- fp$Max_biomass_ugL
  
  output[i-6,1] <- colnames(fp)[i]
  output[i-6,2] <- round(cor(x,y, method = "spearman", use = "complete.obs"),2)
  spearman_p <- cor.test(x,y, method = "spearman", use = "complete.obs", exact = F)
  output[i-6,3] <- spearman_p$p.value
  
}

output[,4] <- p.adjust(output[,3], method = "holm")
output[,2:4] <- as.numeric(output[,2:4])
good_vars <- subset(output, unlist(output[,4]<0.05))
EM_vars <- good_vars[,1]

for (i in 1:length(EM_vars)){
  
  myplot <- ggplot(data = fp, aes_string(x = EM_vars[i], y = "Max_biomass_ugL"))+
    geom_point(size = 2)+
    theme_classic()+
    ggtitle(good_vars[i,2])
  print(myplot)
  
  
}

