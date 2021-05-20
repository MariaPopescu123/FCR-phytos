#3B_EM_manipulation_data_viz
#Author: Mary Lofton
#Date: 17DEC20
options(scipen=999)
pacman::p_load(tidyverse, lubridate, kSamples)
rm(list=ls())

#get data
my.fp.data <- read_csv("./2_Data_analysis/FP_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_DO_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") 
colnames(my.fp.data)
my.fp.data <- my.fp.data[,c(2,1,3:48)]

my.cs.data <- read_csv("./2_Data_analysis/CS_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_DO_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") %>%
  select(Date,shannon:BV_TOTAL) %>%
  filter(!is.na(BV_TOTAL)) 

mydata <- left_join(my.cs.data,my.fp.data, by = "Date") %>%
  mutate(Year = as.factor(Year), MonthDay = format(Date, format="%m-%d"),
         Peak_depth_m = ifelse(Peak_depth_m > 9, NA, Peak_depth_m)) 

colnames(mydata)
vars <- mydata[,c(66,51,50,54,52,53,57,60,59,61,43:48,38,40,42,37,39,41,30:32,34:36,62,63,65,6:14,15:23,3,2,5,4)]

final <- matrix(NA, nrow = length(c(66,51,50,54,52,53,57,60,59,61,43:48,38,40,42,37,39,41,30:32,34:36,62,63,65,6:14,15:23,3,2,5,4)), ncol = 15)

for (i in 1:53){

my.var <- unlist(vars[i])
final[i,1] <- colnames(vars)[i]

em <- my.var[1:35]
no_em <- my.var[36:67]

my.ad <- ad.test(em, no_em)

final[i,6] <- my.ad$ad[[5]]

final[i,2] <- mean(em, na.rm = TRUE)
final[i,3] <- sd(em, na.rm = TRUE)

final[i,4] <- mean(no_em, na.rm = TRUE)
final[i,5] <- sd(no_em, na.rm = TRUE)


d2016 <- my.var[1:20]
d2017 <- my.var[21:35]
d2018 <- my.var[36:50]
d2019 <- my.var[51:67]

final[i,8] <- mean(d2016, na.rm = TRUE)
final[i,9] <- sd(d2016, na.rm = TRUE)

final[i,10] <- mean(d2017, na.rm = TRUE)
final[i,11] <- sd(d2017, na.rm = TRUE)

final[i,12] <- mean(d2018, na.rm = TRUE)
final[i,13] <- sd(d2018, na.rm = TRUE)

final[i,14] <- mean(d2019, na.rm = TRUE)
final[i,15] <- sd(d2019, na.rm = TRUE)

}

final[,7] <- p.adjust(final[,6], method = "holm")
final <- data.frame(final) %>%
  mutate(X2 = as.numeric(X2),
         X3 = as.numeric(X3),
         X4 = as.numeric(X4),
         X5 = as.numeric(X5),
         X6 = as.numeric(X6),
         X7 = as.numeric(X7),
         X8 = as.numeric(X8),
         X9 = as.numeric(X9),
         X10 = as.numeric(X10),
         X11 = as.numeric(X11),
         X12 = as.numeric(X12),
         X13 = as.numeric(X13),
         X14 = as.numeric(X14),
         X15 = as.numeric(X15))
colnames(final) <- c("Driver","mean_EM","SD_EM","mean_no_EM","SD_no_EM","AD_p","adj_AD_p","mean_2016","SD_2016","mean_2017","SD_2017","mean_2018","SD_2018","mean_2019","SD_2019")
write.csv(final, "./2_Data_analysis/Anderson_Darling_results.csv",row.names = FALSE)
final <- final %>%
  filter(adj_AD_p <= 0.05)





###########prep for pairwise comparisons
sig.vars <- mydata[,colnames(mydata) %in% final$Driver]
final.sig <- matrix(NA, nrow = 23, ncol = 21)

for (i in 1:23){
  
  my.var <- unlist(sig.vars[i])
  final.sig[i,1] <- colnames(sig.vars)[i]
  
d2016 <- my.var[1:20]
d2017 <- my.var[21:35]
d2018 <- my.var[36:50]
d2019 <- my.var[51:67]

my.ad.1 <- ad.test(d2016,d2017)
my.ad.2 <- ad.test(d2016,d2018)
my.ad.3 <- ad.test(d2016,d2019)
my.ad.4 <- ad.test(d2017,d2018)
my.ad.5 <- ad.test(d2017,d2019)
my.ad.6 <- ad.test(d2018,d2019)

final.sig[i,2] <- my.ad.1$ad[[5]]
final.sig[i,3] <- my.ad.2$ad[[5]]
final.sig[i,4] <- my.ad.3$ad[[5]]
final.sig[i,5] <- my.ad.4$ad[[5]]
final.sig[i,6] <- my.ad.5$ad[[5]]
final.sig[i,7] <- my.ad.6$ad[[5]]

final.sig[i,8] <- mean(d2016, na.rm = TRUE)
final.sig[i,9] <- sd(d2016, na.rm = TRUE)

final.sig[i,10] <- mean(d2017, na.rm = TRUE)
final.sig[i,11] <- sd(d2017, na.rm = TRUE)

final.sig[i,12] <- mean(d2018, na.rm = TRUE)
final.sig[i,13] <- sd(d2018, na.rm = TRUE)

final.sig[i,14] <- mean(d2019, na.rm = TRUE)
final.sig[i,15] <- sd(d2019, na.rm = TRUE)
}

final.sig[,16] <- p.adjust(final.sig[,2], method = "holm")
final.sig[,17] <- p.adjust(final.sig[,3], method = "holm")
final.sig[,18] <- p.adjust(final.sig[,4], method = "holm")
final.sig[,19] <- p.adjust(final.sig[,5], method = "holm")
final.sig[,20] <- p.adjust(final.sig[,6], method = "holm")
final.sig[,21] <- p.adjust(final.sig[,7], method = "holm")

final.sig <- data.frame(final.sig) %>%
  mutate_at(c("X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20","X21"),as.numeric)
colnames(final.sig) <- c("Driver",
                     "AD_p_16_17",
                     "AD_p_16_18",
                     "AD_p_16_19",
                     "AD_p_17_18",
                     "AD_p_17_19",
                     "AD_p_18_19",
                     "mean_2016",
                     "sd_2016",
                     "mean_2017",
                     "sd_2017",
                     "mean_2018",
                     "sd_2018",
                     "mean_2019",
                     "sd_2019",
                     "adj_AD_p_16_17",
                     "adj_AD_p_16_18",
                     "adj_AD_p_16_19",
                     "adj_AD_p_17_18",
                     "adj_AD_p_17_19",
                     "adj_AD_p_18_19")
final.sig <- final.sig %>%
  mutate(EM_effect = ifelse(adj_AD_p_16_17 > 0.05 & adj_AD_p_18_19 > 0.05 & adj_AD_p_16_18 <= 0.05 & adj_AD_p_16_19 <= 0.05 & adj_AD_p_17_18 <= 0.05 & adj_AD_p_17_19 <= 0.05,"yes",
                            ifelse(adj_AD_p_16_17 > 0.05 & adj_AD_p_18_19 > 0.05 & adj_AD_p_16_18 < adj_AD_p_16_17 & adj_AD_p_16_19 < adj_AD_p_16_17 & adj_AD_p_17_18 < adj_AD_p_16_17 & adj_AD_p_17_19 < adj_AD_p_16_17 & adj_AD_p_16_18 < adj_AD_p_18_19 & adj_AD_p_16_19 < adj_AD_p_18_19 & adj_AD_p_17_18 < adj_AD_p_18_19 & adj_AD_p_17_19 < adj_AD_p_18_19,"maybe","no")))
write.csv(final.sig, "./2_Data_analysis/pairwise_Anderson_Darling_results.csv",row.names = FALSE)
final.sig <- final.sig %>%
  filter(EM_effect %in% c("yes","maybe"))




