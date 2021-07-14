#Table S5
#Author: Mary Lofton
#Date: 14JUL21

####run month-to-month ANOSIMs for 2016####
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

#run ANOSIMs
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
  select(-Year, -EM1, -EM2, -EM3, -perc_light_grab)

#select fp data for 2019
fp2019 <- fp %>% filter(year(Date) == 2019) 

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


####make Table S5####
final.2016$Year <- "2016"
final.2017$Year <- "2017"
final.2018$Year <- "2018"
final.2019$Year <- "2019"
final <- bind_rows(final.2016,final.2017,final.2018,final.2019) %>%
  mutate(Month_A = ifelse(Month_A == 5,"May",
                          ifelse(Month_A == 6, "June",
                                 ifelse(Month_A == 7,"July",
                                        ifelse(Month_A == 8,"August",
                                               ifelse(Month_A == 9,"September",""))))),
         Month_B = ifelse(Month_B == 5,"May",
                          ifelse(Month_B == 6, "June",
                                 ifelse(Month_B == 7,"July",
                                        ifelse(Month_B == 8,"August",
                                               ifelse(Month_B == 9,"September",""))))),
         comparison = paste(Month_A,Month_B, sep = "-"),
         R = round(R,2),
         p = round(p,3)) %>%
  select(Year, comparison, R, p)
write.csv(final,"./2_Data_analysis/intra-annual_pairwise_month_ANOSIM.csv",row.names = FALSE)
