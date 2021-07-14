#Table S7
#Author: Mary Lofton
#Date: 14JUL21

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
env <- my.cs.data[,c(2,5:6,11:13,27,29:30,38,39:40,42)]
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

####run ANOSIMs on all year-month pairings####
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

m.y <- c(1:20)
m.y.c <- c("May 2016","June 2016","July 2016","August 2016","September 2016",
           "May 2017","June 2017","July 2017","August 2017","September 2017",
           "May 2018","June 2018","July 2018","August 2018","September 2018",
           "May 2019","June 2019","July 2019","August 2019","September 2019")
m.y.key <- data.frame(m.y,m.y.c)
colnames(m.y.key) <- c("MonthYear_A","MonthYear_Char_A")

m.y.table <- left_join(final.m.y,m.y.key,by = "MonthYear_A")

m.y.key <- data.frame(m.y,m.y.c)
colnames(m.y.key) <- c("MonthYear_B","MonthYear_Char_B")

m.y.table2 <- left_join(m.y.table,m.y.key,by = "MonthYear_B") %>%
  mutate(comparison = paste(MonthYear_Char_A,MonthYear_Char_B, sep = "-"),
         R = round(R,2),
         p = round(p,3)) %>%
  select(comparison, R, p)
write.csv(m.y.table2,"./2_Data_analysis/Table_S7.csv",row.names = FALSE)
