#Table S4
#Author: Mary Lofton
#Date 06MAY21

####SET-UP####
#load packages
pacman::p_load(tidyverse, lubridate, vegan, Hmisc, picante, indicspecies, vegan3d)
rm(list=ls())

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
phytos <- read_csv("./00_Data_files/EDI_phytos/phytoplankton.csv") %>%
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

# any non-detects filled with 0
phytos2[phytos2 > 0]<- 1

phytos3 <- data.frame(colnames(phytos2),colSums(phytos2)) %>%
  arrange(desc(colSums.phytos2.))

write.csv(phytos3, "./2_Data_analysis/TableS4.csv",row.names = FALSE)

