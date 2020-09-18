#1F_Biodiversity_indices
#Author: Mary Lofton
#Date: 17SEP20

#calculate biodiversity metrics on phytoplankton samples
#1. Shannon diversity
#2. richness
#3. Bray-Curtis (between sampling days)
#4. Jaccard (beween sampling days)

#load packages
#install.packages('pacman')
pacman::p_load(tidyverse, lubridate, rLakeAnalyzer, vegan)
rm(list=ls())

#data wrangling to get input matrix for vegan functions
phytos <- read_csv("./00_Data_files/EDI_phytos/phytoplankton.csv") %>%
  select(Date, Genus, BV_um3mL) %>%
  spread(key = Genus, value = BV_um3mL)

phytos[is.na(phytos)]<- 0

phytos1 <- phytos %>%
  select(-Date)

phytos2 <- as.matrix(phytos1)

#calculate Shannon diversity
shannon <- diversity(phytos2, index = "shannon", base = exp(1))
shannon1 <- data.frame(shannon)

BD <- data.frame(phytos$Date)
colnames(BD) <- "Date"
BD1 <- bind_cols(BD, shannon1)

#calculate richness
richness <- rowSums(phytos2 != 0)
richness1 <- data.frame(richness)

BD2 <- bind_cols(BD1, richness1)

#calculate Bray-Curtis distance
bc <- as.matrix(vegdist(phytos2, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = TRUE))

#ok now only calculate the change from week to week

final <- NULL

for (i in 1:99){
  final <- c(final, bc[i,i+1])
}

#elements that are skipping across a winter and should be NA
#23, 46, 66

final[c(23,25,46,66:70)] <- NA
final <- c(NA,final)

braycurtis <- data.frame(final)
colnames(braycurtis) <- "braycurtis"

BD3 <- bind_cols(BD2, braycurtis)

#calculate Jaccard distance
jac <- as.matrix(vegdist(phytos2, method="jaccard", binary=TRUE, diag=FALSE, upper=FALSE, na.rm = TRUE))

#ok now only calculate the change from week to week

final <- NULL

for (i in 1:99){
  final <- c(final, jac[i,i+1])
}

#elements that are skipping across a winter and should be NA
#23, 46, 66

final[c(23,25,46,66:70)] <- NA
final <- c(NA,final)

jaccard <- data.frame(final)
colnames(jaccard) <- "jaccard"

BD4 <- bind_cols(BD3, jaccard)

write.csv(BD4, "./00_Data_files/Biodiversity.csv",row.names = FALSE)

#visualization
bd <- read_csv("./00_Data_files/Biodiversity.csv") %>%
  mutate(Year = year(Date))

ggplot(data = bd, aes(x = Date, y = shannon))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = bd, aes(x = Date, y = richness))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = bd, aes(x = Date, y = braycurtis))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = bd, aes(x = Date, y = jaccard))+
  facet_wrap(vars(Year), scales = "free")+
  geom_line(size = 1)+
  theme_classic()
