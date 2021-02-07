#Title: Quebec lakes JJS dataset NMDS CORRELATION PLOTS W/ ENVIRONMENTAL VARIABLES
#Author: Mary Lofton
#Date last updated: 17FEB18

setwd("C:/Users/Mary Lofton/Dropbox/Quebec_lakes/NMDS")
library(tidyr)
library(vegan)
library(Hmisc)
library(picante)
library("ggpubr")

#read in environmental data
lake_data <- read.csv("physicochemical_summary.csv")

####READ IN NMDS FROM JUNE AND CREATE CORR PLOTS
#craft site abbreviations that make sense
#June = Ju, July = Jy, September = Se
#Baldwin = Ba
#Bowker = Bo
#Brome = Br
#Brompton = Bt
#DArgent = Da
#DesMonts = Dm
#Fitch = Fi
#Fraser = Fr
#Lovering = Lo
#Lyster = Ly
#Orford = Or
#Parker = Pa
#Simoneau = Si
#StGeorges = St
#Stukely = Su
#Tomcod = To
#TroisLacs = Tr
#Waterloo = Wa
lake_abbrevs <- c("Ba","Bo","Br","Bt","Da","Dm","Fi","Fr",'Lo','Ly','Or','Pa','Si',
                  'St','Su','To','Tr','Wa')
month_abbrevs <- c("Ju",'Jy','Se')

#create column of site abbrevs

abbrevs <- c(paste0(lake_abbrevs,'_',month_abbrevs[1]),
             paste0(lake_abbrevs,'_',month_abbrevs[2]),
             paste0(lake_abbrevs,'_',month_abbrevs[3]))
abbrevs <- sort(abbrevs)

#select the physicochemical variables you would like to include
nmds_data <- data.frame(abbrevs,lake_data[,c(4:7,9,10,12:16,18:19,30)])

#read in NMDS for each month
X <- readRDS("./June_NMDS_results.rds")
X1 <- readRDS("./July_NMDS_results.rds")
X2 <- readRDS("./Sept_NMDS_results.rds")

##JUNE (creates corr. plots for each NMDS axis and each environmental variable)

physics_June <- subset(nmds_data,lake_data$Month == "June")
mydata <- data.frame(X$points,physics_June)

setwd("C:/Users/Mary Lofton/Dropbox/Quebec_lakes/NMDS/June_NMDS_corr_plots")
for (i in 1:3){
for (j in 5:18){
  
nm <- c(colnames(mydata)[1:18])

filename = paste0(nm[i],"_", nm[j], "_June.png")

png(filename = filename,width = 9,height = 6,units = "in",res = 300)

print(ggscatter(mydata, x = nm[i], y = nm[j], 
          add = "reg.line", conf.int = TRUE, label = lake_abbrevs,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = nm[i], ylab = nm[j]))

dev.off()

}
}



##JULY

physics_July <- subset(nmds_data,lake_data$Month == "July")
mydata <- data.frame(X1$points,physics_July)

setwd("C:/Users/Mary Lofton/Dropbox/Quebec_lakes/NMDS/July_NMDS_corr_plots")
for (i in 1:3){
  for (j in 5:18){
    
    nm <- c(colnames(mydata)[1:18])
    
    filename = paste0(nm[i],"_", nm[j], "_July.png")
    
    png(filename = filename,width = 9,height = 6,units = "in",res = 300)
    
    print(ggscatter(mydata, x = nm[i], y = nm[j], 
                    add = "reg.line", conf.int = TRUE, label = lake_abbrevs,
                    cor.coef = TRUE, cor.method = "pearson",
                    xlab = nm[i], ylab = nm[j]))
    
    dev.off()
    
  }
}



##SEPT

physics_Sept <- subset(nmds_data,lake_data$Month == "September")
mydata <- data.frame(X2$points,physics_Sept)

setwd("C:/Users/Mary Lofton/Dropbox/Quebec_lakes/NMDS/Sept_NMDS_corr_plots")
for (i in 1:3){
  for (j in 5:18){
    
    nm <- c(colnames(mydata)[1:18])
    
    filename = paste0(nm[i],"_", nm[j], "_Sept.png")
    
    png(filename = filename,width = 9,height = 6,units = "in",res = 300)
    
    print(ggscatter(mydata, x = nm[i], y = nm[j], 
                    add = "reg.line", conf.int = TRUE, label = lake_abbrevs,
                    cor.coef = TRUE, cor.method = "pearson",
                    xlab = nm[i], ylab = nm[j]))
    
    dev.off()
    
  }
}



