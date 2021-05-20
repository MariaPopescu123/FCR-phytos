library(kSamples)
secchi <- read_csv("C:/Users/Mary Lofton/Downloads/Secchi_depth_2013-2020.csv") %>%
  filter(Reservoir == "FCR" & Site == 50 & date(DateTime) %in% sample_info$Date)
mean(secchi$Secchi_m)
sd(secchi$Secchi_m)

em <- subset(secchi, year(DateTime) %in% c(2016:2017))
no_em <- subset(secchi, year(DateTime) %in% c(2018:2019))

my.ad <- ad.test(em$Secchi_m, no_em$Secchi_m)

my.ad$ad[[5]]
mean(em$Secchi_m)
mean(no_em$Secchi_m)

d2016 <- subset(secchi, year(DateTime) ==2016)
d2017 <- subset(secchi, year(DateTime) ==2017)
d2018 <- subset(secchi, year(DateTime) ==2018)
d2019 <- subset(secchi, year(DateTime) ==2019)
mean(d2016$Secchi_m)
mean(d2017$Secchi_m)
mean(d2018$Secchi_m)
mean(d2019$Secchi_m)

my.ad.1 <- ad.test(d2016$Secchi_m,d2017$Secchi_m)
my.ad.2 <- ad.test(d2016$Secchi_m,d2018$Secchi_m)
my.ad.3 <- ad.test(d2016$Secchi_m,d2019$Secchi_m)
my.ad.4 <- ad.test(d2017$Secchi_m,d2018$Secchi_m)
my.ad.5 <- ad.test(d2017$Secchi_m,d2019$Secchi_m)
my.ad.6 <- ad.test(d2018$Secchi_m,d2019$Secchi_m)

my.ad.1$ad[[5]]
my.ad.2$ad[[5]]
my.ad.3$ad[[5]]
my.ad.4$ad[[5]]
my.ad.5$ad[[5]]
my.ad.6$ad[[5]]
