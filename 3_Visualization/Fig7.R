#Fig. 4
#Author: Mary Lofton
#Date: 23FEB21

####SET-UP####
#load packages
pacman::p_load(tidyverse, lubridate, vegan, Hmisc, picante, indicspecies, vegan3d, cowplot)
rm(list=ls())

#set up vectors of different divisions
cyano <- c("Pseudanabaena","Aphanocapsa","Dolichospermum","Synechococcus",
           "Microcystis","unicell Cyanobacterium","Woronichinia","Aphanocapsa",
           "Merismopedia","Snowella","Dactylococcopsis","Oscillatoria" )
chloro <- c("Chlorophyte sp. 1","Ankistrodesmus" ,"Chlorophyte sp. 2",
            "Oocystis","Chlorophyte sp. 6","Selenastrum","Chlorophyte sp. 4",
            "Dictyosphaerium" ,"Schroederia","unicell Eukaryote","Tetraselmis",
            "Dysmorphococcus","Chlamydomonas","Nephroselmis","Elakatothrix",
            "Pandorina","Eudorina","Oonephris","Monomastix","Platydorina",
            "Pediastrum","Micractinium","Gloeocystis","Chlorophyte sp. 3",
            "Botryococcus","Kirchneriella","Quadrigula","Actinastrum" )
baci <- c("Nitzchia","Synedra","Asterionella","Cyclotella" ,"Fragilaria",
          "Ceratoneis","Tabellaria","Navicula" )
chryso <- c("Synura","Dinobryon","Oochromonas")
dino <- c("Gymnodinium","dinoflagellate cyst","Peridinium","naked dino","Parvodinium","Gloeodinium")
desmid <- c("Spondylosium","Staurastrum","Staurodesmus","Closterium",
            "Desmid spp. 1","Bambusina","Actinotaenium","Cosmarium" )
crypto <- c("Cryptomonas","Rhodomonas" )
eugleno <- c("Trachelomonas","Euglena","Euglena?" )
raphid <- c("Gonyostomum")


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

# prepare phyto data for NMDS input
phytos2 <- phytos1 %>%
  select(-Date,-EM1,-EM2,-EM3,-Month,-Year,-MonthDay) 

phytos3 <- as.matrix(phytos2)


####RUN THE NMDS FOR ALL YEARS####

# #calculate Hellinger distance - FOR COMMUNITY DATA
# WHY WOULD I DO THIS??
# hellinger_vars <- decostand(phytos[,-c(1:3)], method = "hellinger")

#based on results of scree plot, the best choice for k is 3

#run NMDS for all months
set.seed(1)
Q <- metaMDS(phytos3, distance='bray', k=3, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
plot(Q, display=c('sites'), choices=c(1,2), type='p')
plot(Q, display=c('sites'), choices=c(1,3), type='p')
# Q$species
Q$points

#subset important species
Y <- Q$species
colnames(phytos3)
Y.sel <- Y[c(8,14,20,22,24,38,42,49,51,53),]
Y.no <- Y[-c(8,14,20,22,24,38,42,49,51,53),]

#subset divisions
cy <- subset(Y.no, row.names(Y.no) %in% cyano)
chl <- subset(Y.no, row.names(Y.no) %in% chloro)
ba <- subset(Y.no, row.names(Y.no) %in% baci)
chr <- subset(Y.no, row.names(Y.no) %in% chryso)
di <- subset(Y.no, row.names(Y.no) %in% dino)
de <- subset(Y.no, row.names(Y.no) %in% desmid)
cr <- subset(Y.no, row.names(Y.no) %in% crypto)
eu <- subset(Y.no, row.names(Y.no) %in% eugleno)
ra <- subset(Y.no, row.names(Y.no) %in% raphid)

cy.sel <- subset(Y.sel, row.names(Y.sel) %in% cyano)
chl.sel <- subset(Y.sel, row.names(Y.sel) %in% chloro)
ba.sel <- subset(Y.sel, row.names(Y.sel) %in% baci)
chr.sel <- subset(Y.sel, row.names(Y.sel) %in% chryso)
di.sel <- subset(Y.sel, row.names(Y.sel) %in% dino)
de.sel <- subset(Y.sel, row.names(Y.sel) %in% desmid)
cr.sel <- subset(Y.sel, row.names(Y.sel) %in% crypto)
eu.sel <- subset(Y.sel, row.names(Y.sel) %in% eugleno)
ra.sel <- subset(Y.sel, row.names(Y.sel) %in% raphid)

#make plot
tiff(filename = "./3_Visualization/Fig5.tif",width = 5,height = 7.5,units = "in",res = 300)
par(mfrow = c(3,2),mar = c(3,3,1.5,0),mgp = c(1.5,0.5,0))

plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
title("A. all years", adj = 0.02, line = 0.5)
ordisurf(Q, my.cs.data$Peak_depth_m, display = c("sites", "species"),choices = c(1,2), type = "n",main = "",
         levels = c(3.6), add = TRUE, penalty = 1.4, method = "ML", col = "gray39")
ordisurf(Q, my.cs.data$Peak_depth_m, display = c("sites", "species"),choices = c(1,2), type = "n",main = "",
         levels = c(3.1), add = TRUE, penalty = 1.4, method = "ML", col = "gray79")
ordisurf(Q, my.cs.data$Peak_depth_m, display = c("sites", "species"),choices = c(1,2), type = "n",main = "",
         levels = c(4.9), add = TRUE, penalty = 1.4, method = "ML", col = "gray23")
points(cy[,1], cy[,2],  col="cadetblue3", pch = "+")
points(chl[,1], chl[,2],  col="darkgreen", pch = "+")
points(ba[,1], ba[,2],  col="burlywood3", pch = "+")
points(chr[,1], chr[,2],  col="darkgoldenrod2", pch = "+")
points(di[,1], di[,2],  col="brown1", pch = "+")
points(de[,1], de[,2],  col="chartreuse3", pch = "+")
points(cr[,1], cr[,2],  col="chocolate1", pch = "+")
points(eu[,1], eu[,2],  col="gray48", pch = "+")
points(ra[,1], ra[,2],  col="gray8", pch = "+")
text(cy.sel[,1], cy.sel[,2], substr(rownames(cy.sel),1,3), col="cadetblue3", cex = 0.8)
text(chl.sel[,1], chl.sel[,2], substr(rownames(chl.sel),1,3), col="darkgreen", cex = 0.8)
text(ba.sel[,1], ba.sel[,2], substr(rownames(ba.sel),1,3), col="burlywood3", cex = 0.8)
text(di.sel[,1], di.sel[,2], substr(rownames(di.sel),1,3), col="brown1", cex = 0.8)
text(de.sel[,1], de.sel[,2], substr(rownames(de.sel),1,3), col="chartreuse3", cex = 0.8)
text(cr.sel[,1], cr.sel[,2], substr(rownames(cr.sel),1,3), col="chocolate1", cex = 0.8)
text(eu.sel[,1], eu.sel[,2], substr(rownames(eu.sel),1,3), col="gray48", cex = 0.8)
legend("topright",legend = c("k = 0.11"),bty = "n")
legend("bottomleft",legend = c("Cyano","Chloro","Baci","Chryso","Dino","Desmid","Crypto","Eugleno","Raphid"),
       pch = 15, pt.cex = 2,col = c("cadetblue3","darkgreen","burlywood3","darkgoldenrod2","brown1","chartreuse3","chocolate1","gray48","gray8"),
       bty = "n")

plot(Q, display=c('sites','species'),choices=c(1,3), type='n')
title("B. all years", adj = 0.02, line = 0.5)
ordisurf(Q, my.cs.data$Peak_depth_m, display = c("sites", "species"),choices = c(1,3), type = "n",main = "",
         levels = c(3.6), add = TRUE, penalty = 1.4, method = "ML", col = "gray39")
ordisurf(Q, my.cs.data$Peak_depth_m, display = c("sites", "species"),choices = c(1,3), type = "n",main = "",
         levels = c(3.1), add = TRUE, penalty = 1.4, method = "ML", col = "gray79")
ordisurf(Q, my.cs.data$Peak_depth_m, display = c("sites", "species"),choices = c(1,3), type = "n",main = "",
         levels = c(4.9), add = TRUE, penalty = 1.4, method = "ML", col = "gray23")
points(cy[,1], cy[,3],  col="cadetblue3", pch = "+")
points(chl[,1], chl[,3],  col="darkgreen", pch = "+")
points(chr[,1], chr[,3],  col="darkgoldenrod2", pch = "+")
points(di[,1], di[,3],  col="brown1", pch = "+")
points(de[,1], de[,3],  col="chartreuse3", pch = "+")
points(cr[,1], cr[,3],  col="chocolate1", pch = "+")
points(eu[,1], eu[,3],  col="gray48", pch = "+")
points(ra[,1], ra[,3],  col="gray8", pch = "+")
text(cy.sel[,1], cy.sel[,3], substr(rownames(cy.sel),1,3), col="cadetblue3", cex = 0.8)
text(chl.sel[,1], chl.sel[,3], substr(rownames(chl.sel),1,3), col="darkgreen", cex = 0.8)
text(ba.sel[,1], ba.sel[,2], substr(rownames(ba.sel),1,3), col="burlywood3", cex = 0.8)
text(di.sel[,1], di.sel[,3], substr(rownames(di.sel),1,3), col="brown1", cex = 0.8)
text(de.sel[,1], de.sel[,3], substr(rownames(de.sel),1,3), col="chartreuse3", cex = 0.8)
text(cr.sel[,1], cr.sel[,3], substr(rownames(cr.sel),1,3), col="chocolate1", cex = 0.8)
text(eu.sel[,1], eu.sel[,3], substr(rownames(eu.sel),1,3), col="gray48", cex = 0.8)
legend("topright",legend = c("k = 0.11"),bty = "n")


#### NMDS FOR 2016 ONLY ####

#select 2016 data
p2016.data <- phytos1 %>%
  filter(year(Date) == 2016)

#prep data for NMDS input
p2016 <- p2016.data %>%
  select(-Year,-Date,-EM1,-EM2,-EM3,-Month,-MonthDay)
p2016 <- as.matrix(p2016)
pd2016 <- my.cs.data[1:20,]


#based on scree plot results, the best value for k is 2

#run NMDS for all months
set.seed(1)
Q <- metaMDS(p2016, distance='bray', k=2, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
# plot(Q, display=c('sites'), choices=c(1,2), type='p')

#select out environmental variables that are constant in 2016
env2016 <- env[1:20,] %>%
  select(EM3, Temp_C_grab, n2, Max_biomass_ugL, Month)
set.seed(1)
en_16 = envfit(Q, env2016, permutations = 9999, na.rm = TRUE)
row.names(en_16$vectors$arrows) <- c("DEEP2/N2/Mo","T","","Biom","")


#get indicator spp
Y <- Q$species
colnames(p2016)
Y.sel <- Y[c(3,8),]
Y.no <- Y[-c(3,8),]

#subset divisions
cy <- subset(Y.no, row.names(Y.no) %in% cyano)
chl <- subset(Y.no, row.names(Y.no) %in% chloro)
ba <- subset(Y.no, row.names(Y.no) %in% baci)
chr <- subset(Y.no, row.names(Y.no) %in% chryso)
di <- subset(Y.no, row.names(Y.no) %in% dino)
de <- subset(Y.no, row.names(Y.no) %in% desmid)
cr <- subset(Y.no, row.names(Y.no) %in% crypto)
eu <- subset(Y.no, row.names(Y.no) %in% eugleno)
ra <- subset(Y.no, row.names(Y.no) %in% raphid)

cy.sel <- subset(Y.sel, row.names(Y.sel) %in% cyano)
chl.sel <- subset(Y.sel, row.names(Y.sel) %in% chloro)
ba.sel <- subset(Y.sel, row.names(Y.sel) %in% baci)
chr.sel <- subset(Y.sel, row.names(Y.sel) %in% chryso)
di.sel <- subset(Y.sel, row.names(Y.sel) %in% dino)
de.sel <- subset(Y.sel, row.names(Y.sel) %in% desmid)
cr.sel <- subset(Y.sel, row.names(Y.sel) %in% crypto)
eu.sel <- subset(Y.sel, row.names(Y.sel) %in% eugleno)
ra.sel <- subset(Y.sel, row.names(Y.sel) %in% raphid)


plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
title("C. 2016", adj = 0.02, line = 0.5)
ordisurf(Q, pd2016$Peak_depth_m, display = c("sites", "species"),choices = c(1,2), type = "n",main = "",
         levels = c(3.6), add = TRUE, penalty = 1.4, method = "ML", col = "gray39")
ordisurf(Q, pd2016$Peak_depth_m, display = c("sites", "species"),choices = c(1,2), type = "n",main = "",
         levels = c(3.1), add = TRUE, penalty = 1.4, method = "ML", col = "gray79")
ordisurf(Q, pd2016$Peak_depth_m, display = c("sites", "species"),choices = c(1,2), type = "n",main = "",
         levels = c(4.9), add = TRUE, penalty = 1.4, method = "ML", col = "gray23")
points(cy[,1], cy[,2],  col="cadetblue3", pch = "+")
points(chl[,1], chl[,2],  col="darkgreen", pch = "+")
points(ba[,1], ba[,2],  col="burlywood3", pch = "+")
points(chr[,1], chr[,2],  col="darkgoldenrod2", pch = "+")
points(di[,1], di[,2],  col="brown1", pch = "+")
points(de[,1], de[,2],  col="chartreuse3", pch = "+")
points(cr[,1], cr[,2],  col="chocolate1", pch = "+")
points(eu[,1], eu[,2],  col="gray48", pch = "+")
points(ra[,1], ra[,2],  col="gray8", pch = "+")
plot(en_16, col = "black")
text(cy.sel[,1], cy.sel[,2], substr(rownames(cy.sel),1,3), col="cadetblue3", cex = 0.8)
text(chl.sel[,1], chl.sel[,2], substr(rownames(chl.sel),1,3), col="darkgreen", cex = 0.8)
legend("topright",legend = c("k = 0.13"),bty = "n")



#### 2017 ####

#select just 2017 data
p2017.data <- phytos1 %>%
  filter(year(Date) == 2017)

#prep data for NMDS input
p2017 <- p2017.data %>%
  select(-Year,-Date,-EM1,-EM2,-EM3,-Month,-MonthDay)
p2017 <- as.matrix(p2017)
pd2017 <- my.cs.data[21:35,]


#based on scree plot the best value for k is 2

#run NMDS for 2017
set.seed(1)
Q <- metaMDS(p2017, distance='bray', k=2, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
# plot(Q, display=c('sites'), choices=c(1,2), type='p')

#select out environmental variables that are constant in 2016
env2017 <- env[21:35,] %>%
  select(Grab_SRP_ugL, Grab_DOC_mgL, Month)
set.seed(1)
en_17 = envfit(Q, env2017, permutations = 9999, na.rm = TRUE)
row.names(en_17$vectors$arrows) <- c("SRP","DOC","Mo")


#get indicator spp
Y <- Q$species
colnames(p2017)
Y.sel <- Y[c(8,15,18,20,38),]
Y.no <- Y[-c(8,15,18,20,38),]

#subset divisions
cy <- subset(Y.no, row.names(Y.no) %in% cyano)
chl <- subset(Y.no, row.names(Y.no) %in% chloro)
ba <- subset(Y.no, row.names(Y.no) %in% baci)
chr <- subset(Y.no, row.names(Y.no) %in% chryso)
di <- subset(Y.no, row.names(Y.no) %in% dino)
de <- subset(Y.no, row.names(Y.no) %in% desmid)
cr <- subset(Y.no, row.names(Y.no) %in% crypto)
eu <- subset(Y.no, row.names(Y.no) %in% eugleno)
ra <- subset(Y.no, row.names(Y.no) %in% raphid)

cy.sel <- subset(Y.sel, row.names(Y.sel) %in% cyano)
chl.sel <- subset(Y.sel, row.names(Y.sel) %in% chloro)
ba.sel <- subset(Y.sel, row.names(Y.sel) %in% baci)
chr.sel <- subset(Y.sel, row.names(Y.sel) %in% chryso)
di.sel <- subset(Y.sel, row.names(Y.sel) %in% dino)
de.sel <- subset(Y.sel, row.names(Y.sel) %in% desmid)
cr.sel <- subset(Y.sel, row.names(Y.sel) %in% crypto)
eu.sel <- subset(Y.sel, row.names(Y.sel) %in% eugleno)
ra.sel <- subset(Y.sel, row.names(Y.sel) %in% raphid)


plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
title("D. 2017", adj = 0.02, line = 0.5)
ordisurf(Q, pd2017$Peak_depth_m, display = c("sites", "species"),choices = c(1,2), type = "n",main = "",
         levels = c(3.6), add = TRUE, penalty = 1.4, method = "ML", col = "gray39")
ordisurf(Q, pd2017$Peak_depth_m, display = c("sites", "species"),choices = c(1,2), type = "n",main = "",
         levels = c(3.1), add = TRUE, penalty = 1.4, method = "ML", col = "gray79")
ordisurf(Q, pd2017$Peak_depth_m, display = c("sites", "species"),choices = c(1,2), type = "n",main = "",
         levels = c(4.9), add = TRUE, penalty = 1.4, method = "ML", col = "gray23")
points(cy[,1], cy[,2],  col="cadetblue3", pch = "+")
points(chl[,1], chl[,2],  col="darkgreen", pch = "+")
points(ba[,1], ba[,2],  col="burlywood3", pch = "+")
points(chr[,1], chr[,2],  col="darkgoldenrod2", pch = "+")
points(di[,1], di[,2],  col="brown1", pch = "+")
points(de[,1], de[,2],  col="chartreuse3", pch = "+")
points(cr[,1], cr[,2],  col="chocolate1", pch = "+")
points(eu[,1], eu[,2],  col="gray48", pch = "+")
points(ra[,1], ra[,2],  col="gray8", pch = "+")
plot(en_17, col = "black")
text(cy.sel[,1], cy.sel[,2], substr(rownames(cy.sel),1,3), col="cadetblue3", cex = 0.8)
text(chl.sel[,1], chl.sel[,2], substr(rownames(chl.sel),1,3), col="darkgreen", cex = 0.8)
text(ba.sel[,1], ba.sel[,2], substr(rownames(ba.sel),1,3), col="burlywood3", cex = 0.8)
legend("topright",legend = c("k = 0.08"),bty = "n")
#### 2018 ####

#subset just 2018 data
p2018.data <- phytos1 %>%
  filter(year(Date) == 2018)

#prep data for NMDS input
p2018 <- p2018.data %>%
  select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay,-Year)
p2018 <- as.matrix(p2018)
pd2018 <- my.cs.data[36:50,]


#run NMDS for all months
set.seed(1)
Q <- metaMDS(p2018, distance='bray', k=2, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
# plot(Q, display=c('sites'), choices=c(1,2), type='p')

#select out environmental variables that are constant in 2016
env2018 <- env[36:50,] %>%
  select(Grab_DOC_mgL, Temp_C_grab, Month)
set.seed(1)
en_18 = envfit(Q, env2018, permutations = 9999, na.rm = TRUE)
row.names(en_18$vectors$arrows) <- c("DOC","T","Mo")

#get indicator spp
Y <- Q$species
colnames(p2018)
Y.sel <- Y[c(15,16),]
Y.no <- Y[-c(15),]

#subset divisions
cy <- subset(Y.no, row.names(Y.no) %in% cyano)
chl <- subset(Y.no, row.names(Y.no) %in% chloro)
ba <- subset(Y.no, row.names(Y.no) %in% baci)
chr <- subset(Y.no, row.names(Y.no) %in% chryso)
di <- subset(Y.no, row.names(Y.no) %in% dino)
de <- subset(Y.no, row.names(Y.no) %in% desmid)
cr <- subset(Y.no, row.names(Y.no) %in% crypto)
eu <- subset(Y.no, row.names(Y.no) %in% eugleno)
ra <- subset(Y.no, row.names(Y.no) %in% raphid)

cy.sel <- subset(Y.sel, row.names(Y.sel) %in% cyano)
chl.sel <- subset(Y.sel, row.names(Y.sel) %in% chloro)
ba.sel <- subset(Y.sel, row.names(Y.sel) %in% baci)
chr.sel <- subset(Y.sel, row.names(Y.sel) %in% chryso)
di.sel <- subset(Y.sel, row.names(Y.sel) %in% dino)
de.sel <- subset(Y.sel, row.names(Y.sel) %in% desmid)
cr.sel <- subset(Y.sel, row.names(Y.sel) %in% crypto)
eu.sel <- subset(Y.sel, row.names(Y.sel) %in% eugleno)
ra.sel <- subset(Y.sel, row.names(Y.sel) %in% raphid)


plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
title("E. 2018", adj = 0.02, line = 0.5)
ordisurf(Q, pd2018$Peak_depth_m, display = c("sites", "species"),choices = c(1,2), type = "n",main = "",
         levels = c(3.6), add = TRUE, penalty = 1.4, method = "ML", col = "gray39")
ordisurf(Q, pd2018$Peak_depth_m, display = c("sites", "species"),choices = c(1,2), type = "n",main = "",
         levels = c(3.1), add = TRUE, penalty = 1.4, method = "ML", col = "gray79")
ordisurf(Q, pd2018$Peak_depth_m, display = c("sites", "species"),choices = c(1,2), type = "n",main = "",
         levels = c(4.9), add = TRUE, penalty = 1.4, method = "ML", col = "gray23")
points(cy[,1], cy[,2],  col="cadetblue3", pch = "+")
points(chl[,1], chl[,2],  col="darkgreen", pch = "+")
points(ba[,1], ba[,2],  col="burlywood3", pch = "+")
points(chr[,1], chr[,2],  col="darkgoldenrod2", pch = "+")
points(di[,1], di[,2],  col="brown1", pch = "+")
points(de[,1], de[,2],  col="chartreuse3", pch = "+")
points(cr[,1], cr[,2],  col="chocolate1", pch = "+")
points(eu[,1], eu[,2],  col="gray48", pch = "+")
points(ra[,1], ra[,2],  col="gray8", pch = "+")
plot(en_18, col = "black")
text(ba.sel[,1], ba.sel[,2], substr(rownames(ba.sel),1,3), col="burlywood3", cex = 0.8)
legend("topright",legend = c("k = 0.12"),bty = "n")



#### 2019 ####

#subset 2019 data only
p2019.data <- phytos1 %>%
  filter(year(Date) == 2019)

#prep data for NMDS input
p2019 <- p2019.data %>%
  select(-Date,-EM1,-EM2,-EM3,-Month,-MonthDay,-Year)
p2019 <- as.matrix(p2019)
pd2019 <- my.cs.data[51:67,]



#run NMDS for all months
set.seed(1)
Q <- metaMDS(p2019, distance='bray', k=2, trymax=50, autotransform=FALSE, pc=FALSE, plot=FALSE)
Q$stress
# Q$species
# Q$points
# plot(Q, display=c('sites'), choices=c(1,2), type='p')

#select out environmental variables that are constant in 2016
env2019 <- env[51:67,] %>%
  select(Temp_C_grab, Peak_depth_m)
set.seed(1)
en_19 = envfit(Q, env2019, permutations = 9999, na.rm = TRUE)
row.names(en_19$vectors$arrows) <- c("T","PD")

#get indicator spp
Y <- Q$species
colnames(p2019)
Y.sel <- Y[c(4,20),]
Y.no <- Y[-c(4,20),]

#subset divisions
cy <- subset(Y.no, row.names(Y.no) %in% cyano)
chl <- subset(Y.no, row.names(Y.no) %in% chloro)
ba <- subset(Y.no, row.names(Y.no) %in% baci)
chr <- subset(Y.no, row.names(Y.no) %in% chryso)
di <- subset(Y.no, row.names(Y.no) %in% dino)
de <- subset(Y.no, row.names(Y.no) %in% desmid)
cr <- subset(Y.no, row.names(Y.no) %in% crypto)
eu <- subset(Y.no, row.names(Y.no) %in% eugleno)
ra <- subset(Y.no, row.names(Y.no) %in% raphid)

cy.sel <- subset(Y.sel, row.names(Y.sel) %in% cyano)
chl.sel <- subset(Y.sel, row.names(Y.sel) %in% chloro)
ba.sel <- subset(Y.sel, row.names(Y.sel) %in% baci)
chr.sel <- subset(Y.sel, row.names(Y.sel) %in% chryso)
di.sel <- subset(Y.sel, row.names(Y.sel) %in% dino)
de.sel <- subset(Y.sel, row.names(Y.sel) %in% desmid)
cr.sel <- subset(Y.sel, row.names(Y.sel) %in% crypto)
eu.sel <- subset(Y.sel, row.names(Y.sel) %in% eugleno)
ra.sel <- subset(Y.sel, row.names(Y.sel) %in% raphid)


plot(Q, display=c('sites','species'),choices=c(1,2), type='n')
title("F. 2019", adj = 0.02, line = 0.5)
ordisurf(Q, pd2019$Peak_depth_m, display = c("sites", "species"),choices = c(1,2), type = "n",main = "",
         levels = c(3.6), add = TRUE, penalty = 1.4, method = "ML", col = "gray39")
ordisurf(Q, pd2019$Peak_depth_m, display = c("sites", "species"),choices = c(1,2), type = "n",main = "",
         levels = c(3.1), add = TRUE, penalty = 1.4, method = "ML", col = "gray79")
ordisurf(Q, pd2019$Peak_depth_m, display = c("sites", "species"),choices = c(1,2), type = "n",main = "",
         levels = c(4.9), add = TRUE, penalty = 1.4, method = "ML", col = "gray23")
points(cy[,1], cy[,2],  col="cadetblue3", pch = "+")
points(chl[,1], chl[,2],  col="darkgreen", pch = "+")
points(ba[,1], ba[,2],  col="burlywood3", pch = "+")
points(chr[,1], chr[,2],  col="darkgoldenrod2", pch = "+")
points(di[,1], di[,2],  col="brown1", pch = "+")
points(de[,1], de[,2],  col="chartreuse3", pch = "+")
points(cr[,1], cr[,2],  col="chocolate1", pch = "+")
points(eu[,1], eu[,2],  col="gray48", pch = "+")
points(ra[,1], ra[,2],  col="gray8", pch = "+")
plot(en_19, col = "black")
text(cy.sel[,1], cy.sel[,2], substr(rownames(cy.sel),1,3), col="cadetblue3", cex = 0.8)
text(ba.sel[,1], ba.sel[,2], substr(rownames(ba.sel),1,3), col="burlywood3", cex = 0.8)
legend("topright",legend = c("k = 0.04"),bty = "n")


dev.off()



