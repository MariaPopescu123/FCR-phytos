#Title: EM2016 Phytoplankton biovolume calculations
#Author: Mary Lofton
#Date: 15FEB18

#NOTES 25FEB18:
#Still to be done:
#1. Need to step through code as it stands and check for errors - CHECK
#2. Need to run each day and write to file - CHECK
#3. Need to run each day to extract information for MBFGs and create csv - CHECK
#4. Need to assign qualitative traits for MBFGs - CHECK
#5. Need to calc MBFGs - CHECK
#6. Need to calc total biovolume of each MBFG for each sample day - CHECK!!!!
#7. Need to do basic visualizations of biovolume of each MBFG for FCR and BVR - CHECK!!!
#8. Need to run BACI t-test on biovolume of each MBFG for each sample day - CHECK!!!
#PACKAGES YOU NEED
library(readxl)
library(tidyverse)
library(rLakeAnalyzer)

#####################ENTER SAMPLE INFO HERE##########################
sample_volume = 1 #volume of sample counted in mL(remember to sum for total volume if you counted multiple slides)
chamber_type = 'SEDGWICK' #either SETTLING or SEDGWICK
digitized_before_ESA = FALSE #either true or false depending on date count sheets digitized
#####################################################################

##DATA IMPORT----
#read in your non-filamentous taxa data and convert to data frame
non_fila = read_excel("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/FCR_071116.xlsx", 
                      sheet='count calcs')

non_fila <- data.frame(non_fila)
non_fila[is.na(non_fila)] <- 0

#read in your filamentous taxa data and convert to data frame
fila = read_excel("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/FCR_071116.xlsx", 
                      sheet='colonies&filaments')
fila <- data.frame(fila)

#assign counted area based on slide type
if(chamber_type == "SETTLING"){
  volume_counted = 0.06375
}else{
  volume_counted = 0.0125
}

#create columns you will need for further calculations
non_fila$cells_per_mL <- NA
non_fila$cell_biovolume_um3 <- NA
non_fila$cell_SA_um2 <- NA

fila$cells_per_mL <- NA
fila$cell_biovolume_um3 <- NA
fila$cell_SA_um2 <- NA

##TAXON ASSIGNMENT TO SHAPES----
#assign biovolume equations based on taxon names
sphere <- c("Dictyosphaerium","Anabaena","Microcystis","Chromulina","Nephroselmis",
            "Carteria","Chlamydomonas","Eudorina","Golenkinia","pico round",
            "Synechococcus","Volvox","Chlorella","Volvulina","pico round flagellate")
cylinder <- c("Phormidium","pico oblong","filament","Pseudanabaena","Microspora")
prolate_spheroid <- c("Dinobryon","Mallomonas","Synura","Centritractus","Cryptomonas",
                      "Tetraselmis","Oocystis","Pandorina")
prism_on_triangle <- c("Goniochloris")
box_and_2_cylinders <- c("")
prism_on_parallelogram <- c("Nitzschia")
box <- c("Synedra","Cymatopleura","Asterionella","Bacillaria","Tabellaria")
cone_and_half_sphere <- c("Rhodomonas","Chroomonas")
elliptic_prism <- c("Phacus","Mesostigma")
ellipsoid <- c("Gymnodinium","Peridinium","Trachelomonas","Coccomyxa","Cosmarium","Spondylosium")
cone <- c("Pyramimonas")
cylinder_and_2_cones <- c("Ankistrodesmus")
two_cones <- c("Chlorogonium","Closterium")
two_half_ellipsoids <- c("")
sickle_shaped_prism <- c("Selenastrum")
cylinder_and_2_half_spheres <- c("Gomphonema")

#assemble list of taxa already assigned to shapes
names <- c(sphere, cylinder, prolate_spheroid, prism_on_triangle, box_and_2_cylinders, prism_on_parallelogram,
           box, cone_and_half_sphere, elliptic_prism, ellipsoid, cone, cylinder_and_2_cones, two_cones,
           two_half_ellipsoids, sickle_shaped_prism,cylinder_and_2_half_spheres)

#check to make sure there aren't any rogue taxa in data frame that aren't assigned to a shape
for (i in 1:length(non_fila[,1])){
if(non_fila$Taxon[i] %in% names){
  print("yay")
} else {
  print(non_fila$Taxon[i])
}
}

##CALCULATE AVG. LENGTH, WIDTH, AND CELLS PER ML----

if(digitized_before_ESA == TRUE){
  
  for (i in 1:length(non_fila[,1])){ #only use this for early digitized files (pre-03AUG17)
  if(is.na(non_fila$Count[i]) == FALSE){
    non_fila$cells_per_mL[i] <- non_fila$Count[i]*(sample_volume/volume_counted)
  }
    #filamentous taxa
    fila_temp <- matrix(data=NA, ncol=7, nrow=8)
    fila_temp <- data.frame(fila_temp)
    colnames(fila_temp) <- colnames(non_fila)
    fila_temp$Taxon <- c(rep("Microcystis",2),rep("Phormidium",2),rep("Anabaena",2),rep("Pseudanabaena",2))
    fila_temp$varname <- non_fila$varname[1:8]
}
}else{
  #for non-filamentous taxa
  for (i in 1:length(non_fila[,1])){
    #this part of the loop calculates average length and width of counted cells for each taxon
    d <- non_fila[i,c(3:12)]
    d <- data.frame(d)
    d <- d[!is.na(d)]
    if(length(d)>0){
      non_fila$Avg..BV[i] <- mean(d)*25
    }else{
      non_fila$Avg..BV[i]<- NA
    }
    
    if(is.na(non_fila$Count[i]) == FALSE){
      non_fila$cells_per_mL[i] <- non_fila$Count[i]*(sample_volume/volume_counted)
    }
  }
  #filamentous taxa
  fila_temp <- matrix(data=NA, ncol=17, nrow=8)
  fila_temp <- data.frame(fila_temp)
  colnames(fila_temp) <- colnames(non_fila)
  fila_temp$Taxon <- c(rep("Microcystis",2),rep("Phormidium",2),rep("Anabaena",2),rep("Pseudanabaena",2))
  fila_temp$varname <- non_fila$varname[1:8]
}


if("Microcystis" %in% colnames(fila))
{
  fila_temp$Avg..BV[c(1:2)] <- 0.1
  fila_temp$Count[1:2] <- c(sum(fila$X.cells, na.rm = TRUE),0)
} else {
  fila_temp$Avg..BV[c(1:2)] <- 0
  fila_temp$Count[c(1:2)] <- 0
}
if("Phormidium" %in% colnames(fila))
{
  fila_temp$Avg..BV[3] <- mean(fila$length, na.rm = TRUE)
  fila_temp$Avg..BV[c(4)] <- 0.1
  fila_temp$Count[3:4] <- c(sum(!is.na(fila$length)) ,0)
} else {
  fila_temp$Avg..BV[c(3:4)] <- 0
  fila_temp$Count[c(3:4)] <- 0
}
if("Anabaena" %in% colnames(fila))
{
  fila_temp$Avg..BV[5] <- mean(fila$length__1, na.rm = TRUE)
  fila_temp$Avg..BV[6] <- mean(fila$cell.size, na.rm = TRUE)
  fila_temp$Count[5:6] <- c(sum(!is.na(fila$length__1)) ,0)
} else {
  fila_temp$Avg..BV[c(5:6)] <- 0
  fila_temp$Count[c(5:6)] <- 0
}
if("Psuedanabaena" %in% colnames(fila))
{
  fila_temp$Avg..BV[7] <- mean(fila$length__2, na.rm = TRUE)
  fila_temp$Avg..BV[8] <- mean(fila$cell.size_1, na.rm = TRUE)
  fila_temp$Avg..BV[7:8] <- c(sum(!is.na(fila$length__2)) ,0)
} else {
  fila_temp$Avg..BV[c(7:8)] <- 0
  fila_temp$Count[c(7:8)] <- 0
}

if(digitized_before_ESA == TRUE){
fila_temp$Avg..BV <- fila_temp$Avg..BV*25
non_fila$Avg..BV <- non_fila$Avg..BV*25
}

for (i in 1:length(fila_temp[,1])){
if(is.na(fila_temp$Count[i]) == FALSE){
  fila_temp$cells_per_mL[i] <- fila_temp$Count[i]*(sample_volume/volume_counted)
}
}

#bind fila & non_fila taxa
non_fila <- rbind(non_fila,fila_temp)

##CELL BIOVOLUME----
#this loop calculates cell biovolume based on the taxon
for (i in 1:length(non_fila[,1])){
  if (non_fila$Taxon[i] %in% sphere == TRUE){
    non_fila$cell_biovolume_um3[i] <- (pi/6)*(non_fila$Avg..BV[i]^3) 
  }else if (non_fila$Taxon[i] %in% cylinder == TRUE){
    non_fila$cell_biovolume_um3[i] <- (pi/4)*(non_fila$Avg..BV[i+1]^2)*non_fila$Avg..BV[i]
  }else if (non_fila$Taxon[i] %in% prolate_spheroid == TRUE){
    non_fila$cell_biovolume_um3[i] <- (pi/6)*(non_fila$Avg..BV[i+1]^2)*non_fila$Avg..BV[i]
  }else if (non_fila$Taxon[i] %in% prism_on_triangle == TRUE){
    non_fila$cell_biovolume_um3[i] <- 0.5*((sqrt(3)/2)*non_fila$Avg..BV[i+1])*non_fila$Avg..BV[i]
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% box_and_2_cylinders == TRUE){
    non_fila$cell_biovolume_um3[i] <- NA
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% prism_on_parallelogram == TRUE){
    non_fila$cell_biovolume_um3[i] <- 0.5*non_fila$Avg..BV[i]*non_fila$Avg..BV[i+1]*(0.5*non_fila$Avg..BV[i+1])
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% box == TRUE){
    non_fila$cell_biovolume_um3[i] <- non_fila$Avg..BV[i]*(non_fila$Avg..BV[i+1]^2)
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% cone_and_half_sphere == TRUE){
    non_fila$cell_biovolume_um3[i] <- (pi/12)*(non_fila$Avg..BV[i+1]^2)*(0.5*non_fila$Avg..BV[i]+non_fila$Avg..BV[i+1])
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% elliptic_prism == TRUE){
    non_fila$cell_biovolume_um3[i] <- (pi/4)*non_fila$Avg..BV[i]*non_fila$Avg..BV[i+1]*(0.5*non_fila$Avg..BV[i+1])
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% ellipsoid == TRUE){
    non_fila$cell_biovolume_um3[i] <- (pi/6)*(non_fila$Avg..BV[i+1]^2)*non_fila$Avg..BV[i]
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% cone == TRUE){
    non_fila$cell_biovolume_um3[i] <- (pi/2)*non_fila$Avg..BV[i+1]*((non_fila$Avg..BV[i+1]/2)+((sqrt(3)/2))*non_fila$Avg..BV[i+1])
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% cylinder_and_2_cones == TRUE){
    non_fila$cell_biovolume_um3[i] <- (pi/4)*(non_fila$Avg..BV[i+1]^2)*((0.8*non_fila$Avg..BV[i])+((0.1*non_fila$Avg..BV[i])/2))
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% two_cones == TRUE){
    non_fila$cell_biovolume_um3[i] <- (pi/6)*(non_fila$Avg..BV[i+1]^2)*(0.5*non_fila$Avg..BV[i])
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% two_half_ellipsoids == TRUE){
    non_fila$cell_biovolume_um3[i] <- NA
  }else if (non_fila$Taxon[i] %in% sickle_shaped_prism == TRUE){
    non_fila$cell_biovolume_um3[i] <- (pi/4)*non_fila$Avg..BV[i+1]*((non_fila$Avg..BV[i]*2*non_fila$Avg..BV[i+1])-(non_fila$Avg..BV[i]*non_fila$Avg..BV[i+1]))
  } else if (non_fila$Taxon[i] %in% cylinder_and_2_half_spheres == TRUE){
    non_fila$cell_biovolume_um3[i] <- pi*(non_fila$Avg..BV[i+1]^2)*((non_fila$Avg..BV[i]/4)+(non_fila$Avg..BV[i+1]/6))
  }
}

##CELL SURFACE AREA----
#this for-loop calculates cell surface area based on the taxon
for (i in 1:length(non_fila[,1])){
  if (non_fila$Taxon[i] %in% sphere == TRUE){
    non_fila$cell_SA_um2[i] <- pi*(non_fila$Avg..BV[i]^2) 
  }else if (non_fila$Taxon[i] %in% cylinder == TRUE){
    non_fila$cell_SA_um2[i] <- pi*non_fila$Avg..BV[i+1]*((non_fila$Avg..BV[i+1]/2)+non_fila$Avg..BV[i])
  }else if (non_fila$Taxon[i] %in% prolate_spheroid == TRUE && non_fila$Avg..BV[i] == non_fila$Avg..BV[i+1]){
    non_fila$cell_SA_um2[i] <- pi*(non_fila$Avg..BV[i]^2)
  }else if (non_fila$Taxon[i] %in% prolate_spheroid == TRUE && non_fila$Avg..BV[i] > non_fila$Avg..BV[i+1]){
    non_fila$cell_SA_um2[i] <- ((pi*non_fila$Avg..BV[i+1])/2)*(non_fila$Avg..BV[i+1]+((non_fila$Avg..BV[i]^2)/sqrt((non_fila$Avg..BV[i]^2)-(non_fila$Avg..BV[i+1]^2))*asin(sqrt((non_fila$Avg..BV[i]^2)-(non_fila$Avg..BV[i+1]^2))/non_fila$Avg..BV[i])))
  }else if (non_fila$Taxon[i] %in% prism_on_triangle == TRUE){
    non_fila$cell_SA_um2[i] <- (non_fila$Avg..BV[i+1]*((sqrt(3)/2)*non_fila$Avg..BV[i+1]))+3*non_fila$Avg..BV[i+1]*non_fila$Avg..BV[i]
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% box_and_2_cylinders == TRUE){
    non_fila$cell_SA_um2[i] <- NA
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% prism_on_parallelogram == TRUE){
    non_fila$cell_SA_um2[i] <- (non_fila$Avg..BV[i]*non_fila$Avg..BV[i+1])+((sqrt((non_fila$Avg..BV[i]^2)+(non_fila$Avg..BV[i+1]^2))/4)*(0.5*non_fila$Avg..BV[i+1]))
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% box == TRUE){
    non_fila$cell_SA_um2[i] <- (4*non_fila$Avg..BV[i]*non_fila$Avg..BV[i+1])+(2*(non_fila$Avg..BV[i+1]^2))
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% cone_and_half_sphere == TRUE){
    non_fila$cell_SA_um2[i] <- (pi/2)*non_fila$Avg..BV[i+1]*((non_fila$Avg..BV[i]/sqrt(3))+non_fila$Avg..BV[i+1])
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% elliptic_prism == TRUE){
    non_fila$cell_SA_um2[i] <- (pi/2)*((non_fila$Avg..BV[i]*non_fila$Avg..BV[i+1])+((non_fila$Avg..BV[i]+non_fila$Avg..BV[i+1])*0.5*non_fila$Avg..BV[i+1]))
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% ellipsoid == TRUE && non_fila$Avg..BV[i] == non_fila$Avg..BV[i+1]){
      non_fila$cell_SA_um2[i] <- pi*(non_fila$Avg..BV[i]^2)
  }else if (non_fila$Taxon[i] %in% ellipsoid == TRUE && non_fila$Avg..BV[i] > non_fila$Avg..BV[i+1]){
      non_fila$cell_SA_um2[i] <- ((pi*non_fila$Avg..BV[i+1])/2)*(non_fila$Avg..BV[i+1]+((non_fila$Avg..BV[i]^2)/sqrt((non_fila$Avg..BV[i]^2)-(non_fila$Avg..BV[i+1]^2))*asin(sqrt((non_fila$Avg..BV[i]^2)-(non_fila$Avg..BV[i+1]^2))/non_fila$Avg..BV[i])))
    #used SA of prolate spheroid because ellipsoid not included in Hillenbrand et al.
  }else if (non_fila$Taxon[i] %in% cone == TRUE){
    non_fila$cell_SA_um2[i] <- (pi/2)*non_fila$Avg..BV[i+1]*((non_fila$Avg..BV[i+1]/2)+(non_fila$Avg..BV[i]/sqrt(3)))
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% cylinder_and_2_cones == TRUE){
    non_fila$cell_SA_um2[i] <- pi*non_fila$Avg..BV[i+1]*((0.8*non_fila$Avg..BV[i])+(non_fila$Avg..BV[i]/sqrt(3)))
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% two_cones == TRUE){
    non_fila$cell_SA_um2[i] <- pi*non_fila$Avg..BV[i+1]*(non_fila$Avg..BV[i]/sqrt(3))
    #required extra measurement that was estimated
  }else if (non_fila$Taxon[i] %in% two_half_ellipsoids == TRUE){
    non_fila$cell_SA_um2[i] <- NA
  }else if (non_fila$Taxon[i] %in% sickle_shaped_prism == TRUE){
    non_fila$cell_SA_um2[i] <- (pi/4)*((non_fila$Avg..BV[i]*2*non_fila$Avg..BV[i+1])-(non_fila$Avg..BV[i]*non_fila$Avg..BV[i+1])+((non_fila$Avg..BV[i]+2*non_fila$Avg..BV[i+1])*non_fila$Avg..BV[i+1])+((non_fila$Avg..BV[i]+non_fila$Avg..BV[i+1])*non_fila$Avg..BV[i+1]))
  }else if (non_fila$Taxon[i] %in% cylinder_and_2_half_spheres == TRUE){
    non_fila$cell_SA_um2[i] <- pi*non_fila$Avg..BV[i+1]*(non_fila$Avg..BV[i+1]+non_fila$Avg..BV[i])
  }
}

##DATA WRANGLING FOR EXPORT OF ONE SAMPLE DAY----

#select the parts of the data frame that you want to keep for further use
final_non_fila <- subset(non_fila,non_fila$varname == "length")

if(digitized_before_ESA == TRUE){
  final_non_fila <- final_non_fila[,c(1,3:7)]
}else{
  final_non_fila <- final_non_fila[,c(1,13:17)]
}

colnames(final_non_fila)[2] <- "MLD" 

#replace NAs with 0s (appropriate in this case)
final_non_fila[is.na(final_non_fila)] <- 0

#calculate total biovolume of each taxon in the sample
final_non_fila$biovolume_um3_per_mL <- final_non_fila$cells_per_mL*final_non_fila$cell_biovolume_um3

#calculate relative abundance of each taxon in the sample
final_non_fila$rel_abundance <- final_non_fila$biovolume_um3_per_mL/sum(final_non_fila$biovolume_um3_per_mL)

#calculate SA/V ratio of each taxon in the sample
final_non_fila$SA_to_V <- final_non_fila$cell_SA_um2/final_non_fila$cell_biovolume_um3

#remove taxa that were not actually present
final_non_fila <- subset(final_non_fila, final_non_fila$Count != 0)

setwd("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/FCR_processed")
write.csv(final_non_fila, "071116_processed.csv",row.names = FALSE)

##CREATE THE VARIABLES.CSV FILE FOR THE MBFG CALCULATION-----

##prior to running this section, you need to make sure 'names' is defined
setwd("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/FCR_processed")
sample_files <- list.files("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/FCR_processed")

mld <- matrix(data=NA, ncol=1, nrow=length(names))
mld[,1] <- names
mld <- data.frame(mld)
colnames(mld) <- c("Taxon")

vol <- matrix(data=NA, ncol=1, nrow=length(names))
vol[,1] <- names
vol <- data.frame(vol)
colnames(vol) <- c("Taxon")

sa_v <- matrix(data=NA, ncol=1, nrow=length(names))
sa_v[,1] <- names
sa_v <- data.frame(sa_v)
colnames(sa_v) <- c("Taxon")

for (i in 1:length(sample_files)){
  data <- read.csv(sample_files[i])
  s = sample_files[i]
  s1 = unlist(strsplit(s, split='_', fixed=TRUE))[1]
  colnames(data)[c(2,5,9)] <- s1

  mld <- merge(mld,data[,c(1:2)],by = "Taxon",all.x = TRUE)

  vol <- merge(vol,data[,c(1,5)],by = "Taxon",all.x = TRUE)

  sa_v <- merge(sa_v,data[,c(1,9)],by = "Taxon",all.x = TRUE)

}

mld$avg <- rowMeans(mld[,-1],na.rm = TRUE)
vol$avg <- rowMeans(vol[,-1],na.rm = TRUE)
sa_v$avg <- rowMeans(sa_v[,-1],na.rm = TRUE)


final <- data.frame(mld[,1], mld[,16],vol[,16],sa_v[,16])
colnames(final) <- c("Taxon","MLD","Volume","SV")
final <- final[,c(1,3,4,2)]
final <- final[-c(1:2,16,21,34),]
final$Flagella <- NA
final$Silice <- NA
final$Mucilage <- NA
final$Aerotopes <- NA
#ADD SOMETHING HERE TO ELIMINATE ROWS OF TAXA NOT FOUND
write.csv(final,"variables.csv",row.names = FALSE)

#POST-MBFGs - calc. biovolume of each MBFG for each sample day -----

##again, need to have names loaded to run this

setwd("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/MBFGs_FCR")

#read in mbfgs and list of summary files from each sample day
mbfg <- read.csv("./CLASIF_31AUG18.csv")
sample_files <- list.files("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/FCR_processed")
#sample_files <- sample_files[-15]

#create matrix to populate with biovolume
tot_bv <- matrix(data=NA, ncol=1, nrow=length(names))
tot_bv[,1] <- names
tot_bv <- data.frame(tot_bv)
colnames(tot_bv) <- c("Taxon")

#change working directory for for-loop
setwd("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/BVR_processed")

#create matrix of biovolume of each taxon for each sample day
for (i in 1:length(sample_files)){
  data <- read.csv(sample_files[i])
  s = sample_files[i]
  s1 = unlist(strsplit(s, split='_', fixed=TRUE))[1]
  colnames(data)[7] <- s1
  
  tot_bv <- merge(tot_bv,data[,c(1,7)],by = "Taxon",all.x = TRUE)
  
  
}

#combine biovolume matrix with mbfg matrix
colnames(mbfg)[9] <- "Taxon"
mbfg_bv <- merge(mbfg,tot_bv,by = "Taxon",all.x = TRUE)


###########################################
#IF YOU ARE WORKING WITH BVR, RUN THIS CODE!
##also there is some shitty coding squirreliness regarding number of columns through here in FCR vs. BVR
setwd("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/MBFGs_BVR")
mbfg_bv <- read.csv("./CLASIF_31AUG18.csv")

#create matrix to populate with biovolume of each mbfg
final_mbfg_bv <- matrix(data=NA, ncol=6, nrow=7)
dates <- colnames(mbfg_bv)[10:14]
columns <- c("MBFG",dates)
final_mbfg_bv <- data.frame(final_mbfg_bv)
colnames(final_mbfg_bv) <- columns

#subset each MBFG and sum across biovolume for each sample day

for (i in 1:7){
data <- subset(mbfg_bv, mbfg_bv$FUNC_GROUP == i)
data <- data[,-c(1:9)]
sums <- colSums(data,na.rm = TRUE)
row <- c(i,sums)
final_mbfg_bv[i,] <- row
}

#get total biovolume so can do relative abundance
bv <- colSums(final_mbfg_bv[,c(2:6)],na.rm = TRUE)
bv <- c(NA,bv)
final_mbfg_bv[8,] <- bv

#calculate relative abundance
rel_a <- data.frame(lapply(final_mbfg_bv, function(X) X/X[8]))
rel_a[,1] <- c(1:7,NA)

#transpose so data is tidy
final_mbfg_bv = setNames(data.frame(t(final_mbfg_bv[,-1])), final_mbfg_bv[,1])
colnames(final_mbfg_bv)[8] <- "TOTAL"

rel_a = setNames(data.frame(t(rel_a[,-1])), rel_a[,1])


#write data to file
setwd("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/MBFGs_BVR")
write.csv(final_mbfg_bv,"MBFG_biovolume_31AUG18.csv")
write.csv(rel_a,"MBFG_rel_abund_31AUG18.csv")

###crappy coding skill means you still need to go in and manually fix dates to load in as ts for rLakeAnalyzer :(


##PRELIMINARY DATA VISUALIZATION AND COLLATING CSV FOR BACI TEST

#set working directories and import data
#working with rel. abundance here because it's comparable across functional groups

#FCR
setwd("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/FCR_processed/MBFGs_FCR")
FCR <- load.ts("./MBFG_rel_abund.txt")

#BVR
setwd("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/MBFGs_BVR")
BVR <- load.ts("./MBFG_rel_abund.txt")


#basic time series plot for FCR
setwd("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/FCR_processed/MBFGs_FCR")
png("FCR_MBFG_rel_abund.png",width = 8,height = 6,units = "in",res = 500)
par(mgp = c(2.3,1,0))
plot(FCR$datetime,FCR$x1,
     ylim = c(0,1),
     xlim = c(FCR$datetime[1],FCR$datetime[14]+600000),
     col = 1,
     type = "l",
     lwd = 2,
     ylab = "MBFG relative abundance",
     xlab = "",
     cex.lab = 1.5,
     cex.axis = 1.3)
for (i in 3:8){
  points(FCR$datetime,FCR[,i],col = i-1,type = "l",lwd = 2)
}
legend("topright",lty = rep(1,7),lwd = rep(2,7), legend = c(1:7), bty = "n", col = c(1:7),cex = 1.2)
abline(v=as.numeric(FCR$datetime[3]), lwd=2, col='dimgray')
abline(v=as.numeric(FCR$datetime[10]), lwd=2, col='dimgray')
dev.off()

#basic times series plot for FCR - EM2 only
setwd("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/FCR_processed/MBFGs_FCR")
png("FCR_MBFG_rel_abund_EM2.png",width = 8,height = 6,units = "in",res = 500)
par(mgp = c(2.3,1,0))
plot(FCR$datetime[9:13],FCR$x1[9:13],
     ylim = c(0,1),
     xlim = c(FCR$datetime[9],FCR$datetime[13]+300000),
     col = 1,
     type = "l",
     lwd = 2,
     ylab = "MBFG relative abundance",
     xlab = "",
     cex.lab = 1.5,
     cex.axis = 1.3)
for (i in 3:8){
  points(FCR$datetime[9:13],FCR[9:13,i],col = i-1,type = "l",lwd = 2)
}
legend("topright",lty = rep(1,7),lwd = rep(2,7), legend = c(1:7), bty = "n", col = c(1:7),cex = 1.2)
abline(v=as.numeric(BVR$datetime[2]+(0.5*(BVR$datetime[3]-BVR$datetime[2]))), lwd=2, col='dimgray')
dev.off()

#basic times series plot for BVR
setwd("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/MBFGs_BVR")
png("BVR_MBFG_rel_abund.png",width = 8,height = 6,units = "in",res = 500)
par(mgp = c(2.3,1,0))
plot(BVR$datetime,BVR$x1,
     ylim = c(0,1),
     xlim = c(BVR$datetime[1],BVR$datetime[5]+300000),
     col = 1,
     type = "l",
     lwd = 2,
     ylab = "MBFG relative abundance",
     xlab = "",
     cex.lab = 1.5,
     cex.axis = 1.3)
for (i in 3:8){
  points(BVR$datetime,BVR[,i],col = i-1,type = "l",lwd = 2)
}
legend("topright",lty = rep(1,7),lwd = rep(2,7), legend = c(1:7), bty = "n", col = c(1:7),cex = 1.2)
abline(v=as.numeric(BVR$datetime[2]+(0.5*(BVR$datetime[3]-BVR$datetime[2]))), lwd=2, col='dimgray')
dev.off()

#collate data for input into BACI
BVR$Site <- "BVR"
FCR$Site <- "FCR"

FCR_EM2 <- FCR[9:13,]
FCR_EM2$Period <- c(rep("Before",2),rep("After",3))
FCR_EM2$SamplingTime <- c(1:5)
FCR_EM2$SiteClass <- "Impact"
FCR_EM2 <- FCR_EM2[,-9]
FCR_EM2 <- FCR_EM2[,c(11,12,10,1,2:9)]

BVR_EM2 <- BVR
BVR_EM2$Period <- c(rep("Before",2),rep("After",3))
BVR_EM2$SamplingTime <- c(1:5)
BVR_EM2$SiteClass <- "Control"
BVR_EM2 <- BVR_EM2[,-9]
BVR_EM2 <- BVR_EM2[,c(11,12,10,1,2:9)]

final <- rbind(FCR_EM2,BVR_EM2)

for (i in 1:7){
  new_folder <- paste("MBFG",i,sep = "_")
dir.create(paste("./EM2",new_folder,sep = "/"))
write.csv(final,paste("./EM2",new_folder,"BACI_MBFGs.csv",sep = "/"))
}


#same thing except for biovolume


#set working directories and import data

#FCR
setwd("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/FCR_processed/MBFGs_FCR")
FCR <- load.ts("./MBFG_biovolume.txt")
FCR[FCR == 0] <- 0.001

#BVR
setwd("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/MBFGs_BVR")
BVR <- load.ts("./MBFG_biovolume.txt")


#basic time series plot for FCR
setwd("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/FCR_processed/MBFGs_FCR")
png("FCR_MBFG_biovolume.png",width = 8,height = 5,units = "in",res = 500)
par(mgp = c(2.3,1,0))
plot(FCR$datetime,FCR$x1,
     ylim = c(0.0001,max(FCR$total)),
     xlim = c(FCR$datetime[1],FCR$datetime[14]+700000),
     col = 1,
     type = "l",
     lwd = 2,
     ylab = "log (MBFG biovolume um3/mL)",
     xlab = "",
     cex.lab = 1.5,
     cex.axis = 1.3,
     log = "y")
for (i in 3:8){
  points(FCR$datetime,FCR[,i],col = i-1,type = "l",lwd = 2)
}
legend("topright",lty = rep(1,7),lwd = rep(2,7), legend = c(1:7), bty = "n", col = c(1:7),cex = 1.2)
abline(v=as.numeric(FCR$datetime[3]), lwd=2, col='dimgray')
abline(v=as.numeric(FCR$datetime[10]), lwd=2, col='dimgray')
dev.off()

#basic times series plot for FCR - EM2 only
setwd("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/FCR_processed/MBFGs_FCR")
png("FCR_MBFG_biovolume_EM2.png",width = 8,height = 6,units = "in",res = 500)
par(mgp = c(2.3,1,0))
plot(FCR$datetime[9:13],log(FCR$x1[9:13]),
     ylim = c(0,log(max(FCR$total))),
     xlim = c(FCR$datetime[9],FCR$datetime[13]+300000),
     col = 1,
     type = "l",
     lwd = 2,
     ylab = "log (MBFG biovolume um3/mL)",
     xlab = "",
     cex.lab = 1.5,
     cex.axis = 1.3)
for (i in 3:8){
  points(FCR$datetime[9:13],log(FCR[9:13,i]),col = i-1,type = "l",lwd = 2)
}
legend("topright",lty = rep(1,7),lwd = rep(2,7), legend = c(1:7), bty = "n", col = c(1:7),cex = 1.2)
abline(v=as.numeric(BVR$datetime[2]+(0.5*(BVR$datetime[3]-BVR$datetime[2]))), lwd=2, col='dimgray')
dev.off()

#basic times series plot for BVR
setwd("C:/Users/Mary Lofton/Documents/Ch_1/Pelagic_counts/MBFGs_BVR")
png("BVR_MBFG_biovolume.png",width = 8,height = 6,units = "in",res = 500)
par(mgp = c(2.3,1,0))
plot(BVR$datetime,log(BVR$x1),
     ylim = c(0,log(max(BVR$total))),
     xlim = c(BVR$datetime[1],BVR$datetime[5]+300000),
     col = 1,
     type = "l",
     lwd = 2,
     ylab = "log (MBFG biovolume um3/mL)",
     xlab = "",
     cex.lab = 1.5,
     cex.axis = 1.3)
for (i in 3:8){
  points(BVR$datetime,log(BVR[,i]),col = i-1,type = "l",lwd = 2)
}
legend("topright",lty = rep(1,7),lwd = rep(2,7), legend = c(1:7), bty = "n", col = c(1:7),cex = 1.2)
abline(v=as.numeric(BVR$datetime[2]+(0.5*(BVR$datetime[3]-BVR$datetime[2]))), lwd=2, col='dimgray')
dev.off()

#collate data for input into BACI
BVR$Site <- "BVR"
FCR$Site <- "FCR"

FCR_EM2 <- FCR[9:13,]
FCR_EM2$Period <- c(rep("Before",2),rep("After",3))
FCR_EM2$SamplingTime <- c(1:5)
FCR_EM2$SiteClass <- "Impact"
FCR_EM2 <- FCR_EM2[,-9]
FCR_EM2 <- FCR_EM2[,c(11,12,10,1,2:9)]

BVR_EM2 <- BVR
BVR_EM2$Period <- c(rep("Before",2),rep("After",3))
BVR_EM2$SamplingTime <- c(1:5)
BVR_EM2$SiteClass <- "Control"
BVR_EM2 <- BVR_EM2[,-9]
BVR_EM2 <- BVR_EM2[,c(11,12,10,1,2:9)]

final <- rbind(FCR_EM2,BVR_EM2)


for (i in 1:7){
  new_folder <- paste("MBFG_biovolume",i,sep = "_")
  dir.create(paste("./EM2_biovolume",new_folder,sep = "/"))
  write.csv(final,paste("./EM2_biovolume",new_folder,"BACI_MBFGs_biovolume.csv",sep = "/"))
}
