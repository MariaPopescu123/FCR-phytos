#1G_Community_structure
#Author: Mary Lofton
#Date: 17SEP20

#calculate community structure metrics on phytoplankton structure

#1. total BV
#2. total abundance green algae
#3. relative abundance green algae
#4. total abundance brown algae
#5. relative abundance brown algae
#6. total abundance cyanobacteria
#7. relative abuandance cyanobacteria
#8. total abundance cryptophytes
#9. relative abundance cryptophytes

#load packages
#install.packages('pacman')
pacman::p_load(tidyverse, lubridate)
rm(list=ls())

dat <- dir(path = "C:/Users/Mary Lofton/Dropbox/Ch_2/BV_concentration/Pelagic/Site_50_QAQC", pattern = paste0(".csv")) %>% 
  map_df(~ read_csv(file.path(path = "C:/Users/Mary Lofton/Dropbox/Ch_2/BV_concentration/Pelagic/Site_50_QAQC", .), col_types = cols(.default = "c")))

dat1 <- dat %>% 
  mutate(BV_um3mL = as.double(BV_um3mL)) %>%
  select(Sample_date, Site, Depth, Genus, BV_um3mL) %>%
  mutate(Genus = ifelse(Genus == "unicell Chlorophyte spp. 1" | Genus == "unicell Chlorophyte sp. 1" | Genus == "chlorophyte spp. 1" | Genus == "Chlorophyte spp. 1" | Genus == "Chlorella" | Genus == "Chorella","Chlorophyte sp. 1",Genus)) %>%
  mutate(Genus = ifelse(Genus == "unicell Chlorophyte spp. 2" | Genus == "unicell Chlorophyte sp. 2" | Genus == "Chlorophyte spp. 2" | Genus == "Chlorophyte sp. 2" | Genus == "unicell Chlorophyte spp. 6" | Genus == "unicell Chlorophyte sp. 6" | Genus == "Chlorophyte spp. 6" | Genus == "Chlorophyte sp. 6","Rhodomonas",Genus)) %>%
  mutate(Genus = ifelse(Genus == "Chlorophyte spp. 3" | Genus == "unicell Chlorophyte spp. 3","Chlorophyte sp. 3",Genus)) %>%
  mutate(Genus = ifelse(Genus == "unicell Chlorophyte spp. 4","Chlorophyte sp. 4",Genus)) %>%
  mutate(Genus = ifelse(Genus == "chlorophyte spp. 5","Chlorophyte sp. 5",Genus)) %>%
  mutate(Genus = ifelse(Genus == "Euglena?","Euglena",Genus)) %>%
  mutate(Genus = ifelse(Genus == "Desmid spp. 1","Desmid sp. 1",Genus)) %>%
  mutate(Genus = ifelse(Genus == "filament" | Genus == "Psuedanabaena" | Genus == "Pseduanabaena","Pseudanabaena",Genus)) %>%
  mutate(Genus = ifelse(Genus == "dinoflagellate" | Genus == "Gymnodinium?" | Genus == "Gymnodinium??" | Genus == "Gymnodinuim" | Genus == "Gymnodinum" | Genus == "Gymnodinium" | Genus == "naked dino","Parvodinium",Genus)) %>%
  mutate(Genus = ifelse(Genus == "Mononastix","Monomastix",Genus)) %>%
  mutate(Genus = ifelse(Genus == "Woronchinia" | Genus == "Woronochinia","Woronichinia",Genus)) %>%
  mutate(Genus = ifelse(Genus == "Ankistodesmus" ,"Ankistrodesmus",Genus)) %>%
  mutate(Genus = ifelse(Genus == "Aphanothece" ,"Aphanocapsa",Genus)) %>%
  mutate(Genus = ifelse(Genus == "Dolichosphermum","Dolichospermum",Genus)) %>%
  mutate(Genus = ifelse(Genus == "Dysmorphococcus","Gloeodinium",Genus)) %>%
  mutate(Genus = ifelse(Genus == "Selenastrium","Selenastrum",Genus)) %>%
  mutate(Genus = ifelse(Genus == "dinoflagellate cyst","Prorocentrum",Genus)) %>%
  mutate(Genus = ifelse(Genus == "Dicytosphaerium" | Genus == "Dictyospaerium" | Genus == "Dicytospherium","Dictyosphaerium",Genus)) %>%
  mutate(Genus = ifelse(Genus == "Euglenoid/Cryptomonoid" | Genus == "Euglenoid/Cryptomonoid 2" | Genus == "Cryptophyte","Cryptomonas",Genus)) %>%
  mutate(Sample_date = ifelse(Sample_date == "7/10/2017","2017-07-10",Sample_date)) %>%
  mutate(Sample_date = ifelse(Sample_date == "2019-10-18","2019-10-16",Sample_date)) %>%
  filter(!Genus == "unknown")

#calculate total BV for each sample day
total_bv <- dat1 %>% filter(Site == 50) %>%
  mutate(BV_um3mL = as.double(BV_um3mL)) %>%
  group_by(Sample_date) %>%
  summarize(BV_TOTAL = sum(BV_um3mL, na.rm = TRUE)) %>%
  mutate(Date = as.Date(Sample_date)) %>%
  select(-Sample_date)

#spectral groups
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
            "Botryococcus","Kirchneriella","Quadrigula","Actinastrum",
            "Tetraspora","Nephrocytium","Dunaliella","Coccomyxa","Volvulina",
            "Polytomella")
baci <- c("Nitzchia","Synedra","Asterionella","Cyclotella" ,"Fragilaria",
          "Ceratoneis","Tabellaria","Navicula" )
chryso <- c("Synura","Dinobryon","Oochromonas")
dino <- c("Gymnodinium","dinoflagellate cyst","Peridinium","naked dino","Parvodinium","Gloeodinium",
          "Prorocentrum")
desmid <- c("Spondylosium","Staurastrum","Staurodesmus","Closterium",
            "Desmid sp. 1","Bambusina","Actinotaenium","Cosmarium" )
crypto <- c("Cryptomonas","Rhodomonas" )
eugleno <- c("Trachelomonas","Euglena","Euglena?","Lepocinclis" )
raphid <- c("Gonyostomum")

genera <- c(cyano, chloro, baci, chryso, dino, desmid, crypto, eugleno,raphid)
check <- dat1 %>%
  filter(!Genus %in% genera)
bad_genera <- unique(check$Genus)

#separate out and summarize count data by spectral group
green <- dat1 %>% filter(Site == 50 & Genus %in% c(chloro, desmid, eugleno, raphid)) %>%
  mutate(BV_um3mL = as.double(BV_um3mL)) %>%
  arrange(Sample_date) %>%
  group_by(Sample_date) %>%
  summarize(BV_green = sum(BV_um3mL, na.rm = TRUE)) %>%
  mutate(Date = as.Date(Sample_date)) %>%
  select(-Sample_date)

cyano <- dat1 %>% filter(Site == 50 & Genus %in% cyano) %>%
  mutate(BV_um3mL = as.double(BV_um3mL)) %>%
  arrange(Sample_date) %>%
  group_by(Sample_date) %>%
  summarize(BV_cyano = sum(BV_um3mL, na.rm = TRUE)) %>%
  mutate(Date = as.Date(Sample_date))%>%
  select(-Sample_date)

brown <- dat1 %>% filter(Site == 50 & Genus %in% c(baci, dino, chryso)) %>%
  mutate(BV_um3mL = as.double(BV_um3mL)) %>%
  arrange(Sample_date) %>%
  group_by(Sample_date) %>%
  summarize(BV_brown = sum(BV_um3mL, na.rm = TRUE)) %>%
  mutate(Date = as.Date(Sample_date)) %>%
  select(-Sample_date)

crypto <- dat1 %>% filter(Site == 50 & Genus %in% crypto) %>%
  mutate(BV_um3mL = as.double(BV_um3mL)) %>%
  arrange(Sample_date) %>%
  group_by(Sample_date) %>%
  summarize(BV_crypto = sum(BV_um3mL, na.rm = TRUE)) %>%
  mutate(Date = as.Date(Sample_date))%>%
  select(-Sample_date) 

#joining all spectral groups
comm <- left_join(total_bv, green, by = "Date")
comm1 <- left_join(comm, brown, by = "Date")
comm2 <- left_join(comm1, cyano, by = "Date")
comm3 <- left_join(comm2, crypto, by = "Date")
comm3[is.na(comm3)] <- 0

#calculate relative abundance
comm4 <- comm3 %>%
  mutate(rel_abund_green = BV_green/BV_TOTAL,
         rel_abund_brown = BV_brown/BV_TOTAL,
         rel_abund_cyano = BV_cyano/BV_TOTAL,
         rel_abund_crypto = BV_crypto/BV_TOTAL)

comm5 <- comm4[,c(2,1,3:10)]

write.csv(comm5, "./00_Data_files/Community_structure.csv", row.names = FALSE)

#visualization
cs <- read_csv("./00_Data_files/Community_structure.csv") %>%
  mutate(Year = year(Date))

ggplot(data = cs, aes(x = Date, y = BV_TOTAL))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = BV_green))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = BV_brown))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = BV_cyano))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = BV_crypto))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = rel_abund_green))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = rel_abund_brown))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = rel_abund_cyano))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = rel_abund_crypto))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

cs_plot <- cs %>%
  select(Year, Date, rel_abund_green:rel_abund_crypto)%>%
  gather(rel_abund_green:rel_abund_crypto, key = spectral_group, value = rel_abund)

ggplot(cs_plot, aes(x = Date, y = rel_abund, fill = spectral_group)) + 
  geom_area(position = 'stack') +
  facet_wrap(vars(Year), scales = "free_x")+
  theme_classic()
