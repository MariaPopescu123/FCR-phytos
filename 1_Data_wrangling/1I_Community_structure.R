#1I_Community_structure
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
#10. total abundance of divisions/groups
#11. relative abundance of divisions/groups

#load packages
#install.packages('pacman')
pacman::p_load(tidyverse, lubridate)
rm(list=ls())

dat <- dir(path = "C:/Users/Mary Lofton/Dropbox/Ch_2/BV_concentration/Pelagic/Site_50_QAQC", pattern = paste0(".csv")) %>% 
  map_df(~ read_csv(file.path(path = "C:/Users/Mary Lofton/Dropbox/Ch_2/BV_concentration/Pelagic/Site_50_QAQC", .), col_types = cols(.default = "c")))

dat1 <- dat %>% 
  mutate(BV_um3mL = as.double(BV_um3mL)) %>%
  select(Date, Site, Depth, Genus, BV_um3mL) %>%
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
  mutate(Date = ifelse(Date == "7/10/2017","2017-07-10",Date)) %>%
  mutate(Date = ifelse(Date == "2019-10-18","2019-10-16",Date)) %>%
  filter(!Genus == "unknown")

dat1 <- read_csv("./0_Data_files/phytoplankton.csv")


#calculate total BV for each sample day
total_bv <- dat1 %>% filter(Site == 50) %>%
  mutate(BV_um3mL = as.double(BV_um3mL)) %>%
  group_by(Date) %>%
  summarize(BV_TOTAL = sum(BV_um3mL, na.rm = TRUE)) %>%
  mutate(Date = as.Date(Date)) 

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

#get relative abundances of divisions
dat2 <- dat1 %>%
  filter(Site == 50) %>%
  mutate(Phyto_group = ifelse(Genus %in% cyano,"Cyanobacteria",
                              ifelse(Genus %in% chloro,"Chlorophytes",
                                     ifelse(Genus %in% baci, "Bacillaria",
                                            ifelse(Genus %in% chryso, "Chrysophytes",
                                                   ifelse(Genus %in% dino, "Dinoflagellates",
                                                          ifelse(Genus %in% desmid,"Desmids",
                                                                 ifelse(Genus %in% crypto, "Cryptophytes",
                                                                        ifelse(Genus %in% eugleno, "Euglenoids",
                                                                               ifelse(Genus %in% raphid, "Raphids",""))))))))))

dat3 <- dat2 %>%
  group_by(Date, Phyto_group) %>%
  summarize(BV_group = sum(BV_um3mL, na.rm = TRUE)) %>%
  mutate(Date = as.Date(Date))

dates <- unique(dat3$Date)
groups <- unique(dat3$Phyto_group)
combinations <- expand.grid(Date = dates, Phyto_group = groups)

dat4 <- full_join(dat3, combinations, by = c("Date","Phyto_group")) %>%
  mutate(BV_group = ifelse(is.na(BV_group), 0, BV_group)) 

dat5 <- left_join(dat4, total_bv, by = "Date") %>%
  mutate(rel_abund_group = BV_group/BV_TOTAL,
         Year = year(Date))

dat6 <- dat5 %>%
  ungroup() %>%
  select(Date, BV_group, Phyto_group) %>%
  spread(BV_group,key = Phyto_group)

colnames(dat6)[2:10] <- paste("BV", colnames(dat6)[2:10], sep = "_")

dat7 <- dat5 %>%
  ungroup() %>%
  select(Date, rel_abund_group, Phyto_group) %>%
  spread(rel_abund_group,key = Phyto_group)

colnames(dat7)[2:10] <- paste("rel_abund", colnames(dat7)[2:10], sep = "_")

dat8 <- left_join(dat6, dat7, by = "Date")
dat9 <- left_join(dat8, total_bv, by = "Date")

sample_depths <- dat1 %>%
  rename(Date = Date) %>%
  select(Date, Depth_m) %>%
  distinct() %>%
  rename(Phyto_Depth_m = Depth_m) %>%
  mutate(Date = as.Date(Date))

dat10 <- left_join(dat9, sample_depths, by = "Date")

write.csv(dat10, file = "./0_Data_files/Community_structure.csv",row.names = FALSE)

#plot relative abundance of divisions
p1 <- ggplot(dat5, aes(x = Date, y = rel_abund_group, group = Phyto_group, color = Phyto_group, fill = Phyto_group)) + 
  geom_area(position = "stack") +
  facet_wrap(vars(Year), scales = "free_x")+
  scale_color_brewer(palette = "Set3")+
  scale_fill_brewer(palette = "Set3")+
  ylab("Relative abundance")+
  labs(fill = "Phytoplankton group", color = "Phytoplankton group")+
  theme_classic()
ggsave(plot = p1, filename = "C:/Users/Mary Lofton/Dropbox/Ch_2/Exploratory_viz/relabund.png",
       device = "png",height = 4, width = 7, units = "in")

#plot total biovolume over time
p2 <- ggplot(dat5, aes(x = Date, y = BV_TOTAL)) + 
  geom_line() +
  geom_point(size = 2)+
  facet_wrap(vars(Year), scales = "free_x")+
  ylab(expression(paste("Biovolume ","(",mu,m^3,~mL^-1,")")))+
  theme_classic()
p2
ggsave(plot = p2, filename = "C:/Users/Mary Lofton/Dropbox/Ch_2/Exploratory_viz/BV.png",
       device = "png",height = 2.5, width = 6, units = "in")

#further visualization
cs <- read_csv("./0_Data_files/Community_structure.csv") %>%
  mutate(Year = year(Date))

ggplot(data = cs, aes(x = Date, y = BV_TOTAL))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = BV_Bacillaria))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = BV_Chlorophytes))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = BV_Chrysophytes))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = BV_Cryptophytes))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = BV_Cyanobacteria))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = BV_Desmids))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = BV_Dinoflagellates))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = BV_Euglenoids))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = BV_Raphids))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = rel_abund_Bacillaria))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = rel_abund_Chlorophytes))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = rel_abund_Chrysophytes))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = rel_abund_Cryptophytes))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = rel_abund_Cyanobacteria))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = rel_abund_Desmids))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = rel_abund_Dinoflagellates))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = rel_abund_Euglenoids))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

ggplot(data = cs, aes(x = Date, y = rel_abund_Raphids))+
  facet_wrap(vars(Year), scales = "free_x")+
  geom_line(size = 1)+
  theme_classic()

# autocorrelation
cs1 <- cs %>%
  select(Date:BV_Dinoflagellates, BV_TOTAL,Year) %>%
  gather(BV_Bacillaria:BV_TOTAL, key = "cs_metric", value = "value")

yrs <- unique(cs1$Year)
cs_metrics <- unique(cs1$cs_metric)

png(file = "C:/Users/Mary Lofton/Dropbox/Ch_2/Exploratory_viz/CS_pacf1.png",width = 36, height = 16,
    units = "cm",res = 300)
par(mfrow = c(4,8), mgp = c(2,0.5,0),mar = c(4,3,3,1))

for (j in 1:length(yrs)){
  for (k in 1:length(cs_metrics)){
    
    mydata <- cs1 %>%
      filter(Year == yrs[j],
             cs_metric == cs_metrics[k])
    
    myacf <- acf(mydata$value, 
                 type = "partial",
                 plot = FALSE,
                 na.action = na.pass)
    plot(myacf,main = "")
    title(c(yrs[j],cs_metrics[k]),line = 1)
    
  }
  
}

dev.off()

#autocorrelation for relative abundance
cs2 <- cs %>%
  select(Date,rel_abund_Bacillaria:rel_abund_Dinoflagellates,Year) %>%
  gather(rel_abund_Bacillaria:rel_abund_Dinoflagellates, key = "cs_metric", value = "value")

yrs <- unique(cs2$Year)
cs_metrics <- unique(cs2$cs_metric)

png(file = "C:/Users/Mary Lofton/Dropbox/Ch_2/Exploratory_viz/CS_pacf2.png",width = 36, height = 16,
    units = "cm",res = 300)
par(mfrow = c(4,7), mgp = c(2,0.5,0),mar = c(4,3,3,1))

for (j in 1:length(yrs)){
  for (k in 1:length(cs_metrics)){
    
    mydata <- cs2 %>%
      filter(Year == yrs[j],
             cs_metric == cs_metrics[k])
    
    myacf <- acf(mydata$value, 
                 type = "partial",
                 plot = FALSE,
                 na.action = na.pass)
    plot(myacf,main = "")
    title(c(yrs[j],cs_metrics[k]),line = 1)
    
  }
  
}

dev.off()





