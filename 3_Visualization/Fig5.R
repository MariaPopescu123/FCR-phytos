#Fig5
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
pacman::p_load(tidyverse, lubridate, lemon, grid, gridExtra, cowplot)
rm(list=ls())

phytos <- read_csv("./0_Data_files/phytoplankton.csv") 

#calculate total BV for each sample day
total_bv <- phytos %>%
  mutate(BV_um3mL = as.double(BV_um3mL)) %>%
  group_by(Date) %>%
  summarize(BV_TOTAL = sum(BV_um3mL, na.rm = TRUE)) 

my.dates <- read_csv("./2_Data_analysis/CS_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_DO_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") %>%
  select(-MonthDay) %>%
  filter(!is.na(BV_TOTAL)) %>%
  select(Date)

total_bv <- left_join(my.dates, total_bv, by = "Date") 

phytos1 <- left_join(my.dates, phytos, by = "Date")

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
check <- phytos1 %>%
  filter(!Genus %in% genera)

#get relative abundances of divisions
phytos2 <- phytos1 %>%
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

phytos3 <- phytos2 %>%
  group_by(Date, Phyto_group) %>%
  summarize(BV_group = sum(BV_um3mL, na.rm = TRUE))
check <- phytos3 %>%
  filter(is.na(BV_group))

dates <- unique(phytos3$Date)
groups <- unique(phytos3$Phyto_group)
combinations <- expand.grid(Date = dates, Phyto_group = groups)

dat4 <- full_join(phytos3, combinations, by = c("Date","Phyto_group")) %>%
  mutate(BV_group = ifelse(is.na(BV_group), 0, BV_group)) 

dat5 <- left_join(dat4, total_bv, by = "Date") %>%
  mutate(rel_abund_group = BV_group/BV_TOTAL,
         Year = year(Date))

max(log(dat5$BV_TOTAL),na.rm = TRUE)
min(log(dat5$BV_TOTAL),na.rm = TRUE)
min(dat5$BV_TOTAL)

dom <- dat5 %>%
  group_by(Year, Phyto_group) %>%
  summarize(mean_rel_abund = mean(rel_abund_group, na.rm = TRUE))

#plot relative abundance of divisions
p1 <- ggplot(dat5, aes(x = Date, y = rel_abund_group, group = Phyto_group, color = Phyto_group, fill = Phyto_group)) + 
  geom_area(position = "stack") +
  facet_wrap(vars(Year), scales = "free_x")+
  scale_color_manual(values = c("burlywood3","darkgreen","darkgoldenrod2","chocolate1","cadetblue3","chartreuse3","brown1","gray48","plum2"))+
  scale_fill_manual(values = c("burlywood3","darkgreen","darkgoldenrod2","chocolate1","cadetblue3","chartreuse3","brown1","gray48","plum2"))+
  ylab("Relative abundance")+
  labs(fill = "Phytoplankton group", color = "Phytoplankton group")+
  theme_classic()
p1
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

#write plot theme
mytheme1 <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                  legend.key = element_blank(),legend.background = element_blank(),
                  text = element_text(size=16), axis.text.y = element_text(size = 14),
                  panel.border = element_rect(colour = "black", fill = NA),
                  strip.text.x = element_text(face = "bold",hjust = 0),
                  strip.background.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.title.y = element_text(size = 12),
                  plot.margin = unit(c(0, 1, 0, 0), "cm"),
                  panel.spacing = unit(-1, "lines"))

first2 <- dat5 %>%
  filter(year(Date) %in% c(2016:2017)) %>%
  mutate(Date = date(Date))

year_labels <- c(
  "2016" = "A. 2016",
  "2017" = "B. 2017",
  "2018" = "C. 2018",
  "2019" = "D. 2019"
)

p1 <- ggplot(first2, aes(x = Date, y = rel_abund_group, group = Phyto_group, color = Phyto_group, fill = Phyto_group)) + 
  geom_area(position = "stack") +
  facet_rep_wrap(vars(Year), scales = "free_x", labeller = labeller(Year = year_labels))+
  scale_color_manual(values = c("burlywood3","darkgreen","darkgoldenrod2","chocolate1","cadetblue3","chartreuse3","brown1","gray48","plum2"))+
  scale_fill_manual(values = c("burlywood3","darkgreen","darkgoldenrod2","chocolate1","cadetblue3","chartreuse3","brown1","gray48","plum2"))+
  ylab("Relative abundance")+
  labs(fill = "Phytoplankton group", color = "Phytoplankton group")+
  scale_x_date(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  geom_vline(aes(xintercept=date("2016-05-05"),linetype="Storm"), size=1) +
  geom_vline(aes(xintercept=date("2016-05-29"),linetype="EM"), size=1) +
  geom_vline(aes(xintercept=date("2016-06-27"),linetype="EM"), size=1) +
  geom_vline(aes(xintercept=date("2016-07-25"),linetype="EM"), size=1) +
  geom_vline(aes(xintercept=date("2017-05-30"),linetype="EM"), size=1) +
  geom_vline(aes(xintercept=date("2017-07-10"),linetype="EM"), size=1) +
  scale_linetype_manual(name = "Mixing events", values = c(Storm = "solid", EM = "dashed"), guide = FALSE)+
  xlab("")+
  mytheme1
p1

second2 <- dat5 %>%
  filter(year(Date) %in% c(2018:2019)) %>%
  mutate(Date = date(Date))

p2 <- ggplot(second2, aes(x = Date, y = rel_abund_group, fill = Phyto_group)) + 
  geom_area(position = "stack") +
  facet_rep_wrap(vars(Year), scales = "free_x", labeller = labeller(Year = year_labels))+
  # scale_color_manual(values = c("burlywood3","darkgreen","darkgoldenrod2","chocolate1","cadetblue3","chartreuse3","brown1","gray48","plum2"),
  #                    guide = guide_legend(override.aes = list(color = "white")))+
  scale_fill_manual(values = c("burlywood3","darkgreen","darkgoldenrod2","chocolate1","cadetblue3","chartreuse3","brown1","gray48","plum2"), guide = FALSE)+
  ylab("Relative abundance")+
  labs(fill = "Phytoplankton group", color = "Phytoplankton group")+
  scale_x_date(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  geom_vline(aes(xintercept=date("2019-06-08"),linetype="Storm"), size=1) +
  geom_vline(aes(xintercept=date("2016-06-27"),linetype="EM"), size=1) +
    scale_linetype_manual(name = "Mixing events", values = c(Storm = "solid", EM = "dashed"),
                        guide = guide_legend(override.aes = list(color = "white")))+
  xlab("")+
  mytheme1+
  theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")
  )
p2

mytheme2 <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                  legend.key = element_blank(),legend.background = element_blank(),
                  text = element_text(size=16), axis.text.y = element_text(size = 12),
                  axis.text.x = element_text(size = 14),
                  axis.title.y = element_text(size = 12),
                  panel.border = element_rect(colour = "black", fill = NA),
                  strip.text.x = element_blank(),
                  strip.background.x = element_blank(),
                  plot.margin = unit(c(0, 1, 0, 0), "cm"),
                  panel.spacing = unit(-1, "lines"))

year_labels <- c(
  "2016" = "c",
  "2017" = "d",
  "2018" = "g",
  "2019" = "h"
)

#plot total biovolume over time
p3 <- ggplot(data = first2, aes(x = Date, y = log(BV_TOTAL))) + 
  geom_line() +
  geom_point(size = 2)+
  facet_rep_wrap(vars(Year), scales = "free_x", labeller = labeller(Year = year_labels))+
  ylab(expression(paste("log biovolume (",mu,m^3,~mL^-1,")")))+
  xlab("")+
  scale_x_date(expand = c(0,0))+
  coord_cartesian(ylim = c(11,16.5))+
  geom_vline(aes(xintercept=date("2016-05-05"),linetype="Storm"), size=1) +
  geom_vline(aes(xintercept=date("2016-05-29"),linetype="EM"), size=1) +
  geom_vline(aes(xintercept=date("2016-06-27"),linetype="EM"), size=1) +
  geom_vline(aes(xintercept=date("2016-07-25"),linetype="EM"), size=1) +
  geom_vline(aes(xintercept=date("2017-05-30"),linetype="EM"), size=1) +
  geom_vline(aes(xintercept=date("2017-07-10"),linetype="EM"), size=1) +
  scale_linetype_manual(name = "Mixing events", values = c(Storm = "solid", EM = "dashed"))+
  guides(linetype=guide_legend(
    keywidth=0.5,
    keyheight=0.5,
    default.unit="inch")
  )+
  mytheme2
p3

p4 <- ggplot(data = second2, aes(x = Date, y = log(BV_TOTAL))) + 
  geom_line() +
  geom_point(size = 2)+
  facet_rep_wrap(vars(Year), scales = "free_x", labeller = labeller(Year = year_labels))+
  ylab(expression(paste("log biovolume (",mu,m^3,~mL^-1,")")))+
  xlab("")+
  scale_x_date(expand = c(0,0))+
  coord_cartesian(ylim = c(11,16.5))+
  geom_vline(aes(xintercept=date("2019-06-08"),linetype="Storm"), size=1) +
  geom_vline(aes(xintercept=date("2016-06-27"),linetype="EM"), size=1) +
   scale_linetype_manual(name = "Mixing events", values = c(Storm = "solid", EM = "dashed"),
                        guide = guide_legend(override.aes = list(color = "white")))+
  mytheme2+
  theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")
  )
p4

plot<-plot_grid(p1,p3,p2,p4, align='v', vjust=1, scale = 1,
                nrow = 4, ncol = 1,
                rel_heights = c(1.0, 0.5, 1.0, 0.5))
ggsave(plot, filename = "./3_Visualization/Fig_phyto_timeseries.tif",height = 9, width = 10,
       units = "in", dpi = 300, dev = "tiff")

