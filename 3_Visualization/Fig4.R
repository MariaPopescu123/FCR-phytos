#Fig3
#Author: Mary Lofton
#Date: 23FEB21 ####
pacman::p_load(cowplot, generics, ggridges, kSamples, tidyverse, lubridate,viridis,grid)
rm(list=ls())

####SET-UP####
#get data
my.fp.data <- read_csv("./2_Data_analysis/FP_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_DO_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") 
colnames(my.fp.data)
my.fp.data <- my.fp.data[,c(2,1,3:48)]

my.cs.data <- read_csv("./2_Data_analysis/CS_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_DO_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") %>%
  select(Date,shannon:BV_TOTAL) %>%
  filter(!is.na(BV_TOTAL)) 

mydata <- left_join(my.cs.data,my.fp.data, by = "Date") %>%
  mutate(Year = as.factor(Year), MonthDay = format(Date, format="%m-%d"),
         Peak_depth_m = ifelse(Peak_depth_m > 9, NA, Peak_depth_m)) 

colnames(mydata)

#limit to numeric physicochemical and phyto variables
vars <- mydata[,c(28,56,52,53,37,45,42,32,36,64,11,20)] 
colnames(vars)

vars_long <- vars %>%
  gather(thermo.depth:rel_abund_Desmids, key = "var", value = "value") %>%
  mutate(EM2 = ifelse(EM2 == 0, "reference summers","manipulation summers"))

my.cols <- viridis(12) 

mytheme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.background = element_blank(), axis.line.x = element_line(colour = 'black'), 
                 axis.line.y = element_line(colour = 'black'), 
                 axis.text.x=element_text(size=16, colour='black'), 
                 axis.text.y=element_text(size=16, colour='black'),
                 plot.title=element_text(size=18, colour = 'black',face = 'bold'),
                 axis.title.x=element_text(size=18), axis.title.y=element_text(size=18),
                 legend.text = element_text(size=18), legend.title=element_text(size=18),
                 strip.text.x = element_text(margin = margin(.25,0,.25,0, 'cm'), size=18),
                 strip.text.y = element_text(margin = margin(0,.25,0,.25, 'cm'), size=18))

# PLOTTING ####

# set annotations
annotation1 <- annotation_custom(grid::textGrob(label = "***",
                                               x = unit(0.9, "npc"), y = unit(0.9, "npc"),
                                               gp = grid::gpar(cex = 3)))
annotation2 <- annotation_custom(grid::textGrob(label = "**",
                                                x = unit(0.9, "npc"), y = unit(0.9, "npc"),
                                                gp = grid::gpar(cex = 3)))
annotation3 <- annotation_custom(grid::textGrob(label = "*",
                                                x = unit(0.9, "npc"), y = unit(0.9, "npc"),
                                                gp = grid::gpar(cex = 3)))
####NEW FIGURE 3####
# thermocline depth ####
thermo.depth <- ggplot(data=subset(vars_long, var == "thermo.depth"), 
             aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("A.")+
  annotation1 +
  xlab("Thermocline depth (m)")

thermo.depth

# schmidt stability ####
schmidt.stability <- ggplot(data=subset(vars_long, var == "schmidt.stability"), 
                       aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("B.")+
  annotation1 +
  xlab(expression(paste("Schmidt stability ","(",J,~m^-2,")")))

schmidt.stability

# n2 ####
n2 <- ggplot(data=subset(vars_long, var == "n2"), 
                            aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("C.")+
  annotation1 +
  xlab(expression(paste("Buoyancy frequency ","(",~s^-1,")")))
  
n2

# SRPmax_depth_m ####
SRPmax_depth_m <- ggplot(data=subset(vars_long, var == "SRPmax_depth_m"), 
             aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("D.")+
  annotation1 +
  xlab(expression(paste("Depth of maximum SRP (",mu,g,~L^-1,")")))
  
SRPmax_depth_m

# pz_DOC_mean ####
pz_DOC_mean <- ggplot(data=subset(vars_long, var == "pz_DOC_mean"), 
                         aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("E.")+
  annotation1 +
  xlab(expression(paste("Mean photic zone DOC (",mg,~L^-1,")")))

pz_DOC_mean

# DOCmax_mgL ####
DOCmax_mgL <- ggplot(data=subset(vars_long, var == "DOCmax_mgL"), 
                      aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("F.")+
  annotation1 +
  xlab(expression(paste("Maximum photic zone DOC (",mg,~L^-1,")")))

DOCmax_mgL

# Cmax_DOC_mgL ####
Cmax_DOC_mgL <- ggplot(data=subset(vars_long, var == "Cmax_DOC_mgL"), 
                     aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("G.")+
  annotation2 +
  xlab(expression(paste("DOC at peak depth (",mg,~L^-1,")")))

Cmax_DOC_mgL

# Grab_DOC_mgL ####
Grab_DOC_mgL <- ggplot(data=subset(vars_long, var == "Grab_DOC_mgL"), 
                       aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("H.")+
  annotation1 +
  xlab(expression(paste("DOC at grab sample depth (",mg,~L^-1,")")))

Grab_DOC_mgL

# Peak_depth_m ####
Peak_depth_m <- ggplot(data=subset(vars_long, var == "Peak_depth_m"), 
                       aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("I.")+
  annotation2 +
  xlab("Depth of peak biomass (m)")
  
Peak_depth_m

# BV_Desmids ####
BV_Desmids <- ggplot(data=subset(vars_long, var == "BV_Desmids"), 
                       aes(x= log(value), fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("J.")+
  annotation3 +
  xlab(expression(paste("log(Desmid biovolume ","(",mu,m^3,~mL^-1,"))")))
  
BV_Desmids

# rel_abund_Desmids ####
rel_abund_Desmids <- ggplot(data=subset(vars_long, var == "rel_abund_Desmids"), 
                     aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("K.")+
  annotation3 +
  xlab("relative abundance of Desmids")

rel_abund_Desmids

# legend.plot ####
legend.plot <- ggplot(data=subset(vars_long, var == "thermo.depth"), 
                            aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme+
  ylab("")+
  ggtitle("L.")+
  annotation2 +
  xlab("legend.plot")+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0.5, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))

legend.plot

boxplot_legend <- cowplot::get_legend(legend.plot)

grid.newpage()
grid.draw(boxplot_legend)

new_fig <- plot_grid(thermo.depth, schmidt.stability, n2, SRPmax_depth_m, pz_DOC_mean, DOCmax_mgL, 
          Cmax_DOC_mgL, Grab_DOC_mgL, Peak_depth_m, BV_Desmids, rel_abund_Desmids, boxplot_legend,
          nrow = 4, ncol = 3, align = "hv", scale = 1)

ggsave(new_fig, filename = "./3_Visualization/FigR1.tif",height = 13.5, width = 15,
       units = "in", dpi = 300, dev = "tiff" )

####Fig S3####
#limit to numeric physicochemical and phyto variables
vars <- mydata[,c(71,25,28,56,52,53,37,45,42,32,36,64,11,20)] 
colnames(vars)

vars_long <- vars %>%
  gather(thermo.depth:rel_abund_Desmids, key = "var", value = "value") %>%
  mutate(EM2 = ifelse(EM2 == 0, "reference summers","manipulation summers"))


# thermocline depth ####
thermo.depth <- ggplot(data=subset(vars_long, var == "thermo.depth"), 
                       aes(x= MonthDay, y = value, fill=as.factor(Year), group = as.factor(Year))) +
  geom_line(size = 1.5)+
  geom_point(size = 3)+
  scale_fill_manual(values = c(my.cols[1],my.cols[4],my.cols[7],my.cols[10]))+
  mytheme +
  theme(legend.position = 'none') +
  ylab("")+
  ggtitle("A.")+
  xlab("Thermocline depth (m)")

thermo.depth

# schmidt stability ####
schmidt.stability <- ggplot(data=subset(vars_long, var == "schmidt.stability"), 
                            aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_line(size = 1.5)+
  geom_point(size = 3)+
  scale_fill_manual(values = c(my.cols[1],my.cols[4],my.cols[7],my.cols[10]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("B.")+
  annotation1 +
  xlab(expression(paste("Schmidt stability ","(",J,~m^-2,")")))

schmidt.stability

# n2 ####
n2 <- ggplot(data=subset(vars_long, var == "n2"), 
             aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_line(size = 1.5)+
  geom_point(size = 3)+
  scale_fill_manual(values = c(my.cols[1],my.cols[4],my.cols[7],my.cols[10]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("C.")+
  annotation1 +
  xlab(expression(paste("Buoyancy frequency ","(",~s^-1,")")))

n2

# SRPmax_depth_m ####
SRPmax_depth_m <- ggplot(data=subset(vars_long, var == "SRPmax_depth_m"), 
                         aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("D.")+
  annotation1 +
  xlab(expression(paste("Depth of maximum SRP (",mu,g,~L^-1,")")))

SRPmax_depth_m

# pz_DOC_mean ####
pz_DOC_mean <- ggplot(data=subset(vars_long, var == "pz_DOC_mean"), 
                      aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("E.")+
  annotation1 +
  xlab(expression(paste("Mean photic zone DOC (",mg,~L^-1,")")))

pz_DOC_mean

# DOCmax_mgL ####
DOCmax_mgL <- ggplot(data=subset(vars_long, var == "DOCmax_mgL"), 
                     aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("F.")+
  annotation1 +
  xlab(expression(paste("Maximum photic zone DOC (",mg,~L^-1,")")))

DOCmax_mgL

# Cmax_DOC_mgL ####
Cmax_DOC_mgL <- ggplot(data=subset(vars_long, var == "Cmax_DOC_mgL"), 
                       aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("G.")+
  annotation2 +
  xlab(expression(paste("DOC at peak depth (",mg,~L^-1,")")))

Cmax_DOC_mgL

# Grab_DOC_mgL ####
Grab_DOC_mgL <- ggplot(data=subset(vars_long, var == "Grab_DOC_mgL"), 
                       aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("H.")+
  annotation1 +
  xlab(expression(paste("DOC at depth sample (",mg,~L^-1,")")))

Grab_DOC_mgL

# Peak_depth_m ####
Peak_depth_m <- ggplot(data=subset(vars_long, var == "Peak_depth_m"), 
                       aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("I.")+
  annotation2 +
  xlab("Depth of peak biomass (m)")

Peak_depth_m

# BV_Desmids ####
BV_Desmids <- ggplot(data=subset(vars_long, var == "BV_Desmids"), 
                     aes(x= log(value), fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("J.")+
  annotation3 +
  xlab(expression(paste("log(Desmid biovolume ","(",mu,m^3,~mL^-1,"))")))

BV_Desmids

# rel_abund_Desmids ####
rel_abund_Desmids <- ggplot(data=subset(vars_long, var == "rel_abund_Desmids"), 
                            aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme +
  theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab("")+
  ggtitle("K.")+
  annotation3 +
  xlab("relative abundance of Desmids")

rel_abund_Desmids

# legend.plot ####
legend.plot <- ggplot(data=subset(vars_long, var == "thermo.depth"), 
                      aes(x= value, fill=as.factor(EM2), group = as.factor(EM2))) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c(my.cols[5],my.cols[11]))+
  mytheme+
  ylab("")+
  ggtitle("L.")+
  annotation2 +
  xlab("legend.plot")+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0.5, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))

legend.plot

boxplot_legend <- cowplot::get_legend(legend.plot)

grid.newpage()
grid.draw(boxplot_legend)

new_fig <- plot_grid(thermo.depth, schmidt.stability, n2, SRPmax_depth_m, pz_DOC_mean, DOCmax_mgL, 
                     Cmax_DOC_mgL, Grab_DOC_mgL, Peak_depth_m, BV_Desmids, rel_abund_Desmids, boxplot_legend,
                     nrow = 4, ncol = 3, align = "hv", scale = 1)

ggsave(new_fig, filename = "./3_Visualization/FigR1.tif",height = 13.5, width = 15,
       units = "in", dpi = 300, dev = "tiff" )
