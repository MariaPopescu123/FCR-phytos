#Fig4
#Author: Mary Lofton
#Date: 25MAR21

pacman::p_load(tidyverse, lubridate,lemon, grid, gridExtra,cowplot,gtable)
rm(list=ls())


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
         Peak_depth_m = ifelse(Peak_depth_m > 9, NA, Peak_depth_m)) %>%
  select(EM2, pz_depth_m, thermo.depth, Peak_depth_m, SRPmax_depth_m, DINmax_depth_m, DOCmax_depth_m) %>%
  gather(pz_depth_m:DOCmax_depth_m, key = "metric",value = "value")

mydatashort <- left_join(my.cs.data,my.fp.data, by = "Date") %>%
  mutate(Year = as.factor(Year), MonthDay = format(Date, format="%m-%d"),
         Peak_depth_m = ifelse(Peak_depth_m > 9, NA, Peak_depth_m)) %>%
  select(EM2, pz_depth_m, thermo.depth, Peak_depth_m, SRPmax_depth_m, DINmax_depth_m, DOCmax_depth_m)

em <- subset(mydata, mydata$EM2 == 1) %>%
  mutate(metric = factor(metric, levels = c("pz_depth_m","thermo.depth","Peak_depth_m","DOCmax_depth_m","DINmax_depth_m","SRPmax_depth_m"))) %>%
  mutate(sig = ifelse(metric == "thermo.depth","*",
                      ifelse(metric == "Peak_depth_m","*",
                             ifelse(metric == "SRPmax_depth_m","*"," "))),
         sig_y = 0)
no_em <- subset(mydata, mydata$EM2 == 0)%>%
  mutate(metric = factor(metric, levels = c("pz_depth_m","thermo.depth","Peak_depth_m","DOCmax_depth_m","DINmax_depth_m","SRPmax_depth_m")))

em_short <- subset(mydatashort, mydatashort$EM2 == 1)
no_em_short <- subset(mydatashort, mydatashort$EM2 == 0)

my.cols <- c("#FFE699","#ACCCEA","#70AD78","#1F4E79","#1F4E79","#1F4E79")

p1 <- ggplot(data = no_em, aes(x = metric, y = value, group = metric, color = metric))+
  geom_boxplot(alpha = 0.5, size = 1)+
  scale_color_manual(values = my.cols, guide = FALSE)+
  scale_x_discrete(labels = c("pz_depth_m" = "Photic zone","thermo.depth" = "Thermocline","Peak_depth_m" = "Biomass peak","DOCmax_depth_m" = "Maximum DOC","DINmax_depth_m" = "Maximum DIN","SRPmax_depth_m" = "Maximum SRP"))+
  # geom_hline(yintercept = median(no_em_short$pz_depth_m, na.rm = TRUE), color = "#FFE699", size = 1)+
  # geom_hline(yintercept = median(no_em_short$thermo.depth, na.rm = TRUE), color = "#ACCCEA", size = 1)+
  # geom_hline(yintercept = median(no_em_short$Peak_depth_m, na.rm = TRUE), color = "#70AD78", size = 1)+
  # geom_hline(yintercept = median(no_em_short$DOCmax_depth_m, na.rm = TRUE), color = "#1F4E79", size = 1)+
  # geom_hline(yintercept = median(no_em_short$DINmax_depth_m, na.rm = TRUE), color = "#1F4E79", size = 1)+
  # geom_hline(yintercept = median(no_em_short$SRPmax_depth_m, na.rm = TRUE), color = "#1F4E79", size = 1)+
  scale_y_reverse()+
  ylim(9.5,0)+
  xlab("")+
  ylab("Depth (m)")+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18), axis.title = element_text(size = 20), axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA))
p1
ggsave(p1, filename = "./3_Visualization/conceptual_fig_reproduction_no_EM.png",
       device = "png",height = 6.5, width = 6, units = "in", dpi = 300, bg = "transparent")


p2 <- ggplot(data = em, aes(x = metric, y = value, group = metric, color = metric))+
  geom_boxplot(alpha = 0.5, size = 1)+
  scale_color_manual(values = my.cols, guide = FALSE)+
  scale_x_discrete(labels = c("pz_depth_m" = "Photic zone","thermo.depth" = "Thermocline","Peak_depth_m" = "Biomass peak","DOCmax_depth_m" = "Maximum DOC","DINmax_depth_m" = "Maximum DIN","SRPmax_depth_m" = "Maximum SRP"))+
  # geom_hline(yintercept = median(em_short$pz_depth_m, na.rm = TRUE), color = "#FFE699", size = 1)+
  # geom_hline(yintercept = median(em_short$thermo.depth, na.rm = TRUE), color = "#ACCCEA", size = 1)+
  # geom_hline(yintercept = median(em_short$Peak_depth_m, na.rm = TRUE), color = "#70AD78", size = 1)+
  # geom_hline(yintercept = median(em_short$DOCmax_depth_m, na.rm = TRUE), color = "#1F4E79", size = 1)+
  # geom_hline(yintercept = median(em_short$DINmax_depth_m, na.rm = TRUE), color = "#1F4E79", size = 1)+
  # geom_hline(yintercept = median(em_short$SRPmax_depth_m, na.rm = TRUE), color = "#1F4E79", size = 1)+
  geom_text(aes(y = sig_y+0.05, x = metric,label = sig), size = 10, color = "black")+
  scale_y_reverse()+
  ylim(9.5,0)+
  xlab("")+
  ylab("Depth (m)")+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18), axis.title = element_text(size = 20), axis.text.x = element_text(size = 14,angle = 45, hjust = 1),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA))
p2
ggsave(p2, filename = "./3_Visualization/conceptual_fig_reproduction_EM.png",
       device = "png",height = 6.5, width = 6, units = "in", dpi = 300, bg = "transparent")


#showing relation btwn depth and magnitude of max
depth <- left_join(my.cs.data,my.fp.data, by = "Date") %>%
  mutate(Year = as.factor(Year), MonthDay = format(Date, format="%m-%d"),
         Peak_depth_m = ifelse(Peak_depth_m > 9, NA, Peak_depth_m)) %>%
  select(EM2, Date, Peak_depth_m, SRPmax_depth_m, DINmax_depth_m, DOCmax_depth_m) %>%
  gather(Peak_depth_m:DOCmax_depth_m, key = "depth_metric",value = "max_depth")%>%
  mutate(var = ifelse(depth_metric == "Peak_depth_m","Biomass",
                      ifelse(depth_metric == "SRPmax_depth_m","SRP",
                             ifelse(depth_metric == "DINmax_depth_m","DIN",
                                    ifelse(depth_metric == "DOCmax_depth_m","DOC",NA))))) 
  
mag <- left_join(my.cs.data,my.fp.data, by = "Date") %>%
  mutate(Year = as.factor(Year), MonthDay = format(Date, format="%m-%d"),
         Peak_depth_m = ifelse(Peak_depth_m > 9, NA, Peak_depth_m)) %>%
  select(EM2, Date, Max_biomass_ugL, SRPmax_ugL, DINmax_ugL, DOCmax_mgL) %>%
  gather(Max_biomass_ugL:DOCmax_mgL, key = "mag_metric", value = "magnitude")%>%
  mutate(var = ifelse(mag_metric == "Max_biomass_ugL","Biomass",
                      ifelse(mag_metric == "SRPmax_ugL","SRP",
                             ifelse(mag_metric == "DINmax_ugL","DIN",
                                    ifelse(mag_metric == "DOCmax_mgL","DOC",NA)))))
  
depthmag <- full_join(depth, mag, by = c("EM2","Date","var")) %>%
  select(EM2, var, max_depth, magnitude)


em <- subset(depthmag, depthmag$EM2 == 1) %>%
  mutate(var = factor(var, levels = c("Biomass","DOC","DIN","SRP"))) 
no_em <- subset(depthmag, depthmag$EM2 == 0) %>%
  mutate(var = factor(var, levels = c("Biomass","DOC","DIN","SRP")))

# New facet label names for var variable
em$var_label <- factor(em$var, labels = c(expression(paste("Biomass (",mu,g~L^-1,")")),
                                          expression(paste("DOC (",mg,~L^-1,")")), 
                                          expression(paste("DIN (",mu,g,~L^-1,")")),
                                          expression(paste("SRP (",mu,g,~L^-1,")"))))
p3 <- ggplot(data = em, aes(x = magnitude, y = max_depth, color = var))+
  geom_point(size = 2)+
  facet_rep_wrap(vars(var_label),nrow = 1, ncol = 4, scales = "free_x",
                 labeller = label_parsed, strip.position = "bottom")+
  scale_y_reverse()+
  ylim(9.5,0)+
  scale_color_manual(values = my.cols[3:6], guide = FALSE)+
  theme_classic()+
  ylab("")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text = element_text(size = 18),
        strip.background = element_blank(),
        strip.text = element_text(size = 13.5),
        strip.placement = "outside",
        axis.title = element_text(size = 20))
p3
ggsave(p3, filename = "./3_Visualization/conceptual_fig_reproduction_EM_mag.png",
       device = "png",height = 6.25, width = 9, units = "in", dpi = 300)


no_em$var_label <- factor(no_em$var, labels = c(expression(paste("Biomass (",mu,g~L^-1,")")),
                                          expression(paste("DOC (",mg,~L^-1,")")), 
                                          expression(paste("DIN (",mu,g,~L^-1,")")),
                                          expression(paste("SRP (",mu,g,~L^-1,")"))))


p4 <- ggplot(data = no_em, aes(x = magnitude, y = max_depth, color = var))+
  geom_point(size = 2)+
  facet_rep_wrap(vars(var_label),nrow = 1, ncol = 4, scales = "free_x",
                 labeller = label_parsed, strip.position = "bottom")+
  scale_y_reverse()+
  ylim(9.5,0)+
  scale_color_manual(values = my.cols[3:6], guide = FALSE)+
  theme_classic()+
  ylab("")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text = element_text(size = 18),
        strip.background = element_blank(),
        strip.text = element_text(size = 13.5),
        strip.placement = "outside",
        axis.title.y = element_text(size = 20))
p4
ggsave(p4, filename = "./3_Visualization/conceptual_fig_reproduction_no_EM_mag.png",
       device = "png",height = 6.25, width = 9, units = "in", dpi = 300)

