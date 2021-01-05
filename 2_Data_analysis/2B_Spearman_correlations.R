#2B_Correlation_plots
#Author: Mary Lofton
#Date: 06OCT20

#NEED TO GO BACK THROUGH THIS AND ELIMINATE DATES WHERE THE PEAK_DEPTH_M IS 
#MORE THAN 1.2 M FROM THE PHYTO SAMPLE!!

library(tidyverse)
library(forecast)
library(urca)

mega <- read_csv("./2_Data_analysis/megamatrix.csv")%>%
  select(Max_biomass_ugL:BV_TOTAL, Year,Date) %>%
  filter(!Date %in% as.Date(c("2017-01-19","2017-05-01","2017-05-08","2019-01-21")))

yrs <- c(2016:2019)
final <- list(y2016 = matrix(NA,nrow = 8, ncol = 4),
              y2017 = matrix(NA,nrow = 8, ncol = 4),
              y2018 = matrix(NA,nrow = 8, ncol = 4),
              y2019 = matrix(NA,nrow = 8, ncol = 4))

for (i in 1:4){
  yr.data <- mega %>%
    filter(Year == yrs[i])
  
  for (j in 1:4){
  
    FP.metric <- yr.data[,j]
    
    for (k in 1:8){
      
      other.metric <- yr.data[,k+9]
      my.cor <- cor(FP.metric,other.metric,method = "spearman",use = "complete.obs")
      final[[i]][k,j] <- my.cor
      
    }}}

result <- data.frame(rbind(final[[1]],final[[2]],final[[3]],final[[4]]))
result$Year <- rep(2016:2019, each = 8)
result$Comm_metric <- rep(colnames(mega)[10:17],times = 4)
colnames(result)[1:4] <- c("Max_biomass_ugL","Peak_depth_m","Peak_magnitude_ugL","Peak_width_m")

result <- result %>%
  arrange(Comm_metric) %>%
  gather(Max_biomass_ugL:Peak_width_m, key = "FP_metric",value = "cor")

#peak magnitude and BV_green
ggplot(data = result, aes(x = FP_metric, y = cor, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  facet_wrap(vars(Comm_metric),nrow = 2)+
  geom_hline(yintercept = 0.3)+
  geom_hline(yintercept = -0.3)+
  geom_hline(yintercept = 0, lty = 2)+
  theme(axis.text.x = element_text(angle = 90),
        panel.background = element_blank(),
        panel.grid = element_blank())


for (i in 1:4){
  for (j in 1:8){
    my.plot <- ggplot(data = mega, aes_string(x = mega[,i], y = mega[,j+9], group = as.factor(mega$Year), color = as.factor(mega$Year), fill = as.factor(mega$Year)))+
      geom_point()+
      geom_smooth(method = "lm")+
      xlab(colnames(mega)[i])+
      ylab(colnames(mega)[j+9])+
      theme_classic()
    
    print(my.plot)
  }
}

mega <- mega %>%
  mutate(Year = as.factor(Year))

rho1 <- cor(mega$Peak_depth_m, mega$rel_abund_Cyanobacteria, method = "spearman", use = "complete.obs")
p1 <- ggplot(data = mega, aes(x = Peak_depth_m, y = rel_abund_Cyanobacteria))+
  geom_point()+
  #geom_smooth(method = "lm", se=F)+
  xlab("Peak depth (m)")+
  ylab("Relative abundance of cyanobacteria")+
  theme_classic()
p1
ggsave(plot = p1, filename = "C:/Users/Mary Lofton/Dropbox/Ch_2/Exploratory_viz/peakdepth_v_relabundcyano.png",
       device = "png", height = 3, width = 5, units = "in")

rho2 <- cor(mega$Peak_depth_m, mega$rel_abund_Cryptophytes, method = "spearman", use = "complete.obs")
p2 <- ggplot(data = mega, aes(x = Peak_depth_m, y = rel_abund_Cryptophytes))+
  geom_point()+
  #geom_smooth(method = "lm", se=F)+
  xlab("Peak depth (m)")+
  ylab("Relative abundance of cryptophytes")+
  theme_classic()
p2
ggsave(plot = p2, filename = "C:/Users/Mary Lofton/Dropbox/Ch_2/Exploratory_viz/peakdepth_v_relabundcrypto.png",
       device = "png", height = 3, width = 5, units = "in")

#get right color
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 1
cols = gg_color_hue(n)

rho3 <- cor(mega$Peak_magnitude_ugL, mega$richness, method = "spearman", use = "complete.obs")
p3 <- ggplot(data = mega, aes(x = Peak_magnitude_ugL, y = richness))+
  geom_point(color = cols)+
  #stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se = F)+
  xlab(expression(paste("Peak magnitude ","(",mu,g,~L^-1,")")))+
  ylab("Species richness")+
  theme_classic()
p3
ggsave(plot = p3, filename = "C:/Users/Mary Lofton/Dropbox/Ch_2/Exploratory_viz/peakmagnitude_v_richness.png",
       device = "png", height = 2, width = 3, units = "in")

rho4 <- cor(mega$Max_biomass_ugL, mega$richness, method = "spearman", use = "complete.obs")
p4 <- ggplot(data = mega, aes(x = Max_biomass_ugL, y = richness))+
  geom_point(color = cols)+
  #stat_smooth(method = "lm", size = 1, se = F)+
  xlab(expression(paste("Maximum biomass ","(",mu,g,~L^-1,")")))+
  ylab("Species richness")+
  theme_classic()
ggsave(plot = p4, filename = "C:/Users/Mary Lofton/Dropbox/Ch_2/Exploratory_viz/maxbiomass_v_richness.png",
       device = "png", height = 2, width = 3, units = "in")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 4
cols = gg_color_hue(n)

dev.new(width = 4, height = 4)
plot(1:n, pch = 16, cex = 2, col = cols)
