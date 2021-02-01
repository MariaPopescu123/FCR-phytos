#2B_Correlation_plots
#Author: Mary Lofton
#Date: 06OCT20

library(tidyverse)
library(forecast)
library(urca)

#get data
my.fp.data <- read_csv("./2_Data_analysis/FP_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") %>%
  select(-MonthDay)
my.fp.data <- my.fp.data[,c(2,1,3:22)]

my.cs.data <- read_csv("./2_Data_analysis/CS_megamatrix.csv") %>%
  select(-Chem_Depth_m,-Temp_Depth_m) %>%
  mutate(MonthDay = format(Date, format="%m-%d")) %>%
  filter(MonthDay >= "05-01" & MonthDay <= "09-20") %>%
  select(-MonthDay) %>%
  select(Date,shannon:BV_TOTAL)

mydata <- left_join(my.fp.data,my.cs.data, by = "Date")

mega <- mydata %>%
  select(Max_biomass_ugL:BV_TOTAL, Year,Date) 

check <- cor(mega$Peak_magnitude_ugL,mega$Max_biomass_ugL,method = "spearman",use = "complete.obs")

# yrs <- c(2016:2019)
# final <- list(y2016 = matrix(NA,nrow = 9, ncol = 4),
#               y2017 = matrix(NA,nrow = 9, ncol = 4),
#               y2018 = matrix(NA,nrow = 9, ncol = 4),
#               y2019 = matrix(NA,nrow = 9, ncol = 4))
final <- matrix(NA,nrow = 9, ncol = 4)

  
  for (j in 1:4){
    
    FP.metric <- mega[,j]
    
    for (k in 1:9){
      
      other.metric <- mega[,k+4]
      my.cor <- cor(FP.metric,other.metric,method = "spearman",use = "complete.obs")
      final[k,j] <- my.cor
      
    }}

# for (i in 1:4){
#   yr.data <- mega %>%
#     filter(Year == yrs[i])
#   
#   for (j in 1:4){
#   
#     FP.metric <- yr.data[,j]
#     
#     for (k in 1:9){
#       
#       other.metric <- yr.data[,k+4]
#       my.cor <- cor(FP.metric,other.metric,method = "spearman",use = "complete.obs")
#       final[[i]][k,j] <- my.cor
#       
#     }}}

result <- data.frame(final)
colnames(result)[1:4] <- c("Max_biomass_ugL","Peak_depth_m","Peak_magnitude_ugL","Peak_width_m")
#result[abs(result) < 0.4] <- NA
result <- round(result, digits = 2)
result$Comm_metric <- colnames(mega)[5:13]
write.csv(result, "./2_Data_analysis/spearman_results.csv",row.names = FALSE)
# result <- data.frame(rbind(final[[1]],final[[2]],final[[3]],final[[4]]))
# result$Year <- rep(2016:2019, each = 9)
# result$Comm_metric <- rep(colnames(mega)[5:13],times = 4)
# colnames(result)[1:4] <- c("Max_biomass_ugL","Peak_depth_m","Peak_magnitude_ugL","Peak_width_m")

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
  for (j in 1:9){
    my.plot <- ggplot(data = mega, aes_string(x = mega[,i], y = mega[,j+4], group = as.factor(mega$Year), color = as.factor(mega$Year), fill = as.factor(mega$Year)))+
      geom_point()+
      geom_smooth(method = "lm")+
      xlab(colnames(mega)[i])+
      ylab(colnames(mega)[j+4])+
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
