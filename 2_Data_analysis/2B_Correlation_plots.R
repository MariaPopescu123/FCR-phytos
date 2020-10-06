#2B_Correlation_plots
#Author: Mary Lofton
#Date: 06OCT20

mega <- read_csv("./2_Data_analysis/megamatrix.csv")

ggplot(data = mega, aes(x = Peak_depth_m, y = rel_abund_cyano, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_depth_m, y = rel_abund_crypto, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_depth_m, y = rel_abund_green, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_depth_m, y = rel_abund_brown, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_depth_m, y = richness, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_depth_m, y = shannon, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_magnitude_ugL, y = richness, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)+
  theme_classic()
ggplot(data = mega, aes(x = Peak_magnitude_ugL, y = richness))+
  geom_point()+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)+
  theme_classic()
median(mega$Peak_magnitude_ugL, na.rm = TRUE)
mean(mega$Peak_magnitude_ugL, na.rm = TRUE)
sd(mega$Peak_magnitude_ugL, na.rm = TRUE)



ggplot(data = mega, aes(x = Peak_magnitude_ugL, y = shannon, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_magnitude_ugL, y = rel_abund_cyano, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_magnitude_ugL, y = rel_abund_crypto, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_magnitude_ugL, y = rel_abund_green, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_magnitude_ugL, y = rel_abund_brown, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_width_m, y = rel_abund_cyano, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_width_m, y = rel_abund_crypto, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_width_m, y = rel_abund_brown, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_width_m, y = rel_abund_green, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_width_m, y = richness, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()

ggplot(data = mega, aes(x = Peak_width_m, y = shannon, group = as.factor(Year), color = as.factor(Year), fill = as.factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()
