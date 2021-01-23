#3A_ARIMA_results_exploration
#Author: Mary Lofton
#Date: 01DEC20

#List of desired graphics:

#1. Barplots of which covariates are most common in top models (n=20)
#2. Boxplots of coefficient values for covariates in top models (n=20)
#3. Arrange PowerPoint with plots by response var instead of by year (4 slides)
#4. Slide of best model for each response var overall and by year with fit diagnostics (4 slides)

#load packages
pacman::p_load(tidyverse, lubridate)

#read in top models for full timeseries
FP.top.models <- read_csv("./2_Data_analysis/FP_top_models.csv") 
CS.top.models <- read_csv("./2_Data_analysis/CS_top_models.csv")
top.models <- bind_rows(FP.top.models, CS.top.models)

best.models <- top.models %>%
  filter(Rank == 1)

response.vars <- unique(best.models$Response.var)
for (i in 1:length(response.vars)){
  current.var <- response.vars[i]
  current.model <- subset(best.models, best.models$Response.var == current.var)
  print(current.var)
  dat <- current.model[,7:8]
  print(dat)
}

p1 <- ggplot(data = subset(top.models, top.models$Rank == 1), aes(x = Covar, group = Response.var, color = Response.var, fill = Response.var))+
  geom_bar()+
  facet_wrap(vars(Response.var),ncol = 6, nrow = 2)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")
p1
ggsave(plot = p1, filename = "./3_Visualization/top_model_covars_all.png",
       device = "png", height = 6, width = 16, units = "in")

p2 <- ggplot(data = subset(top.models, top.models$Covar != "intercept"), aes(x = Covar.coefs, group = Covar, color = Covar, fill = Covar))+
  geom_boxplot()+
  facet_wrap(vars(Response.var),ncol = 7, nrow = 2, scales = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_vline(xintercept = 0)
p2
ggsave(plot = p2, filename = "./3_Visualization/top_model_covars_coefs_all.png",
       device = "png", height = 6, width = 16, units = "in")

p3 <- ggplot(data = subset(top.models, top.models$Rank == 1), aes(x = Covar))+
  geom_bar()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")
p3
ggsave(plot = p3, filename = "./3_Visualization/top_model_covars_summary.png",
       device = "png", height = 3, width = 4, units = "in")

check <- best.models %>%
  filter(Covar == "Kd")
#read in top models for 2016
top.models.2016 <- read_csv("./2_Data_analysis/top_models_2016.csv")

p1 <- ggplot(data = top.models.2016, aes(x = Covar, group = Response.var, color = Response.var, fill = Response.var))+
  geom_bar()+
  facet_wrap(vars(Response.var),ncol = 2, nrow = 2, scales = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")
p1
ggsave(plot = p1, filename = "./3_Visualization/top_model_covars_2016.png",
       device = "png", height = 4, width = 8, units = "in")

p2 <- ggplot(data = subset(top.models.2016, top.models.2016$Covar != "intercept"), aes(x = Covar.coefs, group = Covar, color = Covar, fill = Covar))+
  geom_boxplot()+
  facet_wrap(vars(Response.var),ncol = 2, nrow = 2, scales = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_vline(xintercept = 0)
p2
ggsave(plot = p2, filename = "./3_Visualization/top_model_covars_coefs_2016.png",
       device = "png", height = 4, width = 8, units = "in")

#read in top models for 2017
top.models.2017 <- read_csv("./2_Data_analysis/top_models_2017.csv")

p1 <- ggplot(data = top.models.2017, aes(x = Covar, group = Response.var, color = Response.var, fill = Response.var))+
  geom_bar()+
  facet_wrap(vars(Response.var),ncol = 2, nrow = 2, scales = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")
p1
ggsave(plot = p1, filename = "./3_Visualization/top_model_covars_2017.png",
       device = "png", height = 4, width = 8, units = "in")

p2 <- ggplot(data = subset(top.models.2017, top.models.2017$Covar != "intercept"), aes(x = Covar.coefs, group = Covar, color = Covar, fill = Covar))+
  geom_boxplot()+
  facet_wrap(vars(Response.var),ncol = 2, nrow = 2, scales = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_vline(xintercept = 0)
p2
ggsave(plot = p2, filename = "./3_Visualization/top_model_covars_coefs_2017.png",
       device = "png", height = 4, width = 8, units = "in")

#read in top models for 2018
top.models.2018 <- read_csv("./2_Data_analysis/top_models_2018.csv")

p1 <- ggplot(data = top.models.2018, aes(x = Covar, group = Response.var, color = Response.var, fill = Response.var))+
  geom_bar()+
  facet_wrap(vars(Response.var),ncol = 2, nrow = 2, scales = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")
p1
ggsave(plot = p1, filename = "./3_Visualization/top_model_covars_2018.png",
       device = "png", height = 4, width = 8, units = "in")

p2 <- ggplot(data = subset(top.models.2018, top.models.2018$Covar != "intercept"), aes(x = Covar.coefs, group = Covar, color = Covar, fill = Covar))+
  geom_boxplot()+
  facet_wrap(vars(Response.var),ncol = 2, nrow = 2, scales = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_vline(xintercept = 0)
p2
ggsave(plot = p2, filename = "./3_Visualization/top_model_covars_coefs_2018.png",
       device = "png", height = 4, width = 8, units = "in")

#read in top models for 2019
top.models.2019 <- read_csv("./2_Data_analysis/top_models_2019.csv")

p1 <- ggplot(data = top.models.2019, aes(x = Covar, group = Response.var, color = Response.var, fill = Response.var))+
  geom_bar()+
  facet_wrap(vars(Response.var),ncol = 2, nrow = 2, scales = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")
p1
ggsave(plot = p1, filename = "./3_Visualization/top_model_covars_2019.png",
       device = "png", height = 4, width = 8, units = "in")

p2 <- ggplot(data = subset(top.models.2019, top.models.2019$Covar != "intercept"), aes(x = Covar.coefs, group = Covar, color = Covar, fill = Covar))+
  geom_boxplot()+
  facet_wrap(vars(Response.var),ncol = 2, nrow = 2, scales = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_vline(xintercept = 0)
p2
ggsave(plot = p2, filename = "./3_Visualization/top_model_covars_coefs_2019.png",
       device = "png", height = 4, width = 8, units = "in")


