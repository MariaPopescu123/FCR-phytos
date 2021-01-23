#Network Correlation Plots


#transform correlation matrix into adjacency matrix for network plot
driver_cor <- data.frame(driver_cor)

dc <- driver_cor %>%
  mutate(source = row.names(driver_cor)) %>%
  gather(Year:Kd, key = "target", value = "importance") %>%
  filter(!importance %in% c(1,0,"")) %>%
  mutate(importance = round(abs(as.numeric(importance)*10),0)) %>%
  filter(!importance == 0) %>%
  filter(importance >=5) %>%
  mutate(importance = (importance-5)*3)

names <- colnames(mydata)[2:17]

#set up nodes
nodes <- data.frame(
  name=names,
  carac=c("year",rep("mgmt",2),rep("nutrients",9),rep("physics",4)))

# build the network graph object
network <- graph_from_data_frame(d=dc, vertices=nodes, directed=F)
# Make a palette of 3 colors
coul  <- brewer.pal(4, "Set1") 
# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(network)$carac))]
#plot
resfactor = 3
png(filename = "./3_Visualization/driver_cor.png", 
    height = 600*resfactor, width = 600*resfactor, res = 72*resfactor)
plot(network, vertex.color=my_color,edge.width=E(network)$importance,
     layout = layout.circle,
     vertex.size = 35,
     vertex.label.color="black",
     vertex.label.cex = 0.9)
dev.off()

data.2016 <- mydata[c(1:23),] %>%
  select(-HOx, -Year) 

#look at correlations in 2016 - cutoff of 0.5
driver_cor.2016 <- cor(data.2016[,2:15],method = "spearman",use = "complete.obs")
driver_cor.2016[lower.tri(driver_cor.2016)] = ""

#DOCmax_depth_m and EM
#schmidt.stability and EM
#DINmax_depth_m and Cmax_SRP_ugL
#DINmax_ugL and Cmax_DOC_mgL
#DOCmax_mgL and Cmax_DOC_mgL
#Temp_C and Cmax_DOC_mgL
#schmidt.stability and SRPmax_ugL
#DINmax_ugL and Kd
#DINmax_ugL and DOCmax_depth_m
#DINmax_ugL and DOCmax_mgL
#DINmax_ugL and Temp_C
#thermo.depth and DOCmax_depth_m
#DOCmax_depth_m and Temp_C
#DOCmax_ugL and Temp_C
#Temp_C and thermo.depth
#DOCmax_ugL and Kd

#transform correlation matrix into adjacency matrix for network plot
driver_cor.2016 <- data.frame(driver_cor.2016)

dc_2016 <- driver_cor.2016 %>%
  mutate(source = row.names(driver_cor.2016)) %>%
  gather(EM:Kd, key = "target", value = "importance") %>%
  filter(!importance %in% c(1,0,"")) %>%
  mutate(importance = round(abs(as.numeric(importance)*10),0)) %>%
  filter(!importance == 0) %>%
  filter(importance >=5) %>%
  mutate(importance = (importance-5)*3)

names <- colnames(data.2016)[2:15]

#set up nodes
nodes <- data.frame(
  name=names,
  carac=c("mgmt",rep("nutrients",9),rep("physics",4)))

# build the network graph object
network.2016 <- graph_from_data_frame(d=dc_2016, vertices=nodes, directed=F)
# Make a palette of 3 colors
coul  <- brewer.pal(3, "Set1") 
# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(network.2017)$carac))]
#plot
resfactor = 3
png(filename = "./3_Visualization/driver_cor_2016.png", 
    height = 600*resfactor, width = 600*resfactor, res = 72*resfactor)
plot(network.2016, vertex.color=my_color,edge.width=E(network.2016)$importance,
     layout = layout.circle,
     vertex.size = 35,
     vertex.label.color="black",
     vertex.label.cex = 0.9)
dev.off()

data.2017 <- mydata[c(25:46),] %>%
  select(-HOx, -Year)

#look at correlations in 2017 - cutoff of 0.5
driver_cor.2017 <- cor(data.2017[,c(2:15)],method = "spearman",use = "complete.obs")
driver_cor.2017[lower.tri(driver_cor.2017)] = ""
#EM and Cmax_SRP_ugL
#Cmax_SRP_ugL and Temp_C
#Cmax_DOC_mgL and Temp_C
#Cmax_DOC_mgL and schmidt.stability
#SRPmax_depth_m and DOCmax_depth_m
#SRPmax_depth_m and Temp_C
#SRPmax_ugL and DOC_max_mgL
#DINmax_depth_m and schmidt.stability
#DINmax_depth_m and thermo.depth
#DINmax_depth_m and Kd
#DINmax_ugL and DOCmax_depth_m
#Temp_C and schmidt.stability
#thermo.depth and Kd

#transform correlation matrix into adjacency matrix for network plot
driver_cor.2017 <- data.frame(driver_cor.2017)

dc_2017 <- driver_cor.2017 %>%
  mutate(source = row.names(driver_cor.2017)) %>%
  gather(EM:Kd, key = "target", value = "importance") %>%
  filter(!importance %in% c(1,0,"")) %>%
  mutate(importance = round(abs(as.numeric(importance)*10),0)) %>%
  filter(!importance == 0) %>%
  filter(importance >=5) %>%
  mutate(importance = (importance-5)*3)

names <- colnames(data.2017)[2:15]

#set up nodes
nodes <- data.frame(
  name=names,
  carac=c("mgmt",rep("nutrients",9),rep("physics",4)))

# build the network graph object
network.2017 <- graph_from_data_frame(d=dc_2017, vertices=nodes, directed=F)
# Make a palette of 3 colors
coul  <- brewer.pal(3, "Set1") 
# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(network.2017)$carac))]
#plot
resfactor = 3
png(filename = "./3_Visualization/driver_cor_2017.png", 
    height = 600*resfactor, width = 600*resfactor, res = 72*resfactor)
plot(network.2017, vertex.color=my_color,edge.width=E(network.2017)$importance,
     layout = layout.circle,
     vertex.size = 35,
     vertex.label.color="black",
     vertex.label.cex = 0.9)
dev.off()

data.2018 <- mydata[c(48:67),] %>%
  select(-EM,-Year)

#look at correlations in 2018 - cutoff of 0.5
driver_cor.2018 <- cor(data.2018[,c(2:15)],method = "spearman",use = "complete.obs")
driver_cor.2018[lower.tri(driver_cor.2018)] = ""
#HOx and Cmax_DIN_ugL
#HOx and Cmax_DOC_mgL
#HOx and SRPmax_depth_m
#HOx and DINmax_depth_m
#HOx and DINmax_ugL
#HOx and DOCmax_depth_m
#HOx and DOCmax_mgL
#HOx and thermo.depth
#Cmax_SRP_ugL and SRPmax_ugL
#Cmax_SRP_ugL and Kd
#Cmax_DIN_ugL and DINmax_ugL
#Cmax_DIN_ugL and Temp_C
#Cmax_DIN_ugL and schmidt.stability
#Cmax_DOC_mgL and DINmax_depth_m
#Cmax_DOC_mgL and DINmax_ugL
#Cmax_DOC_mgL and DOCmax_depth_m
#Cmax_DOC_mgL and DOCmax_mgL
#Cmax_DOC_mgL and Temp_C
#Cmax_DOC_mgL and thermo.depth
#SRP_max_depth_m and DINmax_ugL
#SRP_max_depth_m and schmidt.stability
#SRPmax_ugL and Kd
#DINmax_depth_m and DINmax_ugL
#DINmax_depth_m and DOCmax_depth_m
#DINmax_depth_m and DOCmax_mgL
#DINmax_depth_m and thermo.depth
#DINmax_ugL and DOCmax_mgL
#DINmax_ugL and schmidt.stability
#DINmax_ugL and thermo.depth
#DOCmax_depth_m and DOCmax_mgL
#DOCmax_depth_m and thermo.depth
#DOCmax_mgL and thermo.depth
#schmidt.stability and thermo.depth

#transform correlation matrix into adjacency matrix for network plot
driver_cor.2018 <- data.frame(driver_cor.2018)

dc_2018 <- driver_cor.2018 %>%
  mutate(source = row.names(driver_cor.2018)) %>%
  gather(HOx:Kd, key = "target", value = "importance") %>%
  filter(!importance %in% c(1,0,"")) %>%
  mutate(importance = round(abs(as.numeric(importance)*10),0)) %>%
  filter(!importance == 0) %>%
  filter(importance >=5) %>%
  mutate(importance = (importance-5)*3)

names <- colnames(data.2018)[2:15]

#set up nodes
nodes <- data.frame(
  name=names,
  carac=c("mgmt",rep("nutrients",9),rep("physics",4)))

# build the network graph object
network.2018 <- graph_from_data_frame(d=dc_2018, vertices=nodes, directed=F)
# Make a palette of 3 colors
coul  <- brewer.pal(3, "Set1") 
# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(network.2018)$carac))]
#plot
resfactor = 3
png(filename = "./3_Visualization/driver_cor_2018.png", 
    height = 600*resfactor, width = 600*resfactor, res = 72*resfactor)
plot(network.2018, vertex.color=my_color,edge.width=E(network.2018)$importance,
     layout = layout.circle,
     vertex.size = 35,
     vertex.label.color="black",
     vertex.label.cex = 0.9)
dev.off()

data.2019 <- mydata[c(75:104),] %>%
  select(-EM,-Year)

#look at correlations in 2019 - cutoff of 0.5
driver_cor.2019 <- cor(data.2019[,c(2:15)],method = "spearman",use = "complete.obs")
driver_cor.2019[lower.tri(driver_cor.2019)] = ""
#HOx and SRPmax_ugL
#Cmax_SRP_ugL and Cmax_DOC_mgL
#Cmax_SRP_ugL and SRPmax_depth_m
#Cmax_SRP_ugL and SRPmax_ugL
#Cmax_SRP_ugL and DOCmax_mgL
#Cmax_SRP_ugL and Temp_C
#Cmax_DIN_ugL and SRPmax_depth_m
#Cmax_DIN_ugL and thermo.depth
#Cmax_DOC_mgL and SRPmax_depth_m
#Cmax_DOC_mgL and DOCmax_mgL
#Cmax_DOC_mgL and thermo.depth
#SRPmax_depth_m and DINmax_ugL
#SRPmax_depth_m and thermo.depth
#SRPmax_ugL and Temp_C
#DINmax_ugL and DOCmax_mgL
#DINmax_ugL and thermo.depth
#Temp_C and schmidt.stability
#Temp_C and Kd
#schmidt.stability and Kd

#transform correlation matrix into adjacency matrix for network plot
driver_cor.2019 <- data.frame(driver_cor.2019)

dc_2019 <- driver_cor.2019 %>%
  mutate(source = row.names(driver_cor.2019)) %>%
  gather(HOx:Kd, key = "target", value = "importance") %>%
  filter(!importance %in% c(1,0,"")) %>%
  mutate(importance = round(abs(as.numeric(importance)*10),0)) %>%
  filter(!importance == 0) %>%
  filter(importance >=5) %>%
  mutate(importance = (importance-5)*3)

names <- colnames(data.2019)[2:15]

#set up nodes
nodes <- data.frame(
  name=names,
  carac=c("mgmt",rep("nutrients",9),rep("physics",4)))

# build the network graph object
network.2019 <- graph_from_data_frame(d=dc_2019, vertices=nodes, directed=F)
# Make a palette of 3 colors
coul  <- brewer.pal(3, "Set1") 
# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(network.2019)$carac))]
#plot
resfactor = 3
png(filename = "./3_Visualization/driver_cor_2019.png", 
    height = 600*resfactor, width = 600*resfactor, res = 72*resfactor)
plot(network.2019, vertex.color=my_color,edge.width=E(network.2019)$importance,
     layout = layout.circle,
     vertex.size = 35,
     vertex.label.color="black",
     vertex.label.cex = 0.9)
dev.off()


