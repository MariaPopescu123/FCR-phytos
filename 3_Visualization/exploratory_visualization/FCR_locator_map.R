#Title: 8C_Fig1_locator_map
#History: created 25APR20 by MEL

#load packages
pacman::p_load(ggmap)

#set local directory
my_directory <- "C:/Users/Mary Lofton/Dropbox/Ch_2/Exploratory_viz/"

#load map
se <- c(left = -83.9, bottom = 24, right = -73.7, top = 40)

#make dataframe with FCR location
points <- data.frame(name = c("FCR"),
                     x = c(-79.837125),
                     y = c(37.308104))

#write map
mymap <- get_stamenmap(se, zoom = 5, maptype = "toner-lite") %>%
  ggmap()+
  geom_point(data = points, aes(x = x, y = y),shape = 23, color = "black",fill = "red", size = 3)+
  geom_text(data = points, aes(x = x, y = y, label = name),hjust = -0.1, vjust = 1.3, fontface = "italic",col = "red")+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
mymap

#save map
ggsave(mymap, filename = file.path(my_directory,paste0("locator_map.tif")),device = "tiff",
       height = 4, width = 2, units = "in", scale = 1, dpi = 300)
