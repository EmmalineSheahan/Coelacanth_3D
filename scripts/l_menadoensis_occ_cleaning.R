# Latimeria menadoensis Occurrence cleaning

library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(voluModel)
library(rnaturalearth)
library(rgdal)
library(ggspatial)
library(paletteer)

# creating land for plotting
wanted_crs <- make_EPSG() %>% filter(code == 4326)
target_crs <- wanted_crs$prj4

land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]

# reading in full occurrence
l_menadoensis <- read.csv('./data/L_menadoensis_occurrence.csv')

# 2D points
# removing duplicates
l_menadoensis_2D <- l_menadoensis %>% dplyr::select(Species, Longitude, Latitude)
l_menadoensis_2D <- l_menadoensis_2D[which(!(duplicated(l_menadoensis_2D))),]

# making sure there's only one point per pixel
all_2D_env_indo <- rast('./Envs_Indo/all_2D_env_indo.tif')
l_menadoensis_2D <- downsample(l_menadoensis_2D[,2:3], all_2D_env_indo[[5]])
write.csv(l_menadoensis_2D, file = './data/l_menadoensis_2D_occs.csv')

# 3D points
l_menadoensis_3D <- l_menadoensis %>% dplyr::select(Longitude, Latitude, Depth..m.)

temperature <- rast('./Envs_Indo/temperature_summer_18_Indo.tif')
layerNames <- as.numeric(gsub("X", "", names(temperature)))
l_menadoensis_3D$index <- unlist(lapply(l_menadoensis_3D$Depth..m., 
                                      FUN = function(x) which.min(abs(layerNames - x))))
indices <- unique(l_menadoensis_3D$index)

# Making sure points are unique to pixels by depth
l_menadoensis_down <- data.frame()
for(i in indices){
  tempPoints <- l_menadoensis_3D[l_menadoensis_3D$index==i,]
  tempPoints <- downsample(tempPoints, temperature[[1]])
  tempPoints$depth <- rep(layerNames[[i]], times = nrow(tempPoints))
  l_menadoensis_down <- rbind(l_menadoensis_down, tempPoints)
}
write.csv(l_menadoensis_down, file = './data/l_menadoensis_3D_occs.csv')

# plotting L menadoensis occs
pdf('./Plots/Menadoensis_3D_downsampled_occs_with_depth.pdf')
ggplot(data = land) +
  geom_sf() +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.4, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(110, 150), ylim = c(-15, 10)) +
  geom_point(data = l_menadoensis_down, aes(x = Longitude, y = Latitude, 
                                          color = depth), size = 3.5) +
  scale_colour_paletteer_c("viridis::plasma", trans = "reverse") +
  theme(panel.background = element_rect(fill = "lightblue")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  labs(colour = "Depth (m)")
dev.off()
