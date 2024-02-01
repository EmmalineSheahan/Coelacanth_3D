# Formatting and thinning the full Ceolacanth occurrence data set

library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(voluModel)
library(rnaturalearth)
library(rgdal)
library(ggspatial)
library(paletteer)
library(spThin)

# creating wanted crs and land polygon
wanted_crs <- make_EPSG() %>% filter(code == 4326)
target_crs <- wanted_crs$prj4

land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]

# reading in the dataset
all_occs <- 
  read.csv('./data/coelacanth_complete_occurrence.csv')

# filtering for only Latimeria chalumnae
l_chalumnae <- all_occs %>% filter(Species == "Latimeria chalumnae")

# 3D occurrences
# filtering for only records which have depth information
l_chalumnae <- l_chalumnae %>% filter(!(is.na(Depth..m.)))

# creating spatial points
l_chalumnae_sp <- l_chalumnae
coordinates(l_chalumnae_sp) <- ~Longitude + Latitude
l_chalumnae_sp <- st_as_sf(l_chalumnae_sp)
st_crs(l_chalumnae_sp) <- target_crs

# occurrence plot
plot(l_chalumnae_sp$geometry, col = "red")
plot(land, col = NA, add = T)

# lon lat and depth only
l_chalumnae_occs <- l_chalumnae %>% dplyr::select(c(Longitude, Latitude, 
                                                    Depth..m.))
colnames(l_chalumnae_occs) <- c("longitude", "latitude", "depth")

# reading in a WOA 1/4 degree raster layer for spatial thinning
temperature <- rast('./WOA_18/temperature_summer_18.tif')

# thinning occurrences
# Get the layer index for each occurrence by matching to depth
layerNames <- as.numeric(gsub("X", "", names(temperature)))
l_chalumnae_occs$index <- unlist(lapply(l_chalumnae_occs$depth, 
                          FUN = function(x) which.min(abs(layerNames - x))))
indices <- unique(l_chalumnae_occs$index)

# downsampling occurrences
l_chalumnae_down <- data.frame()
for(i in indices){
  tempPoints <- l_chalumnae_occs[l_chalumnae_occs$index==i,]
  tempPoints <- downsample(tempPoints, temperature[[1]])
  tempPoints$depth <- rep(layerNames[[i]], times = nrow(tempPoints))
  l_chalumnae_down <- rbind(l_chalumnae_down, tempPoints)
}

# creating spatial points for downsampled occs
l_chalumnae_down_sp <- l_chalumnae_down
coordinates(l_chalumnae_down_sp) <- ~longitude + latitude
l_chalumnae_down_sp <- st_as_sf(l_chalumnae_down_sp)
st_crs(l_chalumnae_down_sp) <- target_crs

# occurrence plot
plot(l_chalumnae_down_sp$geometry, col = "red")
plot(land, col = NA, add = T)

pdf('./Plots/Coelacanth_3D_downsampled_occs_with_depth.pdf')
ggplot(data = land) +
  geom_sf() +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(24.875, 47.125), ylim = c(-37.125, -3.125)) +
  geom_point(data = l_chalumnae_down, aes(x = longitude, y = latitude, 
                                          color = depth), size = 3.5) +
  scale_colour_paletteer_c("viridis::plasma", trans = "reverse") +
  theme(panel.background = element_rect(fill = "lightblue")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  labs(colour = "Depth (m)")
dev.off()

# writing downsampled occurrences with indexed depth to file
write.csv(l_chalumnae_down, file = './data/l_chalumnae_down_3D.csv')

# 2D occurrences
# filtering for only Latimeria chalumnae
l_chalumnae <- all_occs %>% filter(Species == "Latimeria chalumnae")

# removing coordinates with large georeferencing error
l_chalumnae <- 
  l_chalumnae[-which(as.numeric(l_chalumnae$Georeference.error..m.) > 10000),]
l_chalumnae_maxent <- data.frame(l_chalumnae$Species, 
                                 as.numeric(l_chalumnae$Longitude),
                                 as.numeric(l_chalumnae$Latitude))
colnames(l_chalumnae_maxent) <- c("species", "longitude", "latitude")

# spatially thinning coordinates
l_chalumnae_maxent_final <- thin(loc.data = l_chalumnae_maxent,
                                 lat.col = "latitude",
                                 long.col = "longitude",
                                 spec.col = "species",
                                 thin.par = 5,
                                 reps = 1, write.files = F, 
                                 locs.thinned.list.return = T)
species <- rep("Latimeria chalumnae", times = nrow(l_chalumnae_maxent_final[[1]]))
l_chalumnae_maxent_final <- data.frame(species, l_chalumnae_maxent_final[[1]])
write.csv(l_chalumnae_maxent_final, file = './data/l_chalumnae_maxent_2D.csv')

