# Formatting and thinning the full Coelacanth occurrence data set

library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(voluModel)
library(rnaturalearth)
library(rgdal)

# creating wanted crs and land polygon
wanted_crs <- make_EPSG() %>% filter(code == 4326)
target_crs <- wanted_crs$prj4

land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]

# reading in the dataset
all_occs <- 
  read.csv('./data/coelacanth_complete_occurrence.csv')

# filtering for only Latimeria chalumnae
l_chalumnae <- all_occs %>% filter(Species == "Latimeria chalumnae")

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

# reading in bathymetry layer
bath <- rast('./coelacanth_data/bathymetry/ETOPO1_Bed_c_geotiff.tif')
bath_crop <- crop(bath, l_chalumnae_sp)

bath_depths <- extract(bath_crop, l_chalumnae_sp)

depth_compare <- cbind(l_chalumnae_occs, abs(bath_depths))
which(depth_compare$depth == depth_compare$ETOPO1_Bed_c_geotiff)

# reading in a WOA 1 degree raster layer for spatial thinning

# thinning occurrences
l_chalumnae_thin <- 
  l_chalumnae %>% dplyr::select(c(Species, Longitude, Latitude))


