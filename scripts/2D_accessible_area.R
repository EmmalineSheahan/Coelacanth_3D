# 2D accessible area

library(sf)
library(sp)
library(raster)
library(terra)
library(voluModel)
library(rgdal)
library(dplyr)
library(rgeos)
library(rnaturalearthhires)
library(ENMeval)

# creating wanted crs and land polygon
wanted_crs <- make_EPSG() %>% filter(code == 4326)
target_crs <- wanted_crs$prj4
cea_crs <- "+proj=cea +lat_ts=0 +lon_0"

land <- st_as_sf(rnaturalearthhires::countries10)[1] %>% st_set_crs(target_crs)
land_poly <- as_Spatial(land)
land_poly_cea <- spTransform(land_poly, CRSobj = cea_crs)

coasts <- st_as_sf(rnaturalearthhires::coastline10)[1] %>% 
  st_set_crs(target_crs)
coast_poly <- as_Spatial(coasts)
coast_poly_cea <- spTransform(coast_poly, CRSobj = cea_crs)

# reading in occs
l_chalumnae_2D <- read.csv('./data/l_chalumnae_maxent_2D.csv')

# creating various accessible areas of buffer distances 2, 4, 6, and 8 times 
# larger than the area of the coordinates, and then running initial models to 
# select the best buffer for each depth slice via AIC
buff_props <- c(2, 4, 6, 8)

# reading in 2D envs
all_2D_env <- rast('./Envs_2D/all_2D_env.tif')
for (i in 1:dim(all_2D_env)[3]) {
  all_2D_env[[i]] <- interpolateRaster(all_2D_env[[i]])
}

# creating polygon and transforming to cea to deal in meters
temp_alph <- marineBackground(l_chalumnae_2D, fraction = 1, 
                              partCount = 2, 
                              clipToOcean = T,
                              buff = 600000)
temp_alph <- spTransform(as_Spatial(st_as_sf(temp_alph)), 
                         CRSobj = cea_crs)

# creating a list of shapes for each buffer proportion
shape_list <- vector("list", length = length(buff_props))
for (i in seq_along(buff_props)) {
  buffDist = (sqrt(buff_props[i]*gArea(temp_alph)) - sqrt(gArea(temp_alph)))/2
  shape_new <- raster::buffer(x = temp_alph, width = buffDist, dissolve = T)
  shape_new_clipped <- gDifference(shape_new, land_poly_cea)
  shape_sf <- st_as_sf(shape_new_clipped)
  shape_final <- st_transform(shape_sf, target_crs)
  shape_list[[i]] <- shape_final
}

# j loop to crop rasters to each polygon, run ENMeval for each polygon, and store
# model results for each polygon
model_results <- vector(length = length(shape_list))
for (i in 1:length(shape_list)) {
  cropped_env <- crop(all_2D_env, shape_list[[i]])
  clipped_env <- mask(cropped_env, shape_list[[i]])
  simple_model <- ENMevaluate(occs = l_chalumnae_2D[,3:4], envs = clipped_env, 
                              algorithm = "maxnet",
                              tune.args = list(fc = c("L"), rm = 1),
                              partitions = "none", n.bg = 10000)
  mod_res <- min(eval.results(simple_model)$AICc)
  model_results[i] <- mod_res
}

# find which polygon produced the lowest AIC
want_shape_num <- which(model_results == min(model_results, na.rm = T))

# mask envs to accessible area and write to file
all_2D_env_final <- mask(all_2D_env, shape_list[[want_shape_num]])
writeRaster(all_2D_env_final, filename = './Envs_2D/all_2D_env_final.tif',
                filetype = 'GTiff', overwrite = T)

for(i in 1:length(names(all_2D_env_final))) {
  values(all_2D_env_final[[i]])[which(is.na(values(all_2D_env_final[[i]])))] <- 0
  writeRaster(raster(all_2D_env_final[[i]]), filename = paste0('./Envs_2D/', 
              names(all_2D_env_final)[i], '_for_maxent.asc'),
              format = 'ascii',
              overwrite = T)
}
