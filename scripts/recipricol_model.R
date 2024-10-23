# running pipeline to produce a reciprocal model

source('./scripts/all_functions_3D_sdm.R')

library(terra)
library(dplyr)
library(maxnet)
library(rnaturalearthhires)
library(predicts)
library(ENMeval)
library(voluModel)
library(sf)
library(usdm)
library(doParallel)
library(ggplot2)
library(tidyterra)
library(viridis)

# create land for plotting, create target CRS
target_crs <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"

land <- st_as_sf(rnaturalearthhires::countries10)[1] %>% st_set_crs(target_crs)
cea_crs <- "+proj=cea +lat_ts=0 +lon_0"
land_cea <- st_transform(land, crs = cea_crs)

# Read in cleaned L menadoensis dataset
l_menadoensis <- read.csv('./data/l_menadoensis_3D_occs.csv')
l_menadoensis <- l_menadoensis[,2:4]

# read in environmental layers, unmasked
# I need to pull the environmental layers unmasked for the presences but masked
# for the absences- many of the presences are excluded with the mask because
# the mask polygon is too coarse
conductivity_winter_indo <- 
  rast('./Envs_Indo/env_indo_3D_conductivity_winter_interpolated.tif')
conductivity_summer_indo <- rast('./Envs_Indo/conductivity_summer_18_indo_interp.tif')
current_north_december_indo <- 
  rast('./Envs_Indo/env_indo_3D_current_north_december_interpolated.tif')
current_north_june_indo <- rast('./Envs_Indo/current_north_june_indo_interp.tif')
current_east_december_indo <- rast('./Envs_Indo/current_east_december_indo_interp.tif')
current_east_june_indo <- rast('./Envs_Indo/current_east_june_indo_interp.tif')
slope_indo <- rast('./Envs_Indo/env_indo_3D_slope_interpolated.tif')
temperature_summer_indo <- rast('./Envs_Indo/env_indo_3D_temperature_summer_interpolated.tif')
temperature_winter_indo <- rast('./Envs_Indo/temperature_winter_18_indo_interp.tif')
salinity_winter_indo <- rast('./Envs_Indo/salinity_winter_18_indo_interp.tif')
salinity_summer_indo <- rast('./Envs_Indo/salinity_summer_18_indo_interp.tif')
density_winter_indo <- rast('./Envs_Indo/density_winter_18_indo_interp.tif')
density_summer_indo <- rast('./Envs_Indo/density_summer_18_indo_interp.tif')

# continents will be masked out when creating the accessible area

# delineate L menadoensis accessible area
# I'm going to test different buffer proportions on the 10th layer, 
# as that would be the only one with enough points
wanted_depths <- c(125, 150, 175)
wanted_test_points <- l_menadoensis[which(l_menadoensis$depth %in% wanted_depths),]
wanted_test_points <- wanted_test_points[,1:2]

# creating various accessible areas of buffer distances 4, 6, 8, and 10 times 
# larger than the area of the coordinates, and then running initial models to 
# select the best buffer for each depth slice via AIC
buff_props <- c(4, 6, 8, 10)

# creating polygon and transforming to cea to deal in meters
temp_alph <- marineBackground(wanted_test_points, fraction = 1, 
                              partCount = 2, 
                              clipToOcean = T,
                              buff = 1200000)
temp_alph <- st_transform(st_as_sf(temp_alph), crs = cea_crs)

# creating a list of shapes for each buffer proportion
shape_list <- vector("list", length = length(buff_props))
for (i in seq_along(buff_props)) {
  buffDist = (sqrt(buff_props[i]*st_area(temp_alph)) - sqrt(st_area(temp_alph)))/2
  shape_new <- st_buffer(x = temp_alph, dist = buffDist)
  shape_new_clipped <- st_difference(x = shape_new, y = st_union(land_cea))
  shape_final <- st_transform(shape_new_clipped, target_crs)
  shape_list[[i]] <- shape_final
}

# j loop to crop rasters to each polygon, run ENMeval for each polygon, and store
# model results for each polygon
# only running with interpolated rasters
wanted_env <- c(temperature_summer_indo[[23]], conductivity_winter_indo[[23]], 
                   current_north_december_indo[[23]], slope_indo[[23]])

model_results <- vector(length = length(shape_list))
for (i in 1:length(shape_list)) {
  cropped_env <- crop(wanted_env, shape_list[[i]])
  clipped_env <- mask(cropped_env, shape_list[[i]])
  simple_model <- ENMevaluate(occs = wanted_test_points, envs = clipped_env, 
                              algorithm = "maxnet",
                              tune.args = list(fc = c("L"), rm = 1),
                              partitions = "none", n.bg = 10000)
  mod_res <- min(eval.results(simple_model)$AICc)
  model_results[i] <- mod_res
}

# find which polygon produced the lowest AIC
want_shape_num <- which(model_results == min(model_results, na.rm = T))

# for loop to create accessible areas per depth slice
wanted_depths <- unique(l_menadoensis$depth)[order(unique(l_menadoensis$depth))]
depth_acc_area_list <- vector("list", length = length(wanted_depths))
for (i in seq_along(wanted_depths)) {
  
  # pull all occurrences for a given depth slice. for each depth slice, we want
  # occurrences from the slice above and the slice below in order to create 
  # connectivity
  if(i == 1) {
    needed_slices <- c(wanted_depths[i], wanted_depths[i+1])
    new_df1 <- l_menadoensis %>% filter(depth == needed_slices[1])
    new_df2 <- l_menadoensis %>% filter(depth == needed_slices[2])
    new_df <- rbind(new_df1, new_df2)
  } else if(i == length(wanted_depths)) {
    needed_slices <- c(wanted_depths[i-1], wanted_depths[i])
    new_df1 <- l_menadoensis %>% filter(depth == needed_slices[1])
    new_df2 <- l_menadoensis %>% filter(depth == needed_slices[2])
    new_df <- rbind(new_df1, new_df2)
  } else {
    needed_slices <- c(wanted_depths[i-1], wanted_depths[i], wanted_depths[i+1])
    new_df1 <- l_menadoensis %>% filter(depth == needed_slices[1])
    new_df2 <- l_menadoensis %>% filter(depth == needed_slices[2])
    new_df3 <- l_menadoensis %>% filter(depth == needed_slices[3])
    new_df <- rbind(new_df1, new_df2, new_df3)
  }
  new_df_sp <- unique(new_df[,1:2])
  
  # creating polygon and transforming to cea to deal in meters
  temp_alph <- marineBackground(new_df_sp, fraction = 1, 
                                partCount = 2, 
                                clipToOcean = T,
                                buff = 1200000)
  temp_alph <- st_transform(st_as_sf(temp_alph), crs = cea_crs)
  
  # buffering shape by a fact of 4 times the area covered by the points
  buffDist = (sqrt(buff_props[want_shape_num]*st_area(temp_alph)) - 
                sqrt(st_area(temp_alph)))/2
  shape_new <- st_buffer(x = temp_alph, dist = buffDist)
  shape_new_clipped <- st_difference(x = shape_new, y = st_union(land_cea))
  shape_final <- st_transform(shape_new_clipped, target_crs)
  depth_acc_area_list[[i]] <- shape_final
}

indo_acc_area <- do.call(rbind, depth_acc_area_list)

st_write(indo_acc_area, dsn = './Accessible_area', 
         layer = '3D_indo_acc_area',
         driver = "ESRI Shapefile",
         append = F)

# cropping and interpolating environmental layers
# creating depth slice list to match between envs and acc areas
envs_list <- list(conductivity_summer_indo, conductivity_winter_indo, temperature_summer_indo,
                  temperature_winter_indo, density_summer_indo, density_winter_indo,
                  salinity_summer_indo, salinity_winter_indo, current_east_december_indo,
                  current_east_june_indo, current_north_december_indo, 
                  current_north_june_indo, slope_indo)

depths <- unique(l_menadoensis$depth)[order(unique(l_menadoensis$depth))]
env_depths <- as.numeric(gsub("X", '', names(envs_list[[1]])))
acc_area_call <- vector("list", length = length(depths))
for(i in 1:length(depths)) {
  if(i == 1) {
    t <- env_depths[which(env_depths <= depths[i])]
    acc_area_call[[i]] <- rep(i, times = length(t))
  }  else if(i == length(depths)) {
    t <- env_depths[which(env_depths > depths[i-1])]
    acc_area_call[[i]] <- rep(i, times = length(t))
  } else {
    t <- env_depths[which(env_depths <= depths[i] & env_depths > depths[i-1])]
    acc_area_call[[i]] <- rep(i, times = length(t))
  }
}
acc_area_call <- unlist(acc_area_call)

envs_new <- vector("list", length = length(envs_list))
for (i in 1:length(envs_list)) {
  envs_new[[i]] <- crop(envs_list[[i]], indo_acc_area[6,])
}

envs_names <- c("conductivity_summer", "conductivity_winter", "temperature_summer", 
                "temperature_winter", "density_summer", "density_winter", "salinity_summer",
                "salinity_winter", "current_east_december", "current_east_june", 
                "current_north_december", "current_north_june", "slope")

# converting list of envs into a list of spatraster stacks, 
# so each list element is a depth layer
envs_trans <- env_stack_transform(envs_new, envs_names = envs_names)

# using crop_stack function to mask each depth layer to the associated accessible
# area
envs_final <- vector("list", length = length(envs_trans))
for(i in 1:length(envs_trans)) {
  envs_final[[i]] <- crop_stack(enviro = envs_trans[[i]], 
                                accarea = indo_acc_area[acc_area_call[i],], 
                                which_interpolate = NA)
}

# creating list of spatial points by depth slice
colnames(l_menadoensis) <- c("longitude", "latitude", "depth")
l_menadoensis_sp_list <- xyzmat_to_3Dsp(depth_slices = env_depths,
                                      occ_mat = l_menadoensis)

# creating list of background points by depth slice
l_menadoensis_bg_list <- bg_list_maker(enviro_stack = envs_final,
                                     depth_slices = env_depths,
                                     wanted_var = "temperature_summer",
                                     wanted_num = 1000,
                                     presences = l_menadoensis_sp_list)

# selecting variables with highest permutation importance and which reduce VIF
# I'm using the unmasked brick to select variables and run the model
# pulling the background from the masked brick should already ensure that all background 
# points are not on land and are withing the accessible area, running it unmasked ensures
# all points have information
envs_for_model <- var_select(wanted_brick = envs_trans,
                             wanted_sp = l_menadoensis_sp_list,
                             wanted_bg = l_menadoensis_bg_list)

# producing 3D models
# extracting values at occs and bg
all_l_menadoensis_df <- brick_extract(wanted_brick = envs_for_model, 
                                    wanted_sp = l_menadoensis_sp_list,
                                    wanted_bg = l_menadoensis_bg_list)

# defining partition schemes
l_menadoensis_k <- partition_3D(l_menadoensis_sp_list, l_menadoensis_bg_list, 
                              all_l_menadoensis_df[[3]],
                              'k.fold',
                              kfolds = 5)

# training and testing models
l_menadoensis_mod_k <- maxent_3D(df_for_maxent = all_l_menadoensis_df[[3]],
                               wanted_fc = c("L", "Q", "LQ", "LQH"),
                               wanted_rm = c(1:4),
                               wanted_partition = l_menadoensis_k,
                               projection_layers = envs_for_model,
                               occ_list = l_menadoensis_sp_list)
want_num <- which(l_menadoensis_mod_k$results$avg.delta.AICc == 
        min(l_menadoensis_mod_k$results$avg.delta.AICc))
writeRaster(l_menadoensis_mod_k$predictions[[want_num]], 
            filename = './Model_results/menadoensis_3D_suitability.tif',
            overwrite = T)
write.csv(l_menadoensis_mod_k$results, 
          file = './Model_results/menadoensis_3D_results.csv')
best_mod <- l_menadoensis_mod_k$models[[want_num]]
save(best_mod, 
     file = './Model_results/menadoensis_3D_best_mod.Rdata')

# thresholding models
l_menadoensis_suit <- l_menadoensis_mod_k$predictions[[want_num]]
l_menadoensis_k_thresh <- threshold_3D(predicted_layers = l_menadoensis_suit,
                                     thresholding_vals = c(0.99, 0.95, 0.9, 0.85, 0.8, 0.75), 
                                     occ_list = l_menadoensis_sp_list,
                                     bg_list = l_menadoensis_bg_list)
writeRaster(l_menadoensis_k_thresh$threshold_layers,
            filename = './Model_results/menadoensis_3D_thresholded.tif',
            overwrite = T)
write.csv(l_menadoensis_k_thresh$tss_results, 
          file = './Model_results/menadoensis_3D_tss.csv')

# plotting suitability all layers
# suitability plot layer by layer
men_depths <- gsub('X', '', names(conductivity_summer_indo))

l_men_suit <- rast('./Model_results/menadoensis_3D_suitability.tif')
l_men <- rast('./Model_results/menadoensis_3D_thresholded.tif')

pdf('./Plots/reciprocal_indo_3D_suit_all_layers.pdf', pointsize = 3)
for(i in 1:dim(l_men_suit)[3]) {
  p <- ggplot() +
    geom_spatraster(data = l_men_suit[[i]]) +
    geom_sf(data = land) +
    coord_sf(xlim = c(100, 160), ylim = c(-21, 26), expand = F) +
    scale_fill_gradientn(colors = c("white", viridis(4)), 
                         values = c(0, 0.25, 0.5, 0.75, 1),
                         na.value = "transparent",
                         name = "Suitability Score") +
    theme(panel.grid.major = element_blank(), 
          panel.background = element_rect(fill = "white")) +
    ggtitle(paste0("Suitability at ", men_depths[i], " m"))
  print(p)
}
dev.off()

# thresholded plot layer by layer
pdf('./Plots/reciprocal_indo_3D_thresh_all_layers.pdf', pointsize = 3)
for(i in 1:dim(l_men)[3]) {
  p <- ggplot() +
    geom_spatraster(data = l_men[[i]]) +
    geom_sf(data = land) +
    coord_sf(xlim = c(100, 160), ylim = c(-21, 26), expand = F) +
    scale_fill_gradientn(colors = c("white", "darkgreen"), 
                         values = c(0, 1),
                         na.value = "transparent") +
    theme(panel.grid.major = element_blank(), 
          panel.background = element_rect(fill = "white"), legend.position = "none") +
    ggtitle(paste0("Thresholded Presence at ", men_depths[i], " m"))
  print(p)
}
dev.off()

