# Coelacanth 3D modelling

source('./scripts/all_functions_3D_sdm.R')

library(raster)
library(dplyr)
library(maxnet)
library(rnaturalearthhires)
library(dismo)
library(ENMeval)
library(voluModel)
library(sf)
library(usdm)
library(doParallel)

# create land for plotting, create target CRS
target_crs <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"

land <- st_as_sf(rnaturalearthhires::countries10)[1] %>% st_set_crs(target_crs)
land_poly <- as_Spatial(land)

# reading in occurrences
l_chalumnae_down <- read.csv('./data/l_chalumnae_down_3D.csv')

# reading in environmental rasters
density_summer_18 <- rast('./WOA_18/density_summer_18.tif')
density_winter_18 <- rast('./WOA_18/density_winter_18.tif')
conductivity_summer_18 <- rast('./WOA_18/conductivity_summer_18.tif')
conductivity_winter_18 <- rast('./WOA_18/conductivity_winter_18.tif')
salinity_summer_18 <- rast('./WOA_18/salinity_summer_18.tif')
salinity_winter_18 <- rast('./WOA_18/salinity_winter_18.tif')
temperature_summer_18 <- rast('./WOA_18/temperature_summer_18.tif')
temperature_winter_18 <- rast('./WOA_18/temperature_winter_18.tif')
current_east_december <- rast('./Currents/current_east_december.tif')
current_east_june <- rast('./Currents/current_east_june.tif')
current_north_december <- rast('./Currents/current_north_december.tif')
current_north_june <- rast('./Currents/current_north_june.tif')
slope <- rast('./Slope/slope.tif')

envs_list <- list(density_summer_18, density_winter_18, conductivity_summer_18,
                  conductivity_winter_18, salinity_summer_18,
                  salinity_winter_18, temperature_summer_18, 
                  temperature_winter_18, current_east_december, current_east_june,
                  current_north_december, current_north_june, slope)
envs_names <- c("density_summer", "density_winter", "conductivity_summer", 
                "conductivity_winter", "salinity_summer",
                "salinity_winter", "temperature_summer", "temperature_winter",
                "current_east_december", "current_east_june", 
                "current_north_december", "current_north_june", "slope")

# read in accessible area 
acc_area <- st_read(dsn = './Accessible_area', layer = '3D_coelacanth_acc_area')

# converting envs_list into list of stacks where each element is a depth layer
envs_new <- env_stack_transform(envs_list, envs_names)

# cropping environmental variables to the accessible area on a per depth slice
# basis

# creating depth slice list to match between envs and acc areas
depths <- unique(l_chalumnae_down$depth)
env_depths <- as.numeric(gsub("X", '', names(envs_list[[1]])))
acc_area_call <- vector("list", length = length(depths))
for(i in 1:length(depths)) {
  if(i == 1) {
    t <- env_depths[which(env_depths <= depths[i])]
    acc_area_call[[i]] <- rep(i, times = length(t))
  }  else if(i == 21) {
    t <- env_depths[which(env_depths > depths[i-1])]
    acc_area_call[[i]] <- rep(i, times = length(t))
  } else {
  t <- env_depths[which(env_depths <= depths[i] & env_depths > depths[i-1])]
  acc_area_call[[i]] <- rep(i, times = length(t))
  }
}
acc_area_call <- unlist(acc_area_call)

# using crop_stack function to mask each depth layer to the associated accessible
# area
envs_final <- vector("list", length = length(envs_new))
for(i in 1:length(envs_new)) {
  envs_final[[i]] <- crop_stack(enviro = envs_new[[i]], 
                                accarea = acc_area[acc_area_call[i],], 
                                which_interpolate = c(1:12))
}

# creating list of spatial points by depth slice
l_chalumnae_sp_list <- xyzmat_to_3Dsp(depth_slices = env_depths,
                                      occ_mat = l_chalumnae_down[,2:4])

# creating list of background points by depth slice
l_chalumnae_bg_list <- bg_list_maker(enviro_stack = envs_final,
                                     depth_slices = env_depths,
                                     wanted_var = "temperature_summer",
                                     wanted_num = 1000,
                                     presences = l_chalumnae_sp_list)

# selecting variables with highest permutation importance and which reduce VIF
envs_for_model <- var_select(wanted_brick = envs_final,
                             wanted_sp = l_chalumnae_sp_list,
                             wanted_bg = l_chalumnae_bg_list)

# producing 3D models
# extracting values at occs and bg
all_l_chalumnae_df <- brick_extract(wanted_brick = envs_for_model, 
                                    wanted_sp = l_chalumnae_sp_list,
                                    wanted_bg = l_chalumnae_bg_list)

# defining partition schemes
l_chalumnae_k <- partition_3D(l_chalumnae_sp_list, l_chalumnae_bg_list, 
                                 all_l_chalumnae_df[[3]],
                                 'k.fold',
                                 kfolds = 5)

l_chalumnae_block <- partition_3D(l_chalumnae_sp_list, l_chalumnae_bg_list, 
                                  all_l_chalumnae_df[[3]],
                                  'block',
                                  orientation = 'lat_lon')

# training and testing models
l_chalumnae_mod_k <- maxent_3D(df_for_maxent = all_l_chalumnae_df[[3]],
                               wanted_fc = c("L", "Q", "LQ", "LQH"),
                               wanted_rm = c(1:4),
                               wanted_partition = l_chalumnae_k,
                               projection_layers = envs_for_model,
                               occ_list = l_chalumnae_sp_list)
writeRaster(l_chalumnae_mod_k$predictions[[14]], 
            filename = './Model_results/coelacanth_3D_suitability.tif',
            overwrite = T)
write.csv(l_chalumnae_mod_k$results, 
          file = './Model_results/coelacanth_3D_results.csv')
best_mod <- l_chalumnae_mod_k$models[[14]]
save(best_mod, 
     file = './Model_results/coelacanth_3D_best_mod.Rdata')

# thresholding models
l_chalumnae_suit <- l_chalumnae_mod_k$predictions[[14]]
l_chalumnae_k_thresh <- threshold_3D(projection_layers = NA, 
                        predicted_layers = l_chalumnae_suit,
                        thresholding_vals = c(0.99, 0.95, 0.90), 
                        occ_list = l_chalumnae_sp_list,
                        bg_list = l_chalumnae_bg_list)
writeRaster(l_chalumnae_k_thresh$threshold_layers,
            filename = './Model_results/coelacanth_3D_thresholded.tif',
            overwrite = T)
write.csv(l_chalumnae_k_thresh$tss_results, 
          file = './Model_results/coelacanth_3D_tss.csv')


save(all_l_chalumnae_df, file = './data/all_l_chalumnae_df.Rdata')
save(envs_for_model, file = './data/envs_for_model.Rdata')
save(l_chalumnae_sp_list, file = './data/l_chalumnae_sp_list.Rdata')
save(l_chalumnae_bg_list, file = './data/l_chalumnae_bg_list.Rdata')

