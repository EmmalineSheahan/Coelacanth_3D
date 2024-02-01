# projecting the best 2D and 3D models onto the Indo-Pacific to recover omission error
# of Latimeria menadoensis

library(dplyr)
library(terra)
library(sf)
library(dismo)
library(raster)
library(rgdal)
library(doParallel)
library(rnaturalearthhires)

# creating target crs
wanted_crs <- make_EPSG() %>% filter(code == 4326)
target_crs <- wanted_crs$prj4

land <- st_as_sf(rnaturalearthhires::countries10)[1] %>% st_set_crs(target_crs)

# 2D model
# reading in necessary components
l_menadoensis_2D <- read.csv('./data/l_menadoensis_2D_occs.csv')

cl <- makeCluster(2)
registerDoParallel(cl, cores = 2)

indo_env_2D <- rast('./Envs_Indo/all_2D_env_indo.tif')
for (i in 1:dim(indo_env_2D)[3]) {
  indo_env_2D[[i]] <- interpolateRaster(indo_env_2D[[i]])
}

stopCluster(cl)

indo_env_2D_masked <- mask(indo_env_2D, land, inverse = T)

writeRaster(rast(indo_env_2D), filename = './data/temporary_2D_indo_env_real.tif',
            overwrite = T)

indo_env_2D <- stack(indo_env_2D)
indo_env_2D_masked <- stack(indo_env_2D_masked)

tss_2D <- read.csv('./Model_results/l_chalumnae_2D_sdm_tss_results.csv')

load('./Model_results/best_sdm_mod_2D.Rdata')

# projecting model
Indo_pacific_pred <- dismo::predict(best_sdm_model, indo_env_2D_masked)
writeRaster(Indo_pacific_pred, 
            filename = './Model_results/indopacific_projection_2D.tif',
            overwrite = T)

# creating thresholded rasters according to L chalumnae
indo_2D_thresh_list <- vector("list", length = length(tss_2D$Sensitivity))
for (i in 1:length(tss_2D$Sensitivity)) {
  indo_2D_thresh_list[[i]] <- reclassify(Indo_pacific_pred, 
                   rcl = c(0, tss_2D$Suitability[i], 0, tss_2D$Suitability[i], 1, 1))
}
writeRaster(indo_2D_thresh_list[[2]], 
            filename = './Model_results/indopacific_threshold_2D.tif')

# creating spatial points in order to extract values
l_menadoensis_2Dsp <- l_menadoensis_2D
coordinates(l_menadoensis_2Dsp) <- ~Longitude+Latitude
proj4string(l_menadoensis_2Dsp) <- target_crs
omm_df_list <- vector("list", length = length(indo_2D_thresh_list))
ommision_rates_2D <- vector(length = 3)
for (i in 1:length(indo_2D_thresh_list)) {
  omm_df_list[[i]] <- extract(indo_2D_thresh_list[[i]], l_menadoensis_2Dsp)
  ommision_rates_2D[i] <- (1 - 
                             sum(omm_df_list[[i]], na.rm = T)/length(omm_df_list[[i]]))
}

tss_2D_new <- cbind(tss_2D, ommision_rates_2D)
write.csv(tss_2D_new, 
          file = './Model_results/coelacanth_2D_tss_with_menaoensis_om.csv')

# 3D model
# reading in environmental variables to match model
conductivity_winter <- rast('./Envs_Indo/conductivity_winter_18_indo.tif')
temperature_summer <- rast('./Envs_Indo/temperature_summer_18_indo.tif')
current_north_december <- rast('./Envs_Indo/current_north_december_indo.tif')
slope <- rast('./Envs_Indo/slope_indo.tif')

env_indo_3D <- list(conductivity_winter, temperature_summer, current_north_december,
                    slope)

cl <- makeCluster(6)
registerDoParallel(cl, cores = 6)

for (i in 1:length(env_indo_3D)) {
  for (j in 1:dim(env_indo_3D[[i]])[3]) {
    env_indo_3D[[i]][[j]] <- interpolateRaster(env_indo_3D[[i]][[j]])
  }
  env_indo_3D[[i]] <- mask(env_indo_3D[[i]], land, inverse = T)
}

stopCluster(cl)

# transforming list format
env_indo_3D <- env_stack_transform(envs_all = env_indo_3D,
                                   envs_names = c("conductivity_winter", 
                                                  "temperature_summer",
                                                  "current_north_december",
                                                  "slope"))

# loading in model
load('./Model_results/coelacanth_3D_best_mod.Rdata')

# reading in tss df
tss_3D <- read.csv('./Model_results/coelacanth_3D_tss.csv')

# projecting model
indo_suit_list_3D <- vector("list", length = length(env_indo_3D))
for (i in 1:length(env_indo_3D)) {
  Indo_suit_list_3D <- dismo::predict(best_mod, indo_env_3D[[i]])
}
indo_suit_3D <- rast(indo_suit_list_3D)
writeRaster(indo_suit_3D, filename = './Model_results/Indopacific_projection_3D.tif',
            overwrite = T)

# thresholding model
indo_3D_thresh_list <- vector("list", length = length(tss_3D$Sensitivity))
for (i in 1:length(tss_3D$Sensitivity)) {
  all_layers <- vector("list", length = dim(indo_suit_3D)[3])
  for (j in 1:dim(indo_suit_3D)[3]) {
    all_layers[[j]] <- reclassify(Indo_suit_3D[[j]], 
              rcl = c(0, tss_3D$Suitability[i], 0, tss_3D$Suitability[i], 1, 1))
  }
  indo_3D_thresh_list[[i]] <- rast(all_layers)
}
witeRaster(indo_3D_thresh_list[[2]], 
           filename = './Model_results/Indopacific_thresholded_3D.tif')

# creating list of spatial points for menadoensis
# reading in occurrences
l_menadoensis_3D <- read.csv('./data/l_menadoensis_3D_occs.csv')
colnames(l_menadoensis_3D) <- c("longitude", "latitude", "depth")
men_depths <- names(conductivity_winter)
l_menadoensis_3D_sp <- xyzmat_to_3Dsp(depth_slices = men_depths, 
                                      occ_mat = l_menadoensis_3D)

# extracting values at occurrence records for each threshold
thresh_val_list <- vector("list", length = length(tss_3D$Sensitivity))
omission_rate <- vector(length = length(tss_3D$Sensitivity))
for (i in 1:length(tss_3D$Sensitivity)) {
  extract_list <- vector("list", length = dim(indo_3D_thresh_list[[i]])[3])
  for (j in 1:dim(indo_3D_thresh_list[[i]])[3]) {
    if(any(is.na(l_menadoensis_3D_sp[[j]]))) {
      extract_list[[j]] <- NA
    } else {
      extract_list[[j]] <- terra::extract(x = rast(indo_3D_thresh_list[[i]]),
                                          y = l_menadoensis_3D_sp[[j]])
    }
  }
  thresh_val_list[[i]] <- do.call(rbind, extract_list)
  omission_rate[i] <- (1 - 
                (sum(thresh_val_list[[i]], na.rm = T)/length(thresh_val_list[[i]])))
}

all_omission <- data.frame(tss_3D$Percentiles, ommission_rates_2D, omission_rate)
colnames(all_omission) <- c("threshold.percentile", "2D.omission", "3D.omission")
write.csv(all_omission, file = './Model_results/omission_rate_2D_and_3D.csv')
