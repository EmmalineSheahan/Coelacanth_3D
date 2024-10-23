# projecting the l menadoensis model onto the Mozambique channel

source('./scripts/all_functions_3D_sdm.R')

library(terra)
library(sf)
library(ggplot2)
library(dplyr)
library(voluModel)
library(rnaturalearthhires)
library(predicts)
library(tidyterra)
library(viridis)
library(ds4psy)
library(gridExtra)

# creating target crs
target_crs <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"
land <- st_as_sf(rnaturalearthhires::countries10)[1] %>% st_set_crs(target_crs)

# Reading in L chalumnae occurrences and environmental variables
# this used different variables, so you will need to re-interpolate
l_chalumnae_occ <- read.csv('./data/l_chalumnae_down_3D.csv')
current_east_june <- rast('./Currents/current_east_june.tif')
slope <- rast('./WOA_18/env_3D_slope_interpolated.tif')
temperature_winter <- rast('./WOA_18/temperature_winter_18.tif')
conductivity_summer <- rast('./WOA_18/conductivity_summer_18.tif')

# interpolating rasters in menadoensis not in chalumnae for projecting
for (i in 1:dim(current_east_june)[3]) {
  current_east_june[[i]] <- interpolateRaster(current_east_june[[i]])
}
current_east_june <- mask(current_east_june, land, inverse = T)

for (i in 1:dim(temperature_winter)[3]) {
  temperature_winter[[i]] <- interpolateRaster(temperature_winter[[i]])
}
temperature_winter <- mask(temperature_winter, land, inverse = T)

for (i in 1:dim(conductivity_summer)[3]) {
  conductivity_summer[[i]] <- interpolateRaster(conductivity_summer[[i]])
}
conductivity_summer <- mask(conductivity_summer, land, inverse = T)

# creating list of spatraster stacks
env_3D <- list(conductivity_summer, temperature_winter, current_east_june, slope)

# transforming list format
env_3D <- env_stack_transform(envs_all = env_3D,
                              envs_names = c("conductivity_summer",
                                             "temperature_winter",
                                             "current_east_june",
                                             "slope"))

# reading in model
load('./Model_results/menadoensis_3D_best_mod.Rdata')

# projecting model
mozam_suit_list_3D <- vector("list", length = length(env_3D))
for (i in 1:length(env_3D)) {
  mozam_suit_list_3D[[i]] <- predict(best_mod, env_3D[[i]])
}
mozam_suit_3D <- rast(mozam_suit_list_3D)
writeRaster(mozam_suit_3D, 
            filename = './Model_results/Mozambique_projection_3D.tif',
            overwrite = T)

# plotting menadoensis projection onto mozambique channel
men_depths <- as.numeric(gsub('X', '', names(temperature_winter)))

pdf('./Plots/menadoensis_in_mozambique_3D_suit_all_layers.pdf', pointsize = 3)
for(i in 1:dim(mozam_suit_3D)[3]) {
  p <- ggplot() +
    geom_spatraster(data = mozam_suit_3D[[i]]) +
    geom_sf(data = land) +
    coord_sf(xlim = c(9.25, 63.75), ylim = c(-48.75, 14.25), expand = F) +
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

# assessing how many L chalumnae points are successfully predicted

# reading in suitability values at threshold levels
menadoensis_tss <- read.csv('./Model_results/menadoensis_3D_tss.csv')
sens <- menadoensis_tss$Sensitivity
suit <- menadoensis_tss$Suitability

# thresholding model
mozam_thresh_3D_list <- vector("list", length = length(sens))
for (i in 1:length(sens)) {
  all_layers <- vector("list", length = dim(mozam_suit_3D)[3])
  for (j in 1:dim(mozam_suit_3D)[3]) {
    rclmat <- matrix(data = c(0, suit[i], 0, 
                              suit[i], 1, 1), nrow = 2, ncol = 3, byrow = T)
    all_layers[[j]] <- classify(mozam_suit_3D[[j]], 
                                  rcl = rclmat)
  }
  mozam_thresh_3D_list[[i]] <- rast(all_layers)
}
writeRaster(mozam_thresh_3D_list[[5]], 
            filename = './Model_results/Mozambique_thresholded_3D.tif',
            overwrite = T)

# creating list of spatial points for chalumnae
# reading in occurrences
l_chalumnae_3D <- l_chalumnae_occ[,2:4]
men_depths <- as.numeric(gsub('X', '', names(temperature_winter)))
l_chalumnae_3D_sp <- xyzmat_to_3Dsp(depth_slices = men_depths, 
                                    occ_mat = l_chalumnae_3D)

# extracting values at occurrence records for each threshold
thresh_val_list <- vector("list", length = length(sens))
omission_rate <- vector(length = length(sens))
for (i in 1:length(sens)) {
  extract_list <- vector("list", length = dim(mozam_thresh_3D_list[[i]])[3])
  for (j in 1:dim(mozam_thresh_3D_list[[i]])[3]) {
    if(any(is.na(l_chalumnae_3D_sp[[j]]))) {
      extract_list[[j]] <- NA
    } else {
      extract_list[[j]] <- terra::extract(x = mozam_thresh_3D_list[[i]][[j]],
                                          y = l_chalumnae_3D_sp[[j]])
    }
  }
  thresh_val_list[[i]] <- do.call(rbind, extract_list)
  omission_rate[i] <- (1 - 
                         (sum(thresh_val_list[[i]][,2], na.rm = T)/nrow(l_chalumnae_3D)))
}
mozam_omission <- cbind(menadoensis_tss, omission_rate)
colnames(mozam_omission) <- c("Sensitivity", "Specificity", "TSS", "Suitability", 
                              "Omission Rate")
write.csv(mozam_omission, file = "./Model_results/mozam_omission.csv")
