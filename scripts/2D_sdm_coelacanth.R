# coelacanth 2D model production and validation

library(dplyr)
library(rnaturalearthhires)
library(sf)
library(terra)
library(ggplot2)
library(ENMeval)
library(dismo)
library(raster)
library(usdm)

# creating land for plotting
land <- st_as_sf(rnaturalearthhires::countries10)[1] %>% 
  st_set_crs("+proj=longlat +datum=WGS84")

# reading in occurrence data
l_chalumnae_2D <- read.csv('./data/l_chalumnae_maxent_2D.csv')

# reading in environmental rasters cropped and masked to accessible area
final_envs_2D <- rast('./Envs_2D/all_2D_env_final.tif')

# generating background points for modelling with partition schemes
l_chalumnae_bg <- randomPoints(raster(final_envs_2D[[1]]), n = 10000, 
                               p = l_chalumnae_2D[,3:4])
colnames(l_chalumnae_bg) <- colnames(l_chalumnae_2D[,3:4])

# var_remove function for single iterative variable removal
# wanted_envs = the clipped environmental spatraster used in the model
# wanted_occs = the two column occurrence matrix
# wanted_bg = the two column background matrix
# returns the spatraster with the least important variable discarded
var_remove <- function(wanted_envs, wanted_occs, wanted_bg) {
  rast_list <- vector("list", length = dim(wanted_envs)[3])
  for (j in 1:dim(wanted_envs)[3]) {
    env_new <- raster(wanted_envs[[j]])
    rast_list[[j]] <- env_new
  }
  initial_envs <- stack(rast_list)
  initial_mod <- maxent(initial_envs, p = wanted_occs, a = wanted_bg)
  initial_results <- data.frame(initial_mod@results)
  measure_name <- row.names(initial_results)
  perm_import <- data.frame(measure_name, 
                            initial_results[,1])
  colnames(perm_import) <- c("measure_name", "result")
  want_pull <- grep('.permutation.importance', perm_import$measure_name)
  final_perm_import <- perm_import[want_pull,]
  perm_remove <- which(final_perm_import$result == min(final_perm_import$result))
  if (length(perm_remove > 1)) {
    perm_remove <- perm_remove[1]
  }
  remove_this <- final_perm_import$measure_name[perm_remove]
  remove_this <- gsub('.permutation.importance', '', remove_this)
  
  # removing variable from analysis
  remove_env <- which(names(wanted_envs) %in% remove_this) 
  new_envs <- wanted_envs[[-remove_env]]
  
  return(new_envs)
}

# var_select function to iterate var_remove until all vifs are below 5
var_select <- function(wanted_envs, wanted_occs, wanted_bg) {
  
  final_env <- wanted_envs
  t <- max(usdm::vif(final_env)$VIF)
  
  while(t > 5) {
    final_env <- var_remove(wanted_envs = final_env, wanted_occs = wanted_occs,
                            wanted_bg = wanted_bg)
    if(dim(final_env)[3] < 3) {
      t <- 4
    } else {
      t <- max(usdm::vif(final_env)$VIF)
    }
  }
  
  return(final_env)
}

# removing variables of low permutation importance
new_envs_2D <- var_select(final_envs_2D, l_chalumnae_2D[,3:4], l_chalumnae_bg)

# creating k fold partitions
l_chalumnae_k_fold <- get.randomkfold(occs = l_chalumnae_2D[,3:4], 
                                      bg = l_chalumnae_bg, kfolds = 5)
user.grp_kfold <- list(occs.grp = l_chalumnae_k_fold$occs.grp, 
                       bg.grp = l_chalumnae_k_fold$bg.grp)

# creating block partitions
l_chalumnae_block <- get.block(occs = l_chalumnae_2D[,3:4], 
                               bg = l_chalumnae_bg)
user.grp_block <- list(occs.grp = l_chalumnae_block$occs.grp, 
                       bg.grp = l_chalumnae_block$bg.grp)

# running ENMeval with random kfold partition
l_chalumnae_mod_kfold <- ENMevaluate(occs = l_chalumnae_2D[,3:4], 
                                     envs = new_envs_2D,
                                     bg = l_chalumnae_bg, algorithm = "maxent.jar",
                                     tune.args = list(fc = c("L","LQ","LQH"), 
                                                      rm = 1:3),
                                     partitions = 'user', user.grp = user.grp_kfold)

mod_k_results <- eval.results(l_chalumnae_mod_kfold)
write.csv(mod_k_results, file = './Model_results/coelacanth_2D_mod_kfold_eval.csv')

# running ENMeval with block partitioning
l_chalumnae_mod_block <- ENMevaluate(occs = l_chalumnae_2D[,3:4], 
                                     envs = new_envs_2D,
                                     bg = l_chalumnae_bg, algorithm = "maxent.jar",
                                     tune.args = list(fc = c("L","LQ","LQH"), 
                                                      rm = 1:3),
                                     partitions = 'user', user.grp = user.grp_block)
mod_block_results <- eval.results(l_chalumnae_mod_block)
write.csv(mod_block_results, file = './Model_results/Coelacanth_2D_mod_block_eval.csv')

# block and kfold partitions return the same result
# selecting best model, thresholding, and writing to file
if (any(eval.results(l_chalumnae_mod_kfold)$auc.val.avg) < 0.7) {
  select_this <- which(eval.results(l_chalumnae_mod_kfold)$auc.val.avg == 
                         max(eval.results(l_chalumnae_mod_kfold)$auc.val.avg, 
                             na.rm = T))
  if(length(select_this) > 1) {
    select_this <- select_this[1]
  }
  best_sdm <- l_chalumnae_mod_kfold@predictions[[select_this]]
  best_sdm_tune <- l_chalumnae_mod_kfold@tune.settings[select_this,]
} else {
  select_this <- which(eval.results(l_chalumnae_mod_kfold)$AICc == 
                         min(eval.results(l_chalumnae_mod_kfold)$AICc, na.rm = T))
  if(length(select_this) > 1) {
    select_this <- select_this[1]
  }
  best_sdm_model <- l_chalumnae_mod_kfold@models[[select_this]]
  best_sdm <- l_chalumnae_mod_kfold@predictions[[select_this]]
  best_sdm_tune <- l_chalumnae_mod_kfold@tune.settings[select_this,]
}
write.table(c(best_sdm_tune, names(new_envs_2D)), 
            file = './Model_results/coelacanth_2D_params.txt',
            row.names = F, col.names = F)
writeRaster(best_sdm, filename = './Model_results/coelacanth_2D_suitability.tif',
            overwrite = T)
save(best_sdm_model, file = './Model_results/best_sdm_mod_2D.Rdata')

# thresholding
percentiles <- c(0.99, 0.95, 0.9)
presences <- l_chalumnae_2D[,3:4]
coordinates(presences) <- ~Longitude+Latitude
proj4string(presences) <- "+proj=longlat +datum=WGS84"
absences <- data.frame(l_chalumnae_bg)
coordinates(absences) <- ~Longitude+Latitude
proj4string(absences) <- "+proj=longlat +datum=WGS84"
sdm_suit <- extract(best_sdm, presences)[order(extract(best_sdm, presences), 
                                               na.last = NA,
                                               decreasing = T)]
thresh_sdm_list <- vector("list", length = length(percentiles))
tss_list <- vector(length = length(percentiles))
thresh_val_list <- vector(length = length(percentiles))
specificity_list <- vector(length = length(percentiles))
for (j in seq_along(percentiles)) {
  thresh <- round(length(sdm_suit)*percentiles[j])
  thresh_val <- sdm_suit[thresh]
  thresh_val_list[j] <- thresh_val
  thresholded_sdm <- reclassify(best_sdm, 
                                rcl = c(0, thresh_val, 0, thresh_val, 1, 1))
  thresh_sdm_list[[j]] <- thresholded_sdm
  real_abs <- extract(thresholded_sdm, absences)
  specificity <- length(which(real_abs == 0))/length(real_abs)
  specificity_list[j] <- specificity
  sensitivity <- percentiles[j]
  tss <- (sensitivity + ((1/3)*specificity)) - 1
  tss_list[j] <- tss
}

choose_final <- which(tss_list == max(tss_list))
final_sdm <- thresh_sdm_list[[choose_final]]

tss_df <- data.frame(percentiles, specificity_list, tss_list, thresh_val_list)
colnames(tss_df) <- c("Sensitivity", "Specificty", "TSS", 
                      "Suitability")
write.csv(tss_df, file = './Model_results/l_chalumnae_2D_sdm_tss_results.csv')

# write raster to file
writeRaster(final_sdm, file = './Model_results/coelacanth_2D_thresholded.tif', 
            format = "GTiff", overwrite = T)
