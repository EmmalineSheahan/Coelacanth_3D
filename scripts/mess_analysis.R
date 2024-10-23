# MESS Analysis

package_list <- c("ggplot2", "modEvA", "tidyterra", "viridis", "doParallel")
install.packages(package_list)

library(voluModel)
library(ggplot2)
library(dplyr)
library(terra)
library(sf)
library(modEvA)
library(tidyterra)
library(viridis)
library(doParallel)
library(rnaturalearthhires)

# creating land for plotting
target_crs <- "+proj=longlat +datum=WGS84"
land <- st_as_sf(rnaturalearthhires::countries10)[1] %>% 
  st_set_crs(target_crs)

# reading in environmental datasets
# for 2D
all_2D_env <- rast('./Envs_2D/all_2D_env_final.tif')
indo_2D_env <- rast('./data/temporary_2D_indo_env_real.tif')

pull_for <- which(names(all_2D_env) %in% names(indo_2D_env))
all_2D_env_final <- all_2D_env[[pull_for]]

# for 3D
conductivity_winter <- rast('./WOA_18/env_3D_conductivity_winter_interpolated.tif')
temperature_summer <- rast('./WOA_18/env_3D_temperature_summer_interpolated.tif')
current_north_december <- 
  rast('./WOA_18/env_3D_current_north_december_interpolated.tif')
slope <- rast('./WOA_18/env_3D_slope_interpolated.tif')
conductivity_winter_indo <- 
  rast('./Envs_Indo/env_indo_3D_conductivity_winter_interpolated.tif')
temperature_summer_indo <- 
  rast('./Envs_Indo/env_indo_3D_temperature_summer_interpolated.tif')
current_north_december_indo <- 
  rast('./Envs_Indo/env_indo_3D_current_north_december_interpolated.tif')
slope_indo <- rast('./Envs_Indo/env_indo_3D_slope_interpolated.tif')

# Creating mess dataframes
# for 2D
want_dat_list <- vector("list", length = length(names(all_2D_env_final)))
for (i in 1:length(names(all_2D_env_final))) {
  want_dat_list[[i]] <- values(all_2D_env_final[[i]])[,1]
}
V_2D <- do.call(cbind, want_dat_list)
colnames(V_2D) <- names(all_2D_env_final)

want_dat_list_indo <- vector("list", length = length(names(indo_2D_env)))
for (i in 1:length(names(indo_2D_env))) {
  want_dat_list_indo[[i]] <- values(indo_2D_env[[i]])[,1]
}
P_2D <- do.call(cbind, want_dat_list_indo)
colnames(P_2D) <- names(indo_2D_env)

# for 3D
envs_3D <- list(conductivity_winter, current_north_december, slope, temperature_summer)
envs_3D_indo <- list(conductivity_winter_indo, current_north_december_indo, 
                     slope_indo, temperature_summer_indo)

all_depths <- as.numeric(gsub("X", '', names(conductivity_winter)))

depth_list <- vector("list", length = length(all_depths))
want_dat_list <- vector("list", length = length(envs_3D))
for (j in seq_along(all_depths)) {
  for (i in 1:length(envs_3D)) {
    want_dat_list[[i]] <- values(envs_3D[[i]][[j]])[,1]
  }
  depth_list[[j]] <- do.call(cbind, want_dat_list)
  colnames(depth_list[[j]]) <- c("conductivity_winter", "current_north_december",
                                 "slope", "temperature_summer")
  depth_list[[j]] <- data.frame(depth_list[[j]])
}
cal_3D_list <- depth_list

depth_list_indo <- vector("list", length = length(all_depths))
want_dat_list_indo <- vector("list", length = length(envs_3D_indo))
for (j in seq_along(all_depths)) {
  for (i in 1:length(envs_3D_indo)) {
    want_dat_list_indo[[i]] <- values(envs_3D_indo[[i]][[j]])[,1]
  }
  depth_list_indo[[j]] <- do.call(cbind, want_dat_list_indo)
  colnames(depth_list_indo[[j]]) <- c("conductivity_winter", "current_north_december",
                                 "slope", "temperature_summer")
  depth_list_indo[[j]] <- data.frame(depth_list_indo[[j]])
}
proj_3D_list <- depth_list_indo

# attempting mess functions
# for 2D
messres_2D <- MESS(V = V_2D, P = P_2D)

xy_2D <- crds(indo_2D_env, df = T)
mess_xyz_2D <- cbind(xy_2D, messres_2D$TOTAL)
colnames(mess_xyz_2D) <- c("x", "y", "z")
mess_2D_spat <- rast(mess_xyz_2D, type = "xyz", crs = target_crs)

values(mess_2D_spat)[which(values(mess_2D_spat) < -125)] <- NA

pdf('./Plots/MESS_2D_indo.pdf', pointsize = 3)
ggplot() +
  geom_spatraster(data = mess_2D_spat) +
  geom_sf(data = land) +
  coord_sf(xlim = c(78.75, 150), ylim = c(-30, 39), expand = F) +
  scale_fill_gradientn(colors = c("white", viridis(4)), 
                       breaks = c(0, -25, -50, -75, -100),
                       trans = "reverse",
                       na.value = "transparent",
                       name = "MESS Total") +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"))
dev.off()

# for 3D
cl <- makeCluster(2)
registerDoParallel(cl, cores = 2)

messres_3D_list <- vector("list", length = length(cal_3D_list))
messspat_3D_list <- vector("list", length = length(cal_3D_list))
for (i in 1:length(cal_3D_list)) {
  temp <- MESS(V = cal_3D_list[[i]], P = proj_3D_list[[i]])
  messres_3D_list[[i]] <- temp
  xy_3D <- crds(slope_indo, df = T)
  mess_xyz_3D <- cbind(xy_3D, temp$TOTAL)
  colnames(mess_xyz_3D) <- c("x", "y", "z")
  messspat_3D_list[[i]] <- rast(mess_xyz_3D, type = "xyz", crs = target_crs)
  values(messspat_3D_list[[i]])[which(values(messspat_3D_list[[i]]) < -125)] <- NA
  print(paste0("completed iteration ", i))
}

stopCluster(cl)

messspat_3D <- rast(messspat_3D_list)
writeRaster(messspat_3D, './Model_results/mess_3D.tif', overwrite = T)

pdf('./Plots/MESS_3D_indo.pdf', pointsize = 3)
for (i in 1:dim(messspat_3D)[3]) {
  p <- ggplot() +
    geom_spatraster(data = messspat_3D[[i]]) +
    geom_sf(data = land) +
    coord_sf(xlim = c(78.75, 150), ylim = c(-30, 39), expand = F) +
    scale_fill_gradientn(colors = c("white", viridis(4)), 
                       breaks = c(0, -25, -50, -75, -100),
                       trans = "reverse",
                       na.value = "transparent",
                       name = "MESS Total") +
    theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white")) +
    ggtitle(paste0("MESS at ", wanted_slices[i], " m"))
  print(p)
}
dev.off()