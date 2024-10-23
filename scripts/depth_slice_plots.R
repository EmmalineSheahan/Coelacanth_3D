# depth Slice plots (Depth vs Longitude)

library(sf)
library(terra)
library(ggplot2)
library(ggspatial)
library(paletteer)
library(dplyr)
library(rnaturalearth)
library(viridis)
library(ds4psy)
library(gridExtra)

# creating land for plotting
target_crs <- "+proj=longlat +datum=WGS84"
land <- st_as_sf(rnaturalearthhires::countries10)[1] %>% 
  st_set_crs(target_crs)

# reading in data
chalumnae_3D_thresh <- rast('./Model_results/coelacanth_3D_thresholded.tif')
indo_3D_thresh <- rast('./Model_results/Indopacific_thresholded_3D.tif')
conductivity <- rast('./WOA_18/env_3D_conductivity_winter_interpolated.tif')
temperature <- rast('./WOA_18/env_3D_temperature_summer_interpolated.tif')
north_current <- rast('./WOA_18/env_3D_current_north_december_interpolated.tif')
slope <- rast('./WOA_18/env_3D_slope_interpolated.tif')
bath_new <- rast('./Slope/ETOPO1_Bed_c_geotiff.tif')

conductivity_indo <- rast('./Envs_Indo/env_indo_3D_conductivity_winter_interpolated.tif')
temperature_indo <- rast('./Envs_Indo/env_indo_3D_temperature_summer_interpolated.tif')
north_current_indo <- rast('./Envs_Indo/env_indo_3D_current_north_december_interpolated.tif')
slope_indo <- rast('./Envs_Indo/env_indo_3D_slope_interpolated.tif')

depths <- as.numeric(gsub("X", '', names(conductivity)))

# rast_to_mat function to transform 3D spatraster stacks into xyz matrices
# to interpret
# wanted_rast = spatraster stack of thresholded or suitability rasters
# depth_slices = vector of depth values, assumes positive
# thresholded = T or F, true if stack is thresholded
rast_to_mat <- function(wanted_rast, depth_slices, thresholded = T) {
  if(thresholded) {
    for(i in 1:dim(wanted_rast)[3]) {
      values(wanted_rast[[i]])[which(is.na(values(wanted_rast[[i]])))] <- 0
    }
    mat_list <- vector("list", length = dim(wanted_rast)[3])
    for (i in 1:dim(wanted_rast)[3]) {
      wanted_coords <-  xyFromCell(wanted_rast[[i]], 
                                   which(wanted_rast[[i]][] == 1))
      z <- rep(-as.numeric(depth_slices[i]), times = nrow(wanted_coords))
      wanted_mat <- cbind(wanted_coords, z)
      mat_list[[i]] <- wanted_mat
    }
    mat_list_all <- do.call(rbind, mat_list)
    mat_list_all <- data.frame(mat_list_all)
  } else {
    mat_list <- vector("list", length = dim(wanted_rast)[3])
    for (i in 1:dim(wanted_rast)[3]) {
      p <- xyFromCell(wanted_rast[[i]], cell = which(!(is.na(wanted_rast[[i]][]))))
      z <- rep(-as.numeric(depth_slices[i]), times = nrow(p))
      s <- extract(wanted_rast[[i]], which(!(is.na(wanted_rast[[i]][]))))
      mat1 <- cbind(p, z, s)
      colnames(mat1) <- c("x", "y", "z", "s")
      mat_list[[i]] <- mat1
    }
    mat_list_all <- do.call(rbind, mat_list)
    mat_list_all <- data.frame(mat_list_all)
    colnames(mat_list_all) <- c("x", "y", "z", "s")
  }
  return(mat_list_all)
}

# chalumnae matrices
chalumnae_thresh_mat <- rast_to_mat(wanted_rast = chalumnae_3D_thresh, depth_slices = depths,
                                    thresholded = T)
conductivity_mat <- rast_to_mat(wanted_rast = conductivity, depth_slices = depths,
                                thresholded = F)
temperature_mat <- rast_to_mat(wanted_rast = temperature, depth_slices = depths,
                               thresholded = F)
current_mat <- rast_to_mat(wanted_rast = north_current, depth_slices =  depths,
                           thresholded = F)
slope_mat <- rast_to_mat(wanted_rast = slope, depth_slices = depths,
                         thresholded = F)

# indo matrices
indo_thresh_mat <- rast_to_mat(wanted_rast = indo_3D_thresh, depth_slices = depths,
                               thresholded = T)
indo_temperature_mat <- rast_to_mat(wanted_rast = temperature_indo, depth_slices = depths,
                               thresholded = F)

# selecting depth slices and filling in for plotting
# comoros
comoros_thresh <- chalumnae_thresh_mat %>% filter(y == -12.125)
s <- rep(1, times = nrow(comoros_thresh))
comoros_thresh <- cbind(comoros_thresh, s)
needed_slice <- depths[21:length(depths)]
extra_slices_list <- vector("list", length = length(needed_slice) - 1)
for(i in 1:(length(needed_slice) - 1)) {
  new_vec <- (needed_slice[i]:needed_slice[i+1])
  new_div <- new_vec/5
  new_fac <- new_vec[which(is_wholenumber(new_div))]
  new_fac_df_list <- vector("list", length = length(new_fac))
  for (j in 2:(length(new_fac) - 1)) {
    df <- comoros_thresh %>% filter(z == -(needed_slice[i+1]))
    df$z <- rep(-new_fac[j], times = nrow(df))
    new_fac_df_list[[j]] <- df
  }
  new_fac_df <- do.call(rbind, new_fac_df_list)
  extra_slices_list[[i]] <- new_fac_df
}
extra_slices <- do.call(rbind, extra_slices_list)
comoros_thresh <- rbind(comoros_thresh, extra_slices)

comoros_temp <- temperature_mat %>% filter(y == -12.125)

# South Africa
sa_thresh <- chalumnae_thresh_mat %>% filter(y == -27.125)
s <- rep(1, times = nrow(sa_thresh))
sa_thresh <- cbind(sa_thresh, s)
needed_slice <- depths[21:length(depths)]
extra_slices_list <- vector("list", length = length(needed_slice) - 1)
for(i in 1:(length(needed_slice) - 1)) {
  new_vec <- (needed_slice[i]:needed_slice[i+1])
  new_div <- new_vec/5
  new_fac <- new_vec[which(is_wholenumber(new_div))]
  new_fac_df_list <- vector("list", length = length(new_fac))
  for (j in 2:(length(new_fac) - 1)) {
    df <- sa_thresh %>% filter(z == -(needed_slice[i+1]))
    df$z <- rep(-new_fac[j], times = nrow(df))
    new_fac_df_list[[j]] <- df
  }
  new_fac_df <- do.call(rbind, new_fac_df_list)
  extra_slices_list[[i]] <- new_fac_df
}
extra_slices <- do.call(rbind, extra_slices_list)
sa_thresh <- rbind(sa_thresh, extra_slices)

sa_temp_reduced <- temperature_mat %>% filter(y == -27.125) %>% filter(x == 33.125)

# Indonesia
indo_thresh <- indo_thresh_mat %>% filter(y == -7.625)
s <- rep(1, times = nrow(indo_thresh))
indo_thresh <- cbind(indo_thresh, s)
needed_slice <- depths[21:length(depths)]
extra_slices_list <- vector("list", length = length(needed_slice) - 1)
for(i in 1:(length(needed_slice) - 1)) {
  new_vec <- (needed_slice[i]:needed_slice[i+1])
  new_div <- new_vec/5
  new_fac <- new_vec[which(is_wholenumber(new_div))]
  new_fac_df_list <- vector("list", length = length(new_fac))
  for (j in 2:(length(new_fac) - 1)) {
    df <- indo_thresh %>% filter(z == -(needed_slice[i+1]))
    df$z <- rep(-new_fac[j], times = nrow(df))
    new_fac_df_list[[j]] <- df
  }
  new_fac_df <- do.call(rbind, new_fac_df_list)
  extra_slices_list[[i]] <- new_fac_df
}
extra_slices <- do.call(rbind, extra_slices_list)
indo_thresh <- rbind(indo_thresh, extra_slices)

indo_temp_reduced <- indo_temperature_mat %>% filter(y == -7.625) %>% filter(x == 125.375)

# bath_depth function to pull bathymetry at a given longitude or latitude line for 
# depth slice plotting
# bath = bathymetry layer
# wanted_projection = spatraster with desired resolution so values match
# wanted_axis = either 'lon' or 'lat'
# wanted_val = the given longitude or latitude
# maxdepth = maximum depth of depth slice plot (negative)
# depth_slices = vector of depths to match (negative)
bath_depth <- function(wanted_projection, bath, wanted_axis, wanted_val, maxdepth, 
                       depth_slices) {
  bath <- terra::project(bath, wanted_projection)
  bathxy <- data.frame(crds(bath))
  bathz <- values(bath)
  bathxyz <- cbind(bathxy, bathz)
  colnames(bathxyz) <- c("x", "y", "z")
  want_bath <- bathxyz %>% filter(z > maxdepth)
  s <- rep(NA, times = nrow(want_bath))
  want_bath <- cbind(want_bath, s)
  if(wanted_axis == 'lat') {
    want_bath_final <- want_bath %>% filter(y == wanted_val)
  } else if(wanted_axis == 'lon') {
    want_bath_final <- want_bath %>% filter(x == wanted_val)
  }
  rowadd_list <- vector("list", length = nrow(want_bath_final))
  for (i in 1:nrow(want_bath_final)) {
    numreprow <- want_bath_final[i,]$z
    if(numreprow > 0) {
      replacedepth <- 0
    } else {
      replacedepth <- depth_slices[which(depth_slices > 
                                         numreprow)][length(which(depth_slices > numreprow))]
    }
    want_bath_final[i, 3] <- replacedepth
    starthere <- which(depth_slices < replacedepth)[1]
    seqalongthese1 <- depth_slices[starthere]:depth_slices[length(depth_slices)]
    fac <- seqalongthese1/5
    seqalongthese <- seqalongthese1[which(is_wholenumber(fac))]
    addline_list <- vector("list", length = length(seqalongthese))
    for (j in seq_along(seqalongthese)) {
      addline <- want_bath_final[i,]
      addline$z <- seqalongthese[j]
      addline_list[[j]] <- addline
    }
    addline_full <- do.call(rbind, addline_list)
    rowadd_list[[i]] <- addline_full
  }
  rowadd <- do.call(rbind, rowadd_list)
  want_bath_done <- rbind(want_bath_final, rowadd)
  return(want_bath_done)
}

comoros_bath <- bath_depth(wanted_projection = chalumnae_3D_thresh[[1]], bath = bath_new, 
                           wanted_axis = 'lat', wanted_val = -12.125, maxdepth = -800,
                           depth_slices = -depths)

sa_bath <- bath_depth(wanted_projection = chalumnae_3D_thresh[[1]], bath = bath_new, 
                      wanted_axis = 'lat', wanted_val = -27.125, maxdepth = -800,
                      depth_slices = -depths)
indo_bath <- bath_depth(wanted_projection = indo_3D_thresh[[1]], bath = bath_new,
                        wanted_axis = 'lat', wanted_val = -7.625, maxdepth = -800,
                        depth_slices = -depths)

comoros_thresh_new <- rbind(comoros_thresh, comoros_bath)
sa_thresh_new <- rbind(sa_thresh, sa_bath)
indo_thresh_new <- rbind(indo_thresh, indo_bath)

# plotting
tiff('./Plots/comoros_SA_depthslice_full_figure.tiff', width = 3400, height = 4400, 
     res = 300)
t1 <- ggplot() +
  geom_sf(data = land, fill = "darkgrey", linewidth = 0.6) +
  coord_sf(xlim = c(10, 55), ylim = c(5, -45)) +
  geom_hline(yintercept = -12.125, color = "#FC6A03") +
  geom_hline(yintercept = -27.125, color = "#BF40BF") +
  theme(panel.background = element_rect(fill = "aliceblue")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9))

t2 <- ggplot(data = comoros_thresh_new, aes(x, z, fill = s)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("white", "#03C04A"), 
                       values = c(0, 1),
                       na.value = "darkgrey") +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_text(size = 15), axis.text = element_text(size = 15),
        legend.position = "none",
        axis.line = element_line(linewidth = 0.9, linetype = "solid", colour = "black")) +
  xlab("Longitude") +
  ylab("Depth (m)") +
  xlim(40, 50)

t3 <- ggplot(data = sa_thresh_new, aes(x, z, fill = s)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("white", "#03C04A"), 
                       values = c(0, 1),
                       na.value = "darkgrey") +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_text(size = 15), axis.text = element_text(size = 15),
        legend.position = "none",
        axis.line = element_line(linewidth = 0.9, linetype = "solid", colour = "black")) +
  xlab("Longitude") +
  ylab("Depth (m)") +
  xlim(32, 37) +
  ylim(-800, 0)

comoros_temp_reduced <- comoros_temp %>% filter(x == 43.875)
t4 <- ggplot(data = comoros_temp_reduced, aes(s, z)) +
  geom_line(colour = "blue") +
  ylab("Depth (m)") +
  xlab("Temperature (°C)") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 15), 
        panel.background = element_rect(fill = "white"), 
        axis.line = element_line(linewidth = 0.9, linetype = "solid", colour = "black")) 

t5 <- ggplot(data = sa_temp_reduced, aes(s, z)) +
  geom_line(colour = "blue") +
  ylab("Depth (m)") +
  xlab("Temperature (°C)") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 15), 
        panel.background = element_rect(fill = "white"), 
        axis.line = element_line(linewidth = 0.9, linetype = "solid", colour = "black"))

first_plot <- grid.arrange(t2, t3, t4, t5, ncol = 2)
grid.arrange(t1, first_plot, nrow = 2)
dev.off()

tiff('./Plots/indo_depthslice_full_figure.tiff', width = 3400, height = 4400, res = 300)
t1 <- ggplot() +
  geom_sf(data = land, fill = "darkgrey", linewidth = 0.6) +
  coord_sf(xlim = c(90, 150), ylim = c(30, -30)) +
  geom_hline(yintercept = -7.625, color = "#FC6A03", size = 0.7) +
  theme(panel.background = element_rect(fill = "aliceblue")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9))

t2 <- ggplot(data = indo_thresh_new, aes(x, z, fill = s)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("white", "#03C04A"), 
                       values = c(0, 1),
                       na.value = "darkgrey") +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_text(size = 15), axis.text = element_text(size = 15),
        legend.position = "none",
        axis.line = element_line(linewidth = 0.9, linetype = "solid", colour = "black"),
        aspect.ratio = 1) +
  xlab("Longitude") +
  ylab("Depth (m)") +
  xlim(110, 140) +
  ylim(-800, 0)

t3 <- ggplot(data = indo_temp_reduced, aes(s, z)) +
  geom_line(colour = "blue") +
  ylab("Depth (m)") +
  xlab("Temperature (°C)") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 15), 
        panel.background = element_rect(fill = "white"), 
        axis.line = element_line(linewidth = 0.9, linetype = "solid", colour = "black"),
        aspect.ratio = 1)
grid.arrange(t1, t2, t3, ncol = 1)
dev.off()


