# bivariate chloropleth plots of 3D suitability

library(sf)
library(terra)
library(ggplot2)
library(ggspatial)
library(paletteer)
library(dplyr)
library(rnaturalearth)

# creating land for plotting
target_crs <- "+proj=longlat +datum=WGS84"
land <- st_as_sf(rnaturalearthhires::countries10)[1] %>% 
  st_set_crs(target_crs)

# reading in data
chalumnae_3D_suit <- rast('./Model_results/coelacanth_3D_suitability.tif')
indo_3D_suit <- rast('./Model_results/Indopacific_projection_3D_notmasked.tif')
conductivity <- rast('./WOA_18/conductivity_summer_18.tif')

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
      mat_list[[i]] <- mat1
    }
    mat_list_all <- do.call(rbind, mat_list)
    mat_list_all <- data.frame(mat_list_all)
    colnames(mat_list_all) <- c("x", "y", "z", "s")
  }
  return(mat_list_all)
}

chalumnae_suit_mat <- rast_to_mat(chalumnae_3D_suit, depth_slices = depths, thresholded = F)
indo_suit_mat <- rast_to_mat(indo_3D_suit, depth_slices = depths, thresholded = F)

# Finding breakpoints
# chalumnae
maxdepth <- min(chalumnae_suit_mat$z)
break1 <- -100
break2 <- -400
depth_ax <- vector(length = length(chalumnae_suit_mat$z))
for (i in 1:length(chalumnae_suit_mat$z)) {
  if(chalumnae_suit_mat$z[i] >= break1) {
    depth_ax[i] <- "A"
  } else if(chalumnae_suit_mat$z[i] < break1 & chalumnae_suit_mat$z[i] > break2) {
    depth_ax[i] <- "B"
  } else {
    depth_ax[i] <- "C"
  }
}

maxsuit <- max(chalumnae_suit_mat$s)
break1 <- maxsuit/3
break2 <- break1*2
suit_ax <- vector(length = length(chalumnae_suit_mat$s))
for (i in 1:length(chalumnae_suit_mat$s)) {
  if(chalumnae_suit_mat$s[i] <= break1) {
    suit_ax[i] <- "1"
  } else if(chalumnae_suit_mat$s[i] > break1 & chalumnae_suit_mat$s[i] < break2) {
    suit_ax[i] <- "2"
  } else {
    suit_ax[i] <- "3"
  }
}

# indo projection
maxdepth <- min(indo_suit_mat$z)
break1 <- -100
break2 <- -400
depth_ax_indo <- vector(length = length(indo_suit_mat$z))
for (i in 1:length(indo_suit_mat$z)) {
  if(indo_suit_mat$z[i] >= break1) {
    depth_ax_indo[i] <- "A"
  } else if(indo_suit_mat$z[i] < break1 & indo_suit_mat$z[i] > break2) {
    depth_ax_indo[i] <- "B"
  } else {
    depth_ax_indo[i] <- "C"
  }
}

maxsuit <- max(indo_suit_mat$s)
break1 <- maxsuit/3
break2 <- break1*2
suit_ax_indo <- vector(length = length(indo_suit_mat$s))
for (i in 1:length(indo_suit_mat$s)) {
  if(indo_suit_mat$s[i] <= break1) {
    suit_ax_indo[i] <- "1"
  } else if(indo_suit_mat$s[i] > break1 & indo_suit_mat$s[i] < break2) {
    suit_ax_indo[i] <- "2"
  } else {
    suit_ax_indo[i] <- "3"
  }
}

chlor <- vector(length = nrow(chalumnae_suit_mat))
for (i in 1:nrow(chalumnae_suit_mat)) {
  chlor[i] <- paste0(depth_ax[i], suit_ax[i])
}

chlor_indo <- vector(length = nrow(indo_suit_mat))
for (i in 1:nrow(indo_suit_mat)) {
  chlor_indo[i] <- paste0(depth_ax_indo[i], suit_ax_indo[i])
}

chalumnae_suit_mat <- cbind(chalumnae_suit_mat, chlor)
indo_suit_mat <- cbind(indo_suit_mat, chlor_indo)

# selecting for the unique XY positions where suitability is highest
coords_chalumnae <- crds(chalumnae_3D_suit[[16]])
save_list <- vector("list", length = nrow(coords_chalumnae))
for (i in 1:nrow(coords_chalumnae)) {
  temp_df <- chalumnae_suit_mat %>% filter(x == coords_chalumnae[i,1]) %>% 
    filter(y == coords_chalumnae[i,2])
  save_this <- temp_df[which(temp_df$s == max(temp_df$s)),]
  save_list[[i]] <- save_this
}
chloro_suit_chalumnae <- do.call(rbind, save_list)

coords_indo <- crds(indo_3D_suit[[16]])
save_list_indo <- vector("list", length = nrow(coords_indo))
for (i in 1:nrow(coords_indo)) {
  temp_df <- indo_suit_mat %>% filter(x == coords_indo[i,1]) %>% 
    filter(y == coords_indo[i,2])
  save_this <- temp_df[which(temp_df$s == max(temp_df$s)),]
  save_list_indo[[i]] <- save_this
}
chloro_suit_indo <- do.call(rbind, save_list_indo)

# Plotting
names(chloro_suit_chalumnae) <- c("lon", "lat", "z", "s", "chlor")

names(chloro_suit_indo) <- c("lon", "lat", "z", "s", "chlor")

tiff('./Plots/bivariate_chloropleth_chalumnae_3D.tiff', width = 3400, height = 3400,
     res = 300)
ggplot() +
  geom_tile(data = chloro_suit_chalumnae, aes(x = lon, y = lat, fill = chlor)) +
  geom_sf(data = land, fill = "darkgrey", linewidth = 0.6) +
  coord_sf(xlim = c(9.25, 63.75), ylim = c(-48.75, 14.25), expand = F) +
  scale_fill_manual(values = c("A1" = "#e8e8e8", "A2" = "#dfb0d6", "A3" = "#be64ac",
                               "B1" = "#ace4e4", "B2" = "#a5add3", "B3" = "#8c62aa",
                               "C1" = "#5ac8c8", "C2" = "#5698b9", "C3" = "#3b4994")) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 21), axis.title.y = element_text(size = 21),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9))
dev.off()

tiff('./Plots/bivariate_chloropleth_indo_3D.tiff', width = 3400, height = 3400, res = 300)
ggplot() +
  geom_tile(data = chloro_suit_indo, aes(x = lon, y = lat, fill = chlor)) +
  geom_sf(data = land, fill = "darkgrey", linewidth = 0.6) +
    coord_sf(xlim = c(90, 150), ylim = c(-30, 30), expand = F) +
  scale_fill_manual(values = c("A1" = "#e8e8e8", "A2" = "#dfb0d6", "A3" = "#be64ac",
                               "B1" = "#ace4e4", "B2" = "#a5add3", "B3" = "#8c62aa",
                               "C1" = "#5ac8c8", "C2" = "#5698b9", "C3" = "#3b4994")) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 21), axis.title.y = element_text(size = 21),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9))
dev.off()
  
