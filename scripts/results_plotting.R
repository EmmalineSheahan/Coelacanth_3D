# all 2D and 3D model plotting

library(ggplot2)
library(sf)
library(rnaturalearthhires)
library(terra)
library(tidyterra)
library(dplyr)
library(raster)
library(voluModel)
library(viridis)
library(plotly)
library(ds4psy)
library(gridExtra)

# creating land for plotting
target_crs <- "+proj=longlat +datum=WGS84"
land <- st_as_sf(rnaturalearthhires::countries10)[1] %>% 
  st_set_crs(target_crs)

# reading in 2D suitability and thresholded SpatRasters
suit_2D <- rast('./Model_results/coelacanth_2D_suitability.tif')
thresh_2D <- rast('./Model_results/coelacanth_2D_thresholded.tif')

# ugly plots
pdf('./plots/unofficial_2D_suit.pdf')
plot(suit_2D, main = "Coelacanth 2D Suitability")
plot(land, add = T, col = NA)
dev.off()

pdf('./plots/unofficial_2D_thresh.pdf')
plot(thresh_2D, main = "Coelacanth 2D Thresholded")
plot(land, add = T, col = NA)
dev.off()

# plotting 2D suitability
tiff('./Plots/l_chalumnae_2D_suit_final.tiff',
     width = 3400,
     height = 3400,
     res = 300)
t1 <- ggplot() +
  geom_spatraster(data = suit_2D) +
  geom_sf(data = land, fill = "darkgrey", linewidth = 0.6) +
  coord_sf(xlim = c(15, 57), ylim = c(-45, 7), expand = F) +
  scale_fill_gradientn(colors = c("white", viridis(4)), 
              values = c(0, 0.25, 0.5, 0.75, 1),
              na.value = "transparent",
              name = "Suitability Score") +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 17), axis.text.y = element_text(size = 17),
        legend.text = element_text(size = 15), legend.title = element_text(size = 17),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
        legend.key.size = unit(0.6, "cm"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black", linewidth = 0.8),
        legend.position = c(0.19, 0.87), title = element_text(size = 18)) +
  ggtitle("Raw Suitability")
print(t1)
dev.off()

# plotting 2D threshold
tiff('./Plots/l_chalumnae_2D_thresh_final.tiff',
     width = 3400,
     height = 3400,
     res = 300)
t2 <- ggplot() +
  geom_spatraster(data = thresh_2D) +
  geom_sf(data = land, fill = "darkgrey", linewidth = 0.6) +
  coord_sf(xlim = c(15, 57), ylim = c(-45, 7), expand = F) +
  scale_fill_gradientn(colors = c("white", "darkgreen"), values = c(0,1),
                       na.value = "transparent") +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.text.x = element_text(size = 17), axis.text.y = element_text(size = 17),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
        title = element_text(size = 18)) +
  ggtitle("Presence/Absence")
print(t2)
dev.off()

# chalumnae suitabilty and threshold together
tiff("./Plots/figure_3_real.tiff", width = 3400,
     height = 3400,
     res = 300)
grid.arrange(t1, t2, nrow = 1)
dev.off()

# reading in 3D suitability and thresholded rasters
suit_3D <- rast('./Model_results/coelacanth_3D_suitability.tif')
thresh_3D <- rast('./Model_results/coelacanth_3D_thresholded.tif')

# plotlayers plot of 3D threshold
t <- rast('./WOA_18/density_winter_18.tif')
depth_list <- gsub('X', '', names(t))
names(thresh_3D) <- paste0(depth_list, " m")

pdf('./Plots/l_chalumnae_3D_plotLayers.pdf', pointsize = 3)
plotLayers(rast = thresh_3D, land = land)
dev.off()

# plotting 3D suitability layer by layer
pdf('./Plots/l_chalumnae_3D_suit_all_layers.pdf', pointsize = 3)
for(i in 1:dim(suit_3D)[3]) {
p <- ggplot() +
  geom_spatraster(data = suit_3D[[i]]) +
  geom_sf(data = land, fill = "darkgrey") +
  coord_sf(xlim = c(9.25, 63.75), ylim = c(-48.75, 14.25), expand = F) +
  scale_fill_gradientn(colors = c("white", viridis(4)), 
                       values = c(0, 0.25, 0.5, 0.75, 1),
                       na.value = "transparent",
                       name = "Suitability Score") +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 13), legend.title = element_text(size = 13),
        plot.title = element_text(size = 17),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_rect(fill = "white", linetype = "solid", 
                                         colour = "black"),
        legend.position = c(0.14, 0.904)) +
  ggtitle(paste0("Suitability at ", depth_list[i], " m"))
print(p)
}
dev.off()

pdf('./Plots/l_chalumnae_3D_thresh_all_layers.pdf', pointsize = 3)
for(i in 1:dim(thresh_3D)[3]) {
  p <- ggplot() +
    geom_spatraster(data = thresh_3D[[i]]) +
    geom_sf(data = land, fill = "darkgrey") +
    coord_sf(xlim = c(9.25, 63.75), ylim = c(-48.75, 14.25), expand = F) +
    scale_fill_gradientn(colors = c("white", "darkgreen"), 
                         values = c(0, 1),
                         na.value = "transparent") +
    theme(panel.grid.major = element_blank(), 
          panel.background = element_rect(fill = "white"), legend.position = "none",
          axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
          legend.text = element_text(size = 15), legend.title = element_text(size = 15),
          plot.title = element_text(size = 17),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
    ggtitle(paste0("Thresholded Presence at ", depth_list[i], " m"))
  print(p)
}
dev.off()

# plotly plot of 3D threshold and suitability

# rast_to_mat function to transform 3D spatraster stacks into xyz matrices for plotly
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

# coast_3D function to make coastline for plotly 
# bath_file = filepath for bathymetry layer
# wanted_extent = extent to crop to
# depth_lim = depth limit
# returns a list with a matrix of just bathymetry values for z plotting, a matrix of 
# x coordinates and a matrix of y coordinates
coast_3D <- function(bath_file, wanted_extent, depth_lim){
  bath <- terra::rast(bath_file)
  bathcrop <- crop(bath, wanted_extent)
  values(bathcrop)[which(values(bathcrop) > 0)] <- 0
  values(bathcrop)[which(values(bathcrop) < -depth_lim)] <- -depth_lim
  bathmat <- as.matrix(bathcrop)
  bathcoords_new <- crds(bathcrop)
  bathx <- matrix(bathcoords_new[,1], nrow = dim(bathmat)[1], 
                      ncol = dim(bathmat)[2],
                      byrow = T)
  bathy <- matrix(bathcoords_new[,2], nrow = dim(bathmat)[1], 
                      ncol = dim(bathmat)[2],
                      byrow = T)
  bath_list <- list(bathmat = bathmat, bathx = bathx, bathy = bathy)
  return(bath_list)
}

bath_obj <- coast_3D(bath_file = "./Slope/ETOPO1_Bed_c_geotiff.tif", 
                   wanted_extent = thresh_2D, depth_lim = 900)

# 3D threshold plot
thresh_mat <- rast_to_mat(wanted_rast = thresh_3D, depth_slices = depth_list,
                          thresholded = T)
thresh_3D_plot <- plot_ly(data = thresh_mat, x = ~x, y = ~y, z = ~z, 
                          type = 'scatter3d',
                          marker = list(color = "blue", opacity = 0.25),
                          name = "Thresholded Coelacanth Model")
thresh_3D_plot %>% add_trace(x = ~bath_obj$bathx, y = ~bath_obj$bathy, 
                             z = ~bath_obj$bathmat, type = "surface", 
                             name = "South African Continent",
                             colorscale = list(c(0, 1), c("gray", "gray")),
                             showscale = F)

# 3D suitability plot
suit_mat <- rast_to_mat(wanted_rast = suit_3D, depth_slices = depth_list,
                        thresholded =  F)
plot_ly(data = suit_mat, x = ~x, y = ~y, z = ~z, 
        type = 'scatter3d', color = ~s, alpha = 0.25)

# 3 panel plots

# contour_mat function to create a 2D contour dataframe of whichever side of the 3D
# matrix you're interested in plotting
# this is to ensure you're only plotting what would be visible from a given perspective
# and adequately accounting for hidden points
# wanted_mat = xyz matrix of coordinates
# wanted_z = character column name of which vector is serving as the z axis
# wanted_x = character column name of which vector is serving as the x axis
# wanted_y = character column name of which vector is serving as the y axis
# orientation = "top", "east", "south", "west", "north", or "bottom". 
# Viewing perspective
# thresholded = T or F, if the matrix is thresholded or not. if it isn't and is 
# a suitability matrix, you must supply wanted_val
# wanted_val = column name of suitability value
contour_mat <- function(wanted_mat, wanted_z, wanted_x, wanted_y, orientation, 
                        thresholded = T, wanted_val) {
  z_vec <- unique(wanted_mat %>% dplyr::select(wanted_z))[,1]
  z_list <- vector("list", length = length(z_vec))
  if(orientation == "top") {
    z_list[[1]] <- wanted_mat %>% filter(eval(parse(text = wanted_z)) 
                                         == max(z_vec, na.rm = T))
  } else if(orientation == "east") {
    z_list[[1]] <- wanted_mat %>% filter(eval(parse(text = wanted_z)) 
                                         == max(z_vec, na.rm = T))
  } else if(orientation == "south") {
    z_list[[1]] <- wanted_mat %>% filter(eval(parse(text = wanted_z)) 
                                         == min(z_vec, na.rm = T))
  } else if(orientation == "west") {
    z_list[[1]] <- wanted_mat %>% filter(eval(parse(text = wanted_z)) 
                                         == min(z_vec, na.rm = T))
  } else if(orientation == "north") {
    z_list[[1]] <- wanted_mat %>% filter(eval(parse(text = wanted_z)) 
                                         == max(z_vec, na.rm = T))
  } else if(orientation == "bottom") {
    z_list[[1]] <- wanted_mat %>% filter(eval(parse(text = wanted_z)) 
                                         == min(z_vec, na.rm = T))
  } else {
    print("not a usable orientation")
  }
  if(orientation == "top" | orientation == "east" | orientation == "north") {
    z_vec <- z_vec[order(z_vec, decreasing = T)]
    for (i in 2:length(z_vec)) {
      searchdf <- wanted_mat %>% filter(eval(parse(text = wanted_z)) 
                                        > as.numeric(z_vec[i]))
      thisdf <- wanted_mat %>% filter(eval(parse(text = wanted_z)) 
                                      == as.numeric(z_vec[i]))
      t <- searchdf %>% dplyr::select(wanted_x, wanted_y)
      p <- thisdf %>% dplyr::select(wanted_x, wanted_y)
      tchar <- vector(length = nrow(t))
      for (j in 1:nrow(t)) {
        tchar[j] <- paste0(t[j,1], t[j,2])
      }
      pchar <- vector(length = nrow(p))
      for (j in 1:nrow(p)) {
        pchar[j] <- paste0(p[j,1], p[j,2])
      }
    droprownum <- which(tchar %in% pchar)
    z_list[[i]] <- thisdf[-droprownum, ]
    }
  } else if(orientation == "bottom" | orientation == "west" | orientation == "south") {
    z_vec <- z_vec[order(z_vec, decreasing = F)]
    for (i in 2:length(z_vec)) {
      searchdf <- wanted_mat %>% filter(eval(parse(text = wanted_z)) 
                                        < as.numeric(z_vec[i]))
      thisdf <- wanted_mat %>% filter(eval(parse(text = wanted_z)) 
                                      == as.numeric(z_vec[i]))
      t <- searchdf %>% dplyr::select(wanted_x, wanted_y)
      p <- thisdf %>% dplyr::select(wanted_x, wanted_y)
      tchar <- vector(length = nrow(t))
      for (j in 1:nrow(t)) {
        tchar[j] <- paste0(t[j,1], t[j,2])
      }
      pchar <- vector(length = nrow(p))
      for (j in 1:nrow(p)) {
        pchar[j] <- paste0(p[j,1], p[j,2])
      }
      droprownum <- which(tchar %in% pchar)
      if(length(droprownum) == 0) {
        z_list[[i]] <- thisdf
      } else {
        z_list[[i]] <- thisdf[-droprownum, ]
      }
    }
  } else {
    print("not a usable orientation")
  }
  contour_df <- do.call(rbind, z_list)
  # code to fill in gaps between different ranged units in an axis for smooth plotting
  if(orientation == "east" | orientation == "west" | orientation == "south" | 
     orientation == "north") {
    cont_x <- (contour_df %>% dplyr::select(wanted_x))[,1]
    cont_x <- unique(cont_x[order(cont_x)])
    facx <- vector(length = length(cont_x))
    for(i in 1:length(cont_x)) {
      if(i == 1) {
        facx[i] <- abs(cont_x[i+1] - cont_x[i])
      } else {
        facx[i] <- abs(cont_x[i] - cont_x[i-1])
      }
    }
    facxu <- unique(facx)
    if(any(facxu > 5) & all(is_wholenumber(facxu))) {
      needed_slice <- cont_x[which(facx > 5)]
      extra_slices_list <- vector("list", length = (length(needed_slice) - 1))
      for(i in 1:(length(needed_slice) - 1)) {
        new_vec <- (needed_slice[i]:needed_slice[i+1])
        new_div <- new_vec/5
        new_fac <- new_vec[which(is_wholenumber(new_div))]
        new_fac_df_list <- vector("list", length = length(new_fac))
        for (j in 2:(length(new_fac) - 1)) {
          df <- contour_df %>% filter(eval(parse(text = wanted_x)) == 
                                        needed_slice[i+1])
          df[,which(colnames(df) == wanted_x)] <- rep(new_fac[j], times = nrow(df))
          new_fac_df_list[[j]] <- df
        }
        new_fac_df <- do.call(rbind, new_fac_df_list)
        extra_slices_list[[i]] <- new_fac_df
      }
      extra_slices <- do.call(rbind, extra_slices_list)
      contour_df <- rbind(contour_df, extra_slices)
    }
    cont_y <- (contour_df %>% dplyr::select(wanted_y))[,1]
    cont_y <- unique(cont_y[order(cont_y)])
    facy <- vector(length = length(cont_y))
    for(i in 1:length(cont_y)) {
      if(i == 1) {
        facy[i] <- abs(cont_y[i+1] - cont_y[i])
      } else {
        facy[i] <- abs(cont_y[i] - cont_y[i-1])
      }
    }
    facyu <- unique(facy)
    if(any(facyu > 5) & all(is_wholenumber(facyu))) {
      needed_slice <- cont_y[which(facy > 5)]
      extra_slices_list <- vector("list", length = (length(needed_slice) - 1))
      for(i in 1:(length(needed_slice) - 1)) {
        new_vec <- (needed_slice[i]:needed_slice[i+1])
        new_div <- new_vec/5
        new_fac <- new_vec[which(is_wholenumber(new_div))]
        new_fac_df_list <- vector("list", length = length(new_fac))
        for (j in 2:(length(new_fac) - 1)) {
          df <- contour_df %>% filter(eval(parse(text = wanted_y)) == 
                                        needed_slice[i+1])
          df[,which(colnames(df) == wanted_y)] <- rep(new_fac[j], times = nrow(df))
          new_fac_df_list[[j]] <- df
        }
        new_fac_df <- do.call(rbind, new_fac_df_list)
        extra_slices_list[[i]] <- new_fac_df
      }
      extra_slices <- do.call(rbind, extra_slices_list)
      contour_df <- rbind(contour_df, extra_slices)
    }
  }
  if(thresholded == T) {
    new_x <- contour_df %>% dplyr::select(wanted_x)
    new_y <- contour_df %>% dplyr::select(wanted_y)
    new_z <- contour_df %>% dplyr::select(wanted_z)
    thresh_val <- rep(1, times = nrow(contour_df))
    contour_final <- data.frame(new_x[,1], new_y[,1], new_z[,1], thresh_val)
    names(contour_final) <- c("x", "y", "z", "thresh_val")
  } else {
    new_x <- contour_df %>% dplyr::select(wanted_x)
    new_y <- contour_df %>% dplyr::select(wanted_y)
    new_z <- contour_df %>% dplyr::select(wanted_z)
    new_s <- contour_df %>% dplyr::select(wanted_val)
    contour_final <- data.frame(new_x[,1], new_y[,1], new_z[,1], new_s)
    names(contour_final) <- c("x", "y", "z", "suit_val")
  }
  return(contour_final)
}

# threshold top down
thresh_topdown <- contour_mat(wanted_mat = thresh_mat, wanted_z = "z", wanted_x = "x",
                              wanted_y = "y", orientation = "top", thresholded = T)
thresh_topdown_df <- thresh_topdown %>% dplyr::select(x, y, thresh_val)
thresh_topdown_rast <- rast(thresh_topdown_df, type = 'xyz', crs = target_crs)

pdf('./Plots/l_chalumnae_3D_contour_topdown.pdf')
p_top <- ggplot() +
  geom_spatraster(data = thresh_topdown_rast) +
  geom_sf(data = land) +
  coord_sf(xlim = c(9.25, 63.75), ylim = c(-48.75, 14.25), expand = F) +
  scale_fill_gradientn(colors = c("white", "darkgreen"), 
                       values = c(0, 1),
                       na.value = "transparent") +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.title = element_text(size = 15)) +
  xlab("Longitude") +
  ylab("Latitude")
print(p_top)
dev.off()

# threshold east of top
thresh_east <- contour_mat(wanted_mat = thresh_mat, wanted_z = 'x', wanted_x = "z",
                           wanted_y = "y", orientation = "east", thresholded = T)
thresh_east_df <- thresh_east %>% dplyr::select(x, y, thresh_val)

pdf('./Plots/l_chalumnae_3D_contour_eastoftop.pdf')
p_east <- ggplot(data = thresh_east_df, aes(x, y)) +
  geom_tile(aes(fill = thresh_val), colour = 'darkgreen') +
  scale_x_reverse() +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.title = element_text(size = 15)) +
  xlab("Depth (m)") +
  ylab("Latitude") +
  ylim(-48.75, 14.25)
print(p_east)
dev.off()

# threshold south of top
thresh_south <- contour_mat(wanted_mat = thresh_mat, wanted_z = "y", wanted_x = 'x',
                            wanted_y = "z", orientation = "south", thresholded = T)
thresh_south_df <- thresh_south %>% dplyr::select(x, y, thresh_val)

pdf('./Plots/l_chalumnae_3D_contour_southoftop.pdf')
p_south <- ggplot(data = thresh_south_df, aes(x, y)) +
  geom_tile(aes(fill = thresh_val), colour = 'darkgreen') +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.title = element_text(size = 15)) +
  xlab("Longitude") +
  ylab("Depth (m)") +
  xlim(9.25, 63.75)
print(p_south)
dev.off()

# full L chalumnae threshold figure
pdf('./Plots/Final_chalumnae_thresh_contour_all.pdf')
p_top <- ggplot() +
  geom_spatraster(data = thresh_topdown_rast) +
  geom_sf(data = land) +
  coord_sf(xlim = c(9.25, 63.75), ylim = c(-48.75, 14.25), expand = F) +
  scale_fill_gradientn(colors = c("white", "darkgreen"), 
                       values = c(0, 1),
                       na.value = "transparent") +
  theme(panel.grid.major = element_blank(), axis.title.x = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.title = element_text(size = 15)) +
  ylab("Latitude")
p_east <- ggplot(data = thresh_east_df, aes(x, y)) +
  geom_tile(aes(fill = thresh_val), colour = 'darkgreen') +
  scale_x_reverse() +
  theme(panel.grid.major = element_blank(), axis.title.y = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.title = element_text(size = 15)) +
  xlab("Depth (m)") +
  ylim(-48.75, 14.25)
grid.arrange(p_top, p_east, p_south, ncol = 2)
dev.off()

# suitability topdown
suit_topdown <- contour_mat(wanted_mat = suit_mat, wanted_z = "z", wanted_x = "x",
                              wanted_y = "y", orientation = "top", thresholded = F,
                            wanted_val = "s")
suit_topdown_df <- suit_topdown %>% dplyr::select(x, y, suit_val)
suit_topdown_rast <- rast(suit_topdown_df, type = 'xyz', crs = target_crs)

pdf('./Plots/l_chalumnae_3D_suit_contour_topdown.pdf')
p_top <- ggplot() +
  geom_spatraster(data = suit_topdown_rast) +
  geom_sf(data = land) +
  coord_sf(xlim = c(9.25, 63.75), ylim = c(-48.75, 14.25), expand = F) +
  scale_fill_gradientn(colors = c("white", viridis(4)), 
                       values = c(0, 0.25, 0.5, 0.75, 1),
                       na.value = "transparent",
                       name = "Suitability Score") +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"),
        axis.title = element_text(size = 15)) +
  xlab("Longitude") +
  ylab("Latitude")
print(p_top)
dev.off()

# suitability east of top
suit_east <- contour_mat(wanted_mat = suit_mat, wanted_z = 'x', wanted_x = "z",
                           wanted_y = "y", orientation = "east", thresholded = F,
                         wanted_val = "s")
suit_east_df <- suit_east %>% dplyr::select(x, y, suit_val)

pdf('./Plots/l_chalumnae_3D_suit_contour_eastoftop.pdf')
p_east <- ggplot(data = suit_east_df, aes(x, y, fill = suit_val)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("white", viridis(4)), 
                       values = c(0, 0.25, 0.5, 0.75, 1),
                       na.value = "transparent",
                       name = "Suitability Score") +
  scale_x_reverse() +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.title = element_text(size = 15)) +
  xlab("Depth (m)") +
  ylab("Latitude") +
  ylim(-48.75, 14.25)
print(p_east)
dev.off()

# suitability south of top
suit_south <- contour_mat(wanted_mat = suit_mat, wanted_z = "y", wanted_x = 'x',
                            wanted_y = "z", orientation = "south", thresholded = F,
                          wanted_val = "s")
suit_south_df <- suit_south %>% dplyr::select(x, y, suit_val)

pdf('./Plots/l_chalumnae_3D_suit_contour_southoftop.pdf')
p_south <- ggplot(data = suit_south_df, aes(x, y, fill = suit_val)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("white", viridis(4)), 
                       values = c(0, 0.25, 0.5, 0.75, 1),
                       na.value = "transparent",
                       name = "Suitability Score") +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.title = element_text(size = 15)) +
  xlab("Longitude") +
  ylab("Depth (m)") +
  xlim(9.25, 63.75)
print(p_south)
dev.off()

# final L chalumnae suitability contour figure
pdf('./Plots/Final_chalumnae_suitability_contour_figure.pdf')
p_top <- ggplot() +
  geom_spatraster(data = suit_topdown_rast, show.legend = F) +
  geom_sf(data = land) +
  coord_sf(xlim = c(9.25, 63.75), ylim = c(-48.75, 14.25), expand = F) +
  scale_fill_gradientn(colors = c("white", viridis(4)), 
                       values = c(0, 0.25, 0.5, 0.75, 1),
                       na.value = "transparent",
                       name = "Suitability Score") +
  theme(panel.grid.major = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title = element_text(size = 15)) +
  ylab("Latitude")
p_east <- ggplot(data = suit_east_df, aes(x, y, fill = suit_val)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("white", viridis(4)), 
                       values = c(0, 0.25, 0.5, 0.75, 1),
                       na.value = "transparent",
                       name = "Suitability Score") +
  scale_x_reverse() +
  theme(panel.grid.major = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.title = element_text(size = 15)) +
  xlab("Depth (m)") +
  ylim(-48.75, 14.25)
grid.arrange(p_top, p_east, p_south, ncol = 2)
dev.off()

# Indopacific Projection 2D
l_menadoensis_2D <- read.csv('./data/l_menadoensis_2D_occs.csv')

# suitability
new_ext <- ext(90, 150, -30, 30)
Indo_suit_2D <- rast('./Model_results/indopacific_projection_2D.tif')
Indo_suit_2D <- crop(Indo_suit_2D, new_ext)

tiff('./Plots/Indo_projection_2D_suit_final.tiff', width = 3400, height = 3400, res = 300)
p1 <- ggplot() +
  geom_spatraster(data = Indo_suit_2D) +
  geom_sf(data = land, fill = "darkgrey", linewidth = 0.6) +
  coord_sf(xlim = c(90, 150), ylim = c(-30, 30), expand = F) +
  scale_fill_gradientn(colors = c("white", viridis(4)), 
                       values = c(0, 0.25, 0.5, 0.75, 1),
                       na.value = "transparent",
                       name = "Suitability Score") +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 17), axis.text.y = element_text(size = 17),
        legend.text = element_text(size = 13), legend.title = element_text(size = 15),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black", linewidth = 0.85),
        legend.position = c(0.17, 0.87), title = element_text(size = 18)) +
  ggtitle("Raw Suitability")
print(p1)
dev.off()

# plotting 2D threshold
indo_thresh_2D <- rast('./Model_results/indopacific_threshold_2D.tif')
indo_thresh_2D_nm <- rast('./Model_results/indopacific_threshold_2D_notmasked.tif')
indo_thresh_2D_nm <- crop(indo_thresh_2D_nm, new_ext)

tiff('./Plots/Indo_projection_2D_thresh_final.tiff', width = 3400, height = 3400,
     res = 300)
p2 <- ggplot() +
  geom_spatraster(data = indo_thresh_2D_nm) +
  geom_sf(data = land, fill = "darkgrey", linewidth = 0.6) +
  coord_sf(xlim = c(90, 150), ylim = c(-30, 30), expand = F) +
  scale_fill_gradientn(colors = c("white", "darkgreen"), values = c(0,1),
                       na.value = "transparent") +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 17), axis.text.y = element_text(size = 17),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
        legend.position = "none", title = element_text(size = 18)) +
  ggtitle("Presence/Absence")
print(p2)
dev.off()

# indo projection 2D suitability and threshold together
tiff("./Plots/figure_5_real.tiff",
     width = 3400, height = 3400, res = 300)
grid.arrange(p1, p2, nrow = 1)
dev.off()

tiff('./Plots/Indo_projection_2D_thresh_final_zoom_points.tiff',
     width = 3400, height = 3400, res = 300)
ggplot() +
  geom_spatraster(data = indo_thresh_2D_nm) +
  geom_sf(data = land, fill = "darkgrey", linewidth = 0.65) +
  geom_point(data = l_menadoensis_2D, aes(x = Longitude, y = Latitude), size = 7,
             color = "#440154FF") +
  coord_sf(xlim = c(120, 137), ylim = c(-5, 5), expand = F) +
  scale_fill_gradientn(colors = c("white", "#03C04A"), values = c(0,1),
                       na.value = "transparent") +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.text = element_text(size = 20),
        axis.title = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
        title = element_text(size = 20)) +
  ggtitle("2D Thresholded Presence Between Sulawesi and Papua New Guinea")
dev.off()

# 3D Indo plots
# reading in env layer for depth list
conductivity_summer_18_indo <- rast('./Envs_Indo/conductivity_summer_18_indo.tif')
men_depths <- as.numeric(gsub("X", '', names(conductivity_summer_18_indo)))

# reading in suitability and threshold layers
indo_3D_suit <- rast('./Model_results/Indopacific_projection_3D.tif')
indo_3D_thresh <- rast('./Model_results/Indopacific_thresholded_3D.tif')
indo_3D_thresh_nm <- rast('./Model_results/Indopacific_thresholded_3D_notmasked.tif')
indo_3D_suit <- crop(indo_3D_suit, new_ext)
indo_3D_thresh <- crop(indo_3D_thresh, new_ext)
indo_3D_thresh_nm <- crop(indo_3D_thresh_nm, new_ext)

# plotlayers plot of suitability
names(indo_3D_thresh) <- men_depths

pdf('./Plots/indo_3D_plotlayers.pdf', pointsize = 3)
plotLayers(rast = indo_3D_thresh, land = land)
dev.off()

# suitability plot layer by layer
pdf('./Plots/indo_3D_suit_all_layers.pdf', pointsize = 3)
for(i in 1:dim(indo_3D_suit)[3]) {
  p <- ggplot() +
    geom_spatraster(data = indo_3D_suit[[i]]) +
    geom_sf(data = land) +
    coord_sf(xlim = c(90, 150), ylim = c(-30, 30), expand = F) +
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
pdf('./Plots/indo_3D_thresh_all_layers.pdf', pointsize = 3)
for(i in 1:dim(indo_3D_thresh_nm)[3]) {
  p <- ggplot() +
    geom_spatraster(data = indo_3D_thresh_nm[[i]]) +
    geom_sf(data = land) +
    coord_sf(xlim = c(90, 150), ylim = c(-30, 30), expand = F) +
    scale_fill_gradientn(colors = c("white", "darkgreen"), 
                         values = c(0, 1),
                         na.value = "transparent") +
    theme(panel.grid.major = element_blank(), 
          panel.background = element_rect(fill = "white"), legend.position = "none") +
    ggtitle(paste0("Thresholded Presence at ", men_depths[i], " m"))
  print(p)
}
dev.off()

# creating zoomed in plots of the layers where there are occurrences
men_3D_occs <- read.csv('./data/l_menadoensis_3D_occs.csv')
unique_depths <- unique(men_3D_occs$depth)
want_layers <- which(men_depths %in% unique_depths)

p_list <- vector("list", length = length(want_layers))
for (i in seq_along(want_layers)) {
  need_dat <- men_3D_occs %>% filter(depth == men_depths[want_layers[i]])
  p <- ggplot() +
    geom_spatraster(data = indo_3D_thresh_nm[[want_layers[i]]]) +
    geom_sf(data = land, fill = "darkgrey", linewidth = 0.6) +
    geom_point(data = need_dat, aes(x = Longitude, y = Latitude), size = 5,
               color = "#440154FF") +
    coord_sf(xlim = c(120, 137), ylim = c(-5, 5), expand = F) +
    scale_fill_gradientn(colors = c("white", "#03C04A"), values = c(0,1),
                         na.value = "transparent") +
    theme(panel.grid.major = element_blank(), 
          panel.background = element_rect(fill = "white"), legend.position = "none",
          axis.text = element_text(size = 15),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
          axis.title = element_blank(),
          title = element_text(size = 15)) +
    ggtitle(paste0("Thresholded Presence at ", men_depths[want_layers[i]], " m"))
  p_list[[i]] <- p
}

tiff('./Plots/indo_zoom_menadoensis_occ_plots.tiff', width = 3400, height = 3400, res = 300)
grid.arrange(p_list[[1]], p_list[[2]], p_list[[3]], p_list[[4]], p_list[[5]], p_list[[6]],
             ncol = 2)
dev.off()

# converting rasters to xyz matrices
indo_suit_mat <- rast_to_mat(indo_3D_suit, men_depths, thresholded = F)
indo_thresh_mat <- rast_to_mat(indo_3D_thresh, men_depths, thresholded = T)

# creating contour plots
# threshold top down
thresh_topdown_indo <- contour_mat(wanted_mat = indo_thresh_mat, wanted_z = "z", 
                              wanted_x = "x",
                              wanted_y = "y", orientation = "top", thresholded = T)
thresh_topdown_df_indo <- thresh_topdown_indo %>% dplyr::select(x, y, thresh_val)
thresh_topdown_rast_indo <- rast(thresh_topdown_df_indo, type = 'xyz', 
                                 crs = target_crs)

pdf('./Plots/Indonesia_3D_contour_topdown.pdf')
p_top <- ggplot() +
  geom_spatraster(data = thresh_topdown_rast_indo) +
  geom_sf(data = land) +
  coord_sf(xlim = c(90, 150), ylim = c(-30, 30), expand = F) +
  scale_fill_gradientn(colors = c("white", "darkgreen"), 
                       values = c(0, 1),
                       na.value = "transparent") +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.title = element_text(size = 15)) +
  xlab("Longitude") +
  ylab("Latitude")
print(p_top)
dev.off()

# threshold east of top
thresh_east_indo <- contour_mat(wanted_mat = indo_thresh_mat, wanted_z = 'x', 
                                wanted_x = "z",
                           wanted_y = "y", orientation = "east", thresholded = T)
thresh_east_df_indo <- thresh_east_indo %>% dplyr::select(x, y, thresh_val)

pdf('./Plots/Indonesia_3D_contour_eastoftop.pdf')
p_east <- ggplot(data = thresh_east_df_indo, aes(x, y)) +
  geom_tile(aes(fill = thresh_val), colour = 'darkgreen') +
  scale_x_reverse() +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.title = element_text(size = 15)) +
  xlab("Depth (m)") +
  ylab("Latitude") +
  ylim(-30, 30)
print(p_east)
dev.off()

# threshold south of top
thresh_south_indo <- contour_mat(wanted_mat = indo_thresh_mat, wanted_z = "y", 
                                 wanted_x = 'x',
                            wanted_y = "z", orientation = "south", thresholded = T)
thresh_south_df_indo <- thresh_south_indo %>% dplyr::select(x, y, thresh_val)

pdf('./Plots/Indonesia_3D_contour_southoftop.pdf')
p_south <- ggplot(data = thresh_south_df_indo, aes(x, y)) +
  geom_tile(aes(fill = thresh_val), colour = 'darkgreen') +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.title = element_text(size = 15)) +
  xlab("Longitude") +
  ylab("Depth (m)") +
  xlim(90, 150)
print(p_south)
dev.off()

# Final chalumnae projected 3D threshold contour figure
pdf('./Plots/Final_chalumnae_projected_indo_thresh_3D.pdf')
p_top <- ggplot() +
  geom_spatraster(data = thresh_topdown_rast_indo) +
  geom_sf(data = land) +
  coord_sf(xlim = c(90, 150), ylim = c(-30, 30), expand = F) +
  scale_fill_gradientn(colors = c("white", "darkgreen"), 
                       values = c(0, 1),
                       na.value = "transparent") +
  theme(panel.grid.major = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.title = element_text(size = 15)) +
  ylab("Latitude")
p_east <- ggplot(data = thresh_east_df_indo, aes(x, y)) +
  geom_tile(aes(fill = thresh_val), colour = 'darkgreen') +
  scale_x_reverse() +
  theme(panel.grid.major = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.title = element_text(size = 15)) +
  xlab("Depth (m)") +
  ylim(-30, 30)
grid.arrange(p_top, p_east, p_south, ncol = 2)
dev.off()

# suitability topdown
suit_topdown_indo <- contour_mat(wanted_mat = indo_suit_mat, wanted_z = "z", 
                                 wanted_x = "x",
                            wanted_y = "y", orientation = "top", thresholded = F,
                            wanted_val = "s")
suit_topdown_df_indo <- suit_topdown_indo %>% dplyr::select(x, y, suit_val)
suit_topdown_rast_indo <- rast(suit_topdown_df_indo, type = 'xyz', crs = target_crs)

pdf('./Plots/Indonesia_3D_suit_contour_topdown.pdf')
p_top <- ggplot() +
  geom_spatraster(data = suit_topdown_rast_indo) +
  geom_sf(data = land) +
  coord_sf(xlim = c(90, 150), ylim = c(-30, 30), expand = F) +
  scale_fill_gradientn(colors = c("white", viridis(4)), 
                       values = c(0, 0.25, 0.5, 0.75, 1),
                       na.value = "transparent",
                       name = "Suitability Score") +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"),
        axis.title = element_text(size = 15)) +
  xlab("Longitude") +
  ylab("Latitude")
print(p_top)
dev.off()

# suitability east of top
suit_east_indo <- contour_mat(wanted_mat = indo_suit_mat, wanted_z = 'x', 
                              wanted_x = "z",
                         wanted_y = "y", orientation = "east", thresholded = F,
                         wanted_val = "s")
suit_east_df_indo <- suit_east_indo %>% dplyr::select(x, y, suit_val)

pdf('./Plots/Indonesia_3D_suit_contour_eastoftop.pdf')
p_east <- ggplot(data = suit_east_df_indo, aes(x, y, fill = suit_val)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("white", viridis(4)), 
                       values = c(0, 0.25, 0.5, 0.75, 1),
                       na.value = "transparent",
                       name = "Suitability Score") +
  scale_x_reverse() +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.title = element_text(size = 15)) +
  xlab("Depth (m)") +
  ylab("Latitude") +
  ylim(-30, 30)
print(p_east)
dev.off()

# suitability south of top
suit_south_indo <- contour_mat(wanted_mat = indo_suit_mat, wanted_z = "y", 
                               wanted_x = 'x',
                          wanted_y = "z", orientation = "south", thresholded = F,
                          wanted_val = "s")
suit_south_df_indo <- suit_south_indo %>% dplyr::select(x, y, suit_val)

pdf('./Plots/Indonesia_3D_suit_contour_southoftop.pdf')
p_south <- ggplot(data = suit_south_df_indo, aes(x, y, fill = suit_val)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("white", viridis(4)), 
                       values = c(0, 0.25, 0.5, 0.75, 1),
                       na.value = "transparent",
                       name = "Suitability Score") +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.title = element_text(size = 15)) +
  xlab("Longitude") +
  ylab("Depth (m)") +
  xlim(90, 150)
print(p_south)
dev.off()

# final chalumnae projected suitability figure
pdf('./Plots/Final_chalumnae_projected_indo_suit_3D.pdf')
p_top <- ggplot() +
  geom_spatraster(data = suit_topdown_rast_indo, show.legend = F) +
  geom_sf(data = land) +
  coord_sf(xlim = c(90, 150), ylim = c(-30, 30), expand = F) +
  scale_fill_gradientn(colors = c("white", viridis(4)), 
                       values = c(0, 0.25, 0.5, 0.75, 1),
                       na.value = "transparent",
                       name = "Suitability Score") +
  theme(panel.grid.major = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title = element_text(size = 15)) +
  ylab("Latitude")
p_east <- ggplot(data = suit_east_df_indo, aes(x, y, fill = suit_val)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("white", viridis(4)), 
                       values = c(0, 0.25, 0.5, 0.75, 1),
                       na.value = "transparent",
                       name = "Suitability Score") +
  scale_x_reverse() +
  theme(panel.grid.major = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        panel.background = element_rect(fill = "white"), legend.position = "none",
        axis.title = element_text(size = 15)) +
  xlab("Depth (m)") +
  ylim(-30, 30)
p <- grid.arrange(p_top, p_east, p_south, ncol = 2)
print(p)
dev.off()

