# preparing environmental data for 3D and 2D modeling

library(sf)
library(terra)
library(dplyr)
library(voluModel)
library(R.utils)
library(rnaturalearth)
library(rnaturalearthhires)
library(rgdal)
library(ncdf4)
library(raster)
library(rgeos)

# creating wanted crs and land polygon
wanted_crs <- make_EPSG() %>% filter(code == 4326)
target_crs <- wanted_crs$prj4
cea_crs <- "+proj=cea +lat_ts=0 +lon_0"

land <- st_as_sf(rnaturalearthhires::countries10)[1] %>% st_set_crs(target_crs)
land_poly <- as_Spatial(land)
land_poly_cea <- spTransform(land_poly, CRSobj = cea_crs)

# reading in occurrences for extent matching
# Indian ocean
l_chalumnae_down <- read.csv('./data/coelacanth_complete_occurrence.csv')
l_chalumnae_down <- l_chalumnae_down %>% 
  filter(Species == "Latimeria chalumnae") %>% 
  dplyr::select(c("Longitude", "Latitude"))
colnames(l_chalumnae_down) <- c("longitude", "latitude")

# Indo-Pacific
l_menadoensis <- read.csv('./data/l_menadoensis_occurrence.csv')
l_menadoensis <- l_menadoensis %>% dplyr::select(c("Longitude", "Latitude"))
l_menadoensis <- l_menadoensis[(which(!(duplicated(l_menadoensis)))),]
colnames(l_menadoensis) <- c("longitude", "latitude")

# creating base raster for standardized extent, cropped to the extent 
# of a shape buffered for all of the points
# Indian Ocean
temp_alph <- marineBackground(l_chalumnae_down, fraction = 1, 
                              partCount = 2, 
                              clipToOcean = T,
                              buff = 600000)
temp_alph <- spTransform(as_Spatial(st_as_sf(temp_alph)), 
                         CRSobj = cea_crs)

# Indo-Pacific
temp_alph_men <- marineBackground(l_menadoensis, fraction = 1, 
                              partCount = 2, 
                              clipToOcean = T,
                              buff = 3500000)
temp_alph_men <- spTransform(as_Spatial(st_as_sf(temp_alph_men)), 
                         CRSobj = cea_crs)

# buffering shape by a fact of 6 times the area covered by the points
# Indian Ocean
buffDist <- (sqrt(6*gArea(temp_alph)) - sqrt(gArea(temp_alph)))/2
shape_new <- raster::buffer(x = temp_alph, width = buffDist, dissolve = T)
shape_new_clipped <- gDifference(shape_new, land_poly_cea)
shape_sf <- st_as_sf(shape_new_clipped)
shape_final <- st_transform(shape_sf, target_crs)
wanted_extent <- ext(shape_final)

# Indo-Pacific
buffDist <- (sqrt(2*gArea(temp_alph_men)) - sqrt(gArea(temp_alph_men)))/2
shape_new <- raster::buffer(x = temp_alph_men, width = buffDist, dissolve = T)
shape_new_clipped <- gDifference(shape_new, land_poly_cea)
shape_sf <- st_as_sf(shape_new_clipped)
shape_final <- st_transform(shape_sf, target_crs)
wanted_extent_men <- ext(shape_final)

# creating raster extent
# Indian Ocean
gunzip('./WOA_18/woa18_decav_t13mn04.csv.gz')
temp_csv <- read.csv('./WOA_18/woa18_decav_t13mn04.csv', header = F)
depth_vals <- c(0, matrix(unlist(temp_csv[2, 4:ncol(temp_csv)]))[,1])
for (i in 1:length(depth_vals)) {
  depth_vals[i] <- paste0("X", depth_vals[i])
}
temp_csv <- temp_csv[-(1:2),]
new_colnames <- c("Latitude", "Longitude", depth_vals)
colnames(temp_csv) <- new_colnames
temp_df <- temp_csv %>% dplyr::select(Longitude, Latitude, depth_vals[11])
temp_rast_ext <- rast(temp_df, type = "xyz", crs = target_crs, 
                      extent = wanted_extent)

# Indo-Pacific
temp_csv <- read.csv('./WOA_18/woa18_decav_t13mn04.csv', header = F)
depth_vals <- c(0, matrix(unlist(temp_csv[2, 4:ncol(temp_csv)]))[,1])
for (i in 1:length(depth_vals)) {
  depth_vals[i] <- paste0("X", depth_vals[i])
}
temp_csv <- temp_csv[-(1:2),]
new_colnames <- c("Latitude", "Longitude", depth_vals)
colnames(temp_csv) <- new_colnames
temp_df <- temp_csv %>% dplyr::select(Longitude, Latitude, depth_vals[11])
temp_rast_ext_men <- rast(temp_df, type = "xyz", crs = target_crs, 
                      extent = wanted_extent_men)

# env_maker function
# wanted_data = filepath of csv with env data
# output_name = filepath of output tif
# wanted_extent = raster with desired extent to crop to
# converts WOA csv to spatraster stack of depthslices
env_maker <- function(wanted_data, output_name, wanted_extent) {
  temp_csv <- read.csv(wanted_data, header = F)
  depth_vals <- c(0, matrix(unlist(temp_csv[2, 4:ncol(temp_csv)]))[,1])
  for (i in 1:length(depth_vals)) {
    depth_vals[i] <- paste0("X", depth_vals[i])
  }
  temp_csv <- temp_csv[-(1:2),]
  new_colnames <- c("Latitude", "Longitude", depth_vals)
  colnames(temp_csv) <- new_colnames
  rast_list <- vector("list", length = length(depth_vals))
  for (i in 1:length(depth_vals)) {
    temp_df <- temp_csv %>% dplyr::select(Longitude, Latitude, depth_vals[i])
    temp_rast <- rast(temp_df, type = "xyz", crs = target_crs, 
                      extent = ext(wanted_extent))
    rast_list[[i]] <- temp_rast
  }
  final_rast <- rast(rast_list)
  final_rast <- final_rast[[c(1:43)]]
  writeRaster(final_rast, filename = output_name, overwrite = T)
}

# Temperature
gunzip('./WOA_18/woa18_decav_t15mn04.csv.gz')

env_maker(wanted_data = './WOA_18/woa18_decav_t13mn04.csv', 
          output_name = './WOA_18/temperature_winter_18.tif',
          wanted_extent = temp_rast_ext)
env_maker(wanted_data = './WOA_18/woa18_decav_t15mn04.csv',
          output_name = './WOA_18/temperature_summer_18.tif',
          wanted_extent = temp_rast_ext)

env_maker(wanted_data = './WOA_18/woa18_decav_t15mn04.csv',
          output_name = './Envs_Indo/temperature_summer_18_indo.tif',
          wanted_extent = temp_rast_ext_men)

# Salinity
gunzip('./WOA_18/woa18_decav_s13mn04.csv.gz')
gunzip('./WOA_18/woa18_decav_s15mn04.csv.gz')

env_maker(wanted_data = './WOA_18/woa18_decav_s13mn04.csv', 
          output_name = './WOA_18/salinity_winter_18.tif',
          wanted_extent = temp_rast_ext)
env_maker(wanted_data = './WOA_18/woa18_decav_s15mn04.csv',
          output_name = './WOA_18/salinity_summer_18.tif',
          wanted_extent = temp_rast_ext)

env_maker(wanted_data = './WOA_18/woa18_decav_s13mn04.csv',
          output_name = './Envs_Indo/salinity_winter_18_indo.tif',
          wanted_extent = temp_rast_ext_men)

# Density
gunzip('./WOA_18/woa18_decav_I13mn04.csv.gz')
gunzip('./WOA_18/woa18_decav_I15mn04.csv.gz')

env_maker(wanted_data = './WOA_18/woa18_decav_I13mn04.csv', 
          output_name = './WOA_18/density_winter_18.tif',
          wanted_extent = temp_rast_ext)
env_maker(wanted_data = './WOA_18/woa18_decav_I15mn04.csv',
          output_name = './WOA_18/density_summer_18.tif',
          wanted_extent = temp_rast_ext)

env_maker(wanted_data = './WOA_18/woa18_decav_I15mn04.csv',
          output_name = './Envs_Indo/density_summer_18_indo.tif',
          wanted_extent = temp_rast_ext_men)

# Conductivity
gunzip('./WOA_18/woa18_A5B7_C13mn04.csv.gz')
gunzip('./WOA_18/woa18_A5B7_C15mn04.csv.gz')

env_maker(wanted_data = './WOA_18/woa18_A5B7_C13mn04.csv', 
          output_name = './WOA_18/conductivity_winter_18.tif',
          wanted_ext = temp_rast_ext)
env_maker(wanted_data = './WOA_18/woa18_A5B7_C15mn04.csv',
          output_name = './WOA_18/conductivity_summer_18.tif',
          wanted_extent = temp_rast_ext)

env_maker(wanted_data = './WOA_18/woa18_A5B7_C13mn04.csv', 
          output_name = './Envs_Indo/conductivity_winter_18_indo.tif',
          wanted_ext = temp_rast_ext_men)
env_maker(wanted_data = './WOA_18/woa18_A5B7_C15mn04.csv', 
          output_name = './Envs_Indo/conductivity_summer_18_indo.tif',
          wanted_ext = temp_rast_ext_men)

# creating a SpatRaster of slopes calculated from the bathymetry layer subset by 
# depth slice
# Indian Ocean
temperature <- rast('./WOA_18/temperature_summer_18.tif')
depth_slices <- names(temperature)
depth_slices <- as.numeric(gsub('X', '', depth_slices))

bath <- rast('./Slope/ETOPO1_Bed_c_geotiff.tif')
values(bath)[which(values(bath) > 0)] <- NA
bath_slope <- rast(terrain(raster(bath), opt = "slope", unit = "degrees", 
                      neighbors = 4))
bath_slope <- terra::project(bath_slope, bath)
values(bath) <- abs(values(bath))

slope_list <- vector("list", length = length(depth_slices))
for (i in 1:length(depth_slices)) {
  if(i < length(depth_slices)) {
    temp_ras <- bath_slope
    if(length(which(values(bath) < depth_slices[i])) > 0) {
      values(temp_ras)[which(values(bath) < depth_slices[i])] <- 0
    }
    values(temp_ras)[which(values(bath) >= depth_slices[i+1])] <- 0
  }
  if(i == length(depth_slices)) {
    temp_ras <- bath_slope 
    values(temp_ras)[which(values(bath) >= (depth_slices[i] + 50))] <- 0
    values(temp_ras)[which(values(bath) < depth_slices[i])] <- 0
  }
  temp_ras <- terra::project(temp_ras, temperature)
  slope_list[[i]] <- temp_ras
}

slope_stack <- rast(slope_list)
names(slope_stack) <- names(temperature)
writeRaster(slope_stack, filename = './Slope/slope.tif', overwrite = T)

# Indo-Pacific
temperature_indo <- rast('./Envs_Indo/temperature_summer_18_indo.tif')
depth_slices_indo <- names(temperature_indo)
depth_slices_indo <- as.numeric(gsub('X', '', depth_slices_indo))

bath <- rast('./Slope/ETOPO1_Bed_c_geotiff.tif')
values(bath)[which(values(bath) > 0)] <- NA
bath_slope <- rast(terrain(raster(bath), opt = "slope", unit = "degrees", 
                           neighbors = 4))
bath_slope <- terra::project(bath_slope, bath)
values(bath) <- abs(values(bath))

slope_list <- vector("list", length = length(depth_slices_indo))
for (i in 1:length(depth_slices_indo)) {
  if(i < length(depth_slices_indo)) {
    temp_ras <- bath_slope
    if(length(which(values(bath) < depth_slices_indo[i])) > 0) {
      values(temp_ras)[which(values(bath) < depth_slices_indo[i])] <- 0
    }
    values(temp_ras)[which(values(bath) >= depth_slices_indo[i+1])] <- 0
  }
  if(i == length(depth_slices_indo)) {
    temp_ras <- bath_slope 
    values(temp_ras)[which(values(bath) >= (depth_slices_indo[i] + 50))] <- 0
    values(temp_ras)[which(values(bath) < depth_slices_indo[i])] <- 0
  }
  temp_ras <- terra::project(temp_ras, temperature_indo)
  slope_list[[i]] <- temp_ras
}

slope_stack <- rast(slope_list)
names(slope_stack) <- names(temperature_indo)
writeRaster(slope_stack, filename = './Envs_indo/slope_indo.tif', overwrite = T)

# Currents
# pull_current function to generate SpatRaster stacks of currents by depth
# nc_filename = name of the netcdf containing data
# txt_filename = name of metadata output text file of netcdf
# wanted_direction = "uo_mean" for East and "vo_mean" for North
# raster_name = final spatraster filename
# wanted_extent = raster with desired extent to crop to
pull_currents <- function(nc_filename, txt_filename, 
                          wanted_direction, raster_name, wanted_extent) {
  nc_current <- nc_open(paste0('./Currents/', nc_filename))
  {
    sink(paste0('./Currents/', txt_filename))
    print(nc_current)
    sink()
  }
  lon_current <- ncvar_get(nc_current, "longitude")
  lat_current <- ncvar_get(nc_current, "latitude")
  depth_current <- ncvar_get(nc_current, "depth")
  speed_current <- ncvar_get(nc_current, wanted_direction)
  
  depth_slices <- names(wanted_extent)
  depth_slices <- as.numeric(gsub('X', '', depth_slices))

  current_list <- vector("list", length = length(depth_slices))
  for (i in 1:length(depth_slices)) {
    if(i < (length(depth_slices) - 1)) {
      lower_bound <- which(depth_current > depth_slices[i])
      upper_bound <- which(depth_current[lower_bound] < depth_slices[i+2])
    } else {
      lower_bound <- which(depth_current > depth_slices[i])
      upper_bound <- which(depth_current[lower_bound] < (depth_slices[i] + 100))
    }
    if (length(upper_bound) == 0) {
      whole_depth_slice <- wanted_extent[[1]]
      values(whole_depth_slice) <- NA
      names(whole_depth_slice) <- depth_slices[i]
      current_list[[i]] <- whole_depth_slice
    } else {
      depth_slice_list <- vector("list", length = length(upper_bound))
      for (j in seq_along(upper_bound)) {
        current_slice <- speed_current[, , upper_bound[j]]
        r <- raster(t(current_slice), xmn=min(lon_current), 
                xmx=max(lon_current), 
                ymn=min(lat_current), 
                ymx=max(lat_current), 
                crs=CRS(target_crs))
        r_flip <- flip(r, direction = 'y')
        r_new <- rast(r_flip)
        r_final <- terra::project(r_new, wanted_extent)
        depth_slice_list[[j]] <- r_final
    }
    single_slice <- rast(depth_slice_list)
    whole_depth_slice <- terra::app(single_slice, mean)
    names(whole_depth_slice) <- depth_slices[i]
    current_list[[i]] <- whole_depth_slice
    }
  }
  final_current <- rast(current_list)
  names(final_current) <- names(wanted_extent)
  writeRaster(final_current, filename = raster_name,
              overwrite = T)
}

# Indian Ocean
pull_currents("grepv1_gl1_uv_202006.nc", "current_june_east.txt", 
              "uo_mean", "./Currents/current_east_june.tif", temperature)

pull_currents("grepv1_gl1_uv_202006.nc", "current_june_north.txt", 
              "vo_mean", "./Currents/current_north_june.tif", temperature)

pull_currents("grepv1_gl1_uv_202012.nc", "current_december_east.txt",
              "uo_mean", "./Currents/current_east_december.tif", temperature)

pull_currents("grepv1_gl1_uv_202012.nc", "current_december_north.txt",
              "vo_mean", "./Currents/current_north_december.tif", temperature)

# Indo-Pacific
pull_currents("grepv1_gl1_uv_202006.nc", "current_june_east.txt", 
              "uo_mean", "./Envs_Indo/current_east_june_indo.tif", temperature_indo)

pull_currents("grepv1_gl1_uv_202012.nc", "current_december_east.txt",
              "uo_mean", "./Envs_Indo/current_east_december_indo.tif", 
              temperature_indo)

pull_currents("grepv1_gl1_uv_202012.nc", "current_june_north.txt",
              "vo_mean", "./Envs_Indo/current_north_june_indo.tif", 
              temperature_indo)

pull_currents("grepv1_gl1_uv_202012.nc", "current_december_north.txt",
              "vo_mean", "./Envs_Indo/current_north_december_indo.tif", 
              temperature_indo)

# creating 2D rasters by averaging across depths and stacking into a single
# SpatRaster for 2D modelling

# env_flatten function creates a mean raster using the 3 dimmensional tif
# wanted_3D_tif = file path to 3D tif
# env_name = name of environmental variable
env_flatten <- function(wanted_3D_tif, env_name) {
  new_ras <- rast(wanted_3D_tif)
  mean_ras <- terra::app(new_ras, fun = 'mean', na.rm = T)
  mean_ras <- mask(mean_ras, land, inverse = T)
  names(mean_ras) <- env_name
  return(mean_ras)
}

# Indian Ocean
density_summer <- env_flatten('./WOA_18/density_summer_18.tif', "density_summer")
density_winter <- env_flatten('./WOA_18/density_winter_18.tif', "density_winter")
conductivity_summer <- env_flatten('./WOA_18/conductivity_summer_18.tif', 
                              "conductivity_summer")
conductivity_winter <- env_flatten('./WOA_18/conductivity_winter_18.tif', 
                              "conductivity_winter")
salinity_summer <- env_flatten('./WOA_18/salinity_summer_18.tif', 
                               "salinity_summer")
salinity_winter <- env_flatten('./WOA_18/salinity_winter_18.tif', 
                               "salinity_winter")
temperature_summer <- env_flatten('./WOA_18/temperature_summer_18.tif', 
                                  "temperature_summer")
temperature_winter <- env_flatten('./WOA_18/temperature_winter_18.tif', 
                                  "temperature_winter")
current_east_june <- env_flatten('./Currents/current_east_june.tif', 
                                 "current_east_june.tif")
current_north_june <- env_flatten('./Currents/current_north_june.tif', 
                                 "current_north_june.tif")
current_east_december <- env_flatten('./Currents/current_east_december.tif', 
                                 "current_east_december.tif")
current_north_december <- env_flatten('./Currents/current_north_december.tif', 
                                 "current_north_december.tif")
slope_mean <- env_flatten('./Slope/slope.tif', "slope")

all_2D_env <- c(density_summer, density_winter, conductivity_summer, 
                conductivity_winter,
                salinity_summer, salinity_winter, temperature_summer,
                temperature_winter, current_east_june, current_north_june,
                current_east_december, current_north_december, slope_mean)
names(all_2D_env) <- gsub('.tif', '', names(all_2D_env))
writeRaster(all_2D_env, filename = './Envs_2D/all_2D_env.tif', overwrite = T)

# Indo-Pacific
temperature_summer_indo <- env_flatten('./Envs_Indo/temperature_summer_18_indo.tif',
                                       "temperature_summer")
density_summer_indo <- env_flatten('./Envs_Indo/density_summer_18_indo.tif',
                                   "density_summer")
conductivity_summer_indo <- env_flatten('./Envs_Indo/conductivity_summer_18_indo.tif',
                                        "conductivity_summer")
conductivity_winter_indo <- env_flatten('./Envs_Indo/conductivity_winter_18_indo.tif',
                                        "conductivity_winter")
salinity_winter_indo <- env_flatten('./Envs_Indo/salinity_winter_18_indo.tif',
                                    "salinity_winter")
current_east_june_indo <- env_flatten('./Envs_Indo/current_east_june_indo.tif',
                                      "current_east_june")
current_east_december_indo <- env_flatten('./Envs_Indo/current_east_december_indo.tif',
                                      "current_east_december")
current_north_june_indo <- env_flatten('./Envs_Indo/current_north_june_indo.tif',
                                      "current_north_june")
current_north_december_indo <- 
  env_flatten('./Envs_Indo/current_north_december_indo.tif',
                                          "current_north_december")
slope_mean_indo <- env_flatten('./Envs_Indo/slope_indo.tif', "slope")

all_2D_env <- c(density_summer_indo, conductivity_summer_indo, 
                conductivity_winter_indo, current_east_june_indo,
                current_north_june_indo, current_east_december_indo,
                current_north_december_indo, slope_mean_indo)
names(all_2D_env) <- gsub('.tif', '', names(all_2D_env))
names(all_2D_env) <- gsub('_indo', '', names(all_2D_env))
writeRaster(all_2D_env, filename = './Envs_Indo/all_2D_env_indo.tif', overwrite = T)
