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

# creating wanted crs and land polygon
wanted_crs <- make_EPSG() %>% filter(code == 4326)
target_crs <- wanted_crs$prj4

land <- st_as_sf(rnaturalearthhires::countries10)[1] %>% st_set_crs(target_crs)

# env_maker function
# wanted_data = filepath of csv with env data
# output_name = filepath of output tif
# converts WOA csv to spatraster stack of depthslices
env_maker <- function(wanted_data, output_name) {
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
    temp_rast <- rast(temp_df, type = "xyz", crs = target_crs)
    rast_list[[i]] <- temp_rast
  }
  final_rast <- rast(rast_list)
  writeRaster(final_rast, filename = output_name, overwrite = T)
}

# Temperature
gunzip('./WOA_18/woa18_decav_t13mn01.csv.gz')
gunzip('./WOA_18/woa18_decav_t15mn01.csv.gz')

env_maker(wanted_data = './WOA_18/woa18_decav_t13mn01.csv', 
          output_name = './WOA_18/temperature_winter_18.tif')
env_maker(wanted_data = './WOA_18/woa18_decav_t15mn01.csv',
          output_name = './WOA_18/temperature_summer_18.tif')

# Salinity
gunzip('./WOA_18/woa18_decav_s13mn01.csv.gz')
gunzip('./WOA_18/woa18_decav_s15mn01.csv.gz')

env_maker(wanted_data = './WOA_18/woa18_decav_s13mn01.csv', 
          output_name = './WOA_18/salinity_winter_18.tif')
env_maker(wanted_data = './WOA_18/woa18_decav_s15mn01.csv',
          output_name = './WOA_18/salinity_summer_18.tif')

# Density
gunzip('./WOA_18/woa18_decav_I13mn01.csv.gz')
gunzip('./WOA_18/woa18_decav_I15mn01.csv.gz')

env_maker(wanted_data = './WOA_18/woa18_decav_I13mn01.csv', 
          output_name = './WOA_18/density_winter_18.tif')
env_maker(wanted_data = './WOA_18/woa18_decav_I15mn01.csv',
          output_name = './WOA_18/density_summer_18.tif')

# Oxygen
gunzip('./WOA_18/woa18_all_o13mn01.csv.gz')
gunzip('./WOA_18/woa18_all_o15mn01.csv.gz')

env_maker(wanted_data = './WOA_18/woa18_all_o13mn01.csv', 
          output_name = './WOA_18/oxygen_winter_18.tif')
env_maker(wanted_data = './WOA_18/woa18_all_o15mn01.csv',
          output_name = './WOA_18/oxygen_summer_18.tif')

# Phosphate
gunzip('./WOA_18/woa18_all_p13mn01.csv.gz')
gunzip('./WOA_18/woa18_all_p15mn01.csv.gz')

env_maker(wanted_data = './WOA_18/woa18_all_p13mn01.csv', 
          output_name = './WOA_18/phosphate_winter_18.tif')
env_maker(wanted_data = './WOA_18/woa18_all_p15mn01.csv',
          output_name = './WOA_18/phosphate_summer_18.tif')

# Nitrate
gunzip('./WOA_18/woa18_all_n13mn01.csv.gz')
gunzip('./WOA_18/woa18_all_n15mn01.csv.gz')

env_maker(wanted_data = './WOA_18/woa18_all_n13mn01.csv', 
          output_name = './WOA_18/nitrate_winter_18.tif')
env_maker(wanted_data = './WOA_18/woa18_all_n15mn01.csv',
          output_name = './WOA_18/nitrate_summer_18.tif')

# creating a SpatRaster of slopes calculated from the bathymetry layer subset by 
# depth slice
nitrate <- rast('./WOA_18/nitrate_summer_18.tif')
depth_slices <- names(nitrate)
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
  temp_ras <- terra::project(temp_ras, nitrate)
  slope_list[[i]] <- temp_ras
}

slope_stack <- rast(slope_list)
writeRaster(slope_stack, filename = './Slope/slope.tif')

# Currents
# pull_current function to generate SpatRaster stacks of currents by depth
# nc_filename = name of the netcdf containing data
# txt_filename = name of metadata output text file of netcdf
# wanted_direction = "uo_mean" for East and "vo_mean" for North
# raster_name = final spatraster filename
pull_currents <- function(nc_filename, txt_filename, 
                          wanted_direction, raster_name) {
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

  current_list <- vector("list", length = length(depth_slices))
  for (i in 1:length(depth_slices)) {
    if(i < length(depth_slices)) {
      lower_bound <- which(depth_current > depth_slices[i])
      upper_bound <- which(depth_current[lower_bound] < depth_slices[i+1])
    } else {
      lower_bound <- which(depth_current > depth_slices[i])
      upper_bound <- which(depth_current[lower_bound] < (depth_slices[i] + 50))
    }
    if (length(upper_bound) == 0) {
      whole_depth_slice <- nitrate[[1]]
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
        r_final <- terra::project(r_new, nitrate)
        depth_slice_list[[j]] <- r_final
    }
    single_slice <- rast(depth_slice_list)
    whole_depth_slice <- terra::app(single_slice, mean)
    names(whole_depth_slice) <- depth_slices[i]
    current_list[[i]] <- whole_depth_slice
    }
  }
  final_current <- rast(current_list)
  writeRaster(final_current, filename = paste0('./Currents/', raster_name))
}

pull_currents("grepv1_gl1_uv_202006.nc", "current_june_east.txt", 
              "uo_mean", "current_east_june.tif")

pull_currents("grepv1_gl1_uv_202006.nc", "current_june_north.txt", 
              "vo_mean", "current_north_june.tif")

pull_currents("grepv1_gl1_uv_202012.nc", "current_december_east.txt",
              "uo_mean", "current_east_december.tif")

pull_currents("grepv1_gl1_uv_202012.nc", "current_december_north.txt",
              "vo_mean", "current_north_december.tif")
