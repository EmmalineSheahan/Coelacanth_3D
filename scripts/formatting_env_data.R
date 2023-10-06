# preparing environmental data for 3D and 2D modeling

library(sf)
library(terra)
library(dplyr)
library(voluModel)
library(R.utils)
library(rnaturalearth)
library(rgdal)

# creating wanted crs and land polygon
wanted_crs <- make_EPSG() %>% filter(code == 4326)
target_crs <- wanted_crs$prj4

land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]

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
