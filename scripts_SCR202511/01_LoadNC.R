
library(tidyverse)
library(ncdf4)
# library(cmocean)
# library(gganimate)
# library(sf)
# library(terra)

# Task: ensembling the models using the mean values for each cell/month/year for each of the SSPs

# nc_file <- nc_open("data/othr_blackcorals_2015_2099_cmip22_sorall_identity.nc")  # old version with smaller domain
nc_file <- nc_open("data/nafo_fishingfoot_2015_2099_cmip22_sorall_identity.nc")
nc_file <- nc_open("data/nafo_fishingfoot_2015_2099_cmip22_sorall_identity_mavg_only.nc")
nc_dep <- nc_open("data/nafo_fishingfoot_1993_2014_glorys.nc")

# Extract dimension information
lon <- ncvar_get(nc_file, "lon")
lat <- ncvar_get(nc_file, "lat")
time <- ncvar_get(nc_file, "time")
ssp <- ncvar_get(nc_file, "lev")  # ssp = climate scenarios
ens <- ncvar_get(nc_file, "ens")
vars <- names(nc_file$var)[!(names(nc_file$var) == "areavg")]
# model_names <- ncatt_get(nc_file, 0, "ens_proxy")$value  # 0 corresponds to global attribute rather than variable
# model_names <- str_split_1(model_names, pattern = " ")[-1] %>%
#   str_remove_all("[:digit:]=")

# Convert time to proper dates
time_units <- ncatt_get(nc_file, "time", "units")$value
time_origin <- as.POSIXct("1900-01-01 00:00:00", tz = "UTC")
date <- time_origin + time * 3600  # Convert hours to seconds
rm(time_origin, time_units)

# Function to select appropriate cmocean palette based on variable
get_cmocean_palette <- function(variable) {
  palette_map <- list(
    "tosavg" = "thermal",    # Sea surface temperature
    "tobavg" = "thermal",    # Bottom temperature
    "sosavg" = "haline",     # Sea surface salinity
    "sobavg" = "haline",     # Bottom salinity
    "mldavg" = "deep",       # Mixed layer depth
    "wobavg" = "speed",      # Bottom water velocity
    "bstress" = "amp"        # Bottom stress
  )
  
  return(palette_map[[variable]] %||% "thermal")  # Default to thermal if not found
}

# Get bathymetry data from NOAA for contours
# bathy_noaa <- marmap::getNOAA.bathy(lon1 = min(lon)-0.09, lon2 = max(lon)+0.09,
#                       lat1 = min(lat)-0.09, lat2 = max(lat)+0.09,
#                       resolution = 4) %>%
#   marmap::as.xyz()
# colnames(bathy_noaa) <- c("x","y","z")
# saveRDS(bathy_noaa, "data/bathy_noaa.rds")
bathy_noaa <- readRDS("data/bathy_noaa.rds")

# Get study area extent and create spatial mask ----
sa <- terra::rast("data/NRA_BNAM_b_tmp_mean.tif") %>%
  terra::as.polygons(.) %>%
  sf::st_as_sf(.) %>%
  sf::st_transform(4326) %>%
  mutate(NRA_BNAM_b_tmp_mean = 1) %>%
  group_by(NRA_BNAM_b_tmp_mean) %>%
  summarise(geometry = sf::st_union(geometry))

sa_lims <- sf::st_coordinates(sa) %>% as.data.frame %>%
  select(X,Y)
sa_lims <- c(min(sa_lims$X), max(sa_lims$X),min(sa_lims$Y),max(sa_lims$Y))

# Create spatial mask for study area
create_spatial_mask <- function(study_area_sf) {
  
  coords <- expand.grid(lon = as.vector(lon), lat = as.vector(lat))
  coords_info <- list(lon = lon, lat = lat, coords = coords)
  
  # Convert coordinates to sf points
  coords_sf <- sf::st_as_sf(coords_info$coords, 
                        coords = c("lon", "lat"), 
                        crs = sf::st_crs(study_area_sf))
  
  # Check which points are within study area
  within_study_area <- sf::st_within(coords_sf, study_area_sf, sparse = FALSE)
  mask <- rowSums(within_study_area) > 0
  
  # Convert back to matrix form matching netCDF dimensions
  mask_matrix <- matrix(mask, 
                        nrow = length(coords_info$lon), 
                        ncol = length(coords_info$lat))
  
  return(mask_matrix)
}
sa_mask <- create_spatial_mask(sa)
rm(create_spatial_mask)

# Close the file
nc_close(nc_file)


