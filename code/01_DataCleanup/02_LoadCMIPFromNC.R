library(tidyverse)
cmip_path <- R.utils::readWindowsShortcut("data/raw/NAFO_CC_Projections_Data.lnk")$pathname

# Read in glorys dataset for depth (will use to calculate bstress) ----
nc_dep <- ncdf4::nc_open(file.path(cmip_path, "nafo_fishingfoot_1993_2014_glorys.nc"))
dep <- ncdf4::ncvar_get(nc_dep, "depavg")
dep_df <- as.data.frame.table(dep)
colnames(dep_df) <- c("lon_idx", "lat_idx", "dep")
dep_df$lon <- round(ncdf4::ncvar_get(nc_dep, "lon")[dep_df$lon_idx], 5)
dep_df$lat <- round(ncdf4::ncvar_get(nc_dep, "lat")[dep_df$lat_idx], 5)
dep_df <- select(dep_df, lon, lat, dep)

# Load baseline/reference CMIP data from 1993-2014 and transform into usable format ----
nc_file <- ncdf4::nc_open(file.path(cmip_path, "nafo_fishingfoot_1993_2014_cmip22_sorall_identity.nc"))
vars <- names(nc_file$var)[!(names(nc_file$var) == "areavg")]

nc_data <- lapply(vars, function(var) {
  data <- ncdf4::ncvar_get(nc_file, var)
  data[data == -9999] <- NA

  # Summarise over CMIP models
  data <- rowMeans(data, dims = 3, na.rm = TRUE)

  # Convert to dataframe
  dimnames(data) <- list(
    lon = round(ncdf4::ncvar_get(nc_file, "lon"), 5),
    lat  = round(ncdf4::ncvar_get(nc_file, "lat"), 5),
    time = ncdf4::ncvar_get(nc_file, "time")
  )
  data <- as.data.frame.table(data, responseName = var) %>%
    mutate(
      across(lon:time, ~ as.numeric(as.character(.x))),
      # Time conversion to date
      time = as.POSIXct("1900-01-01 00:00:00", tz = "UTC") + time * 3600
  ) %>%
    rename(date = time)
}) %>% 
  reduce(left_join, by = c("lon", "lat", "date")) %>%
  set_names("lon", "lat", "date", vars) %>%
  mutate(
    month = as.integer(substr(date, 6, 7)),
    season = case_when(
      month %in% c(1,2,3) ~ "W",
      month %in% c(4,5,6) ~ "Sp",
      month %in% c(7,8,9) ~ "Su",
      month %in% c(10,11,12) ~ "F"
    )
  ) %>%
  pivot_wider(names_from = season, values_from = mldavg, names_glue = {"{.value}_{season}"}) %>%
  left_join(dep_df, by = c("lon", "lat")) %>%
  mutate(
    mldavg = coalesce(mldavg_W,mldavg_F,mldavg_Su,mldavg_Sp),
    pressure = gsw::gsw_p_from_z(-dep, lat),
    abs_sal = gsw::gsw_SA_from_SP(sobavg, pressure, lon, lat),
    con_temp = gsw::gsw_CT_from_t(abs_sal, tobavg, pressure),
    sw_dens = gsw::gsw_rho(abs_sal, con_temp, pressure),
    bstress = 3.5 * 10^-3 * sw_dens * wobavg^2
  ) %>%
  select(-c(pressure, abs_sal, con_temp, sw_dens, dep, date, month))

saveRDS(nc_data, "data/processed/cmip_ens_1993_2014_df.rds")

# Load projected CMIP data ----
nc_file <- ncdf4::nc_open(file.path(cmip_path, "nafo_fishingfoot_2015_2099_cmip22_sorall_identity_mavg_only.nc"))  # can't use this for uncertainty
vars <- names(nc_file$var)[!(names(nc_file$var) == "areavg")]

# Apply spatial mask to only keep cells within study area ----
sa <- terra::rast("data/raw/Bathy_Layers/GEBCO2024_FS005.tif") %>%
  terra::as.polygons(.) %>%
  sf::st_as_sf(.) %>%
  sf::st_transform(4326) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  sf::st_make_valid()

sa_lims <- sf::st_coordinates(sa) %>% as.data.frame %>%
  select(X,Y)
sa_lims <- c(min(sa_lims$X), max(sa_lims$X),min(sa_lims$Y),max(sa_lims$Y))

create_spatial_mask <- function() {
  
  lon <- ncdf4::ncvar_get(nc_file, "lon")
  lat <- ncdf4::ncvar_get(nc_file, "lat")
  coords <- expand.grid(lon = as.vector(lon), lat = as.vector(lat))
  coords_info <- list(lon = lon, lat = lat, coords = coords)
  
  # Convert coordinates to sf points
  coords_sf <- sf::st_as_sf(coords_info$coords, coords = c("lon", "lat"), crs = sf::st_crs(sa))
  
  # Check which points are within study area
  within_study_area <- sf::st_within(coords_sf, sa, sparse = FALSE)
  mask <- rowSums(within_study_area) > 0
  
  # Convert back to matrix form matching netCDF dimensions
  mask_matrix <- matrix(
    mask, 
    nrow = length(coords_info$lon), 
    ncol = length(coords_info$lat))
  
  return(mask_matrix)
}
sa_mask <- create_spatial_mask()
rm(create_spatial_mask)

extract_data_constrained <- function(variable) {
  data <- ncdf4::ncvar_get(nc_file, variable)
  data[data == -9999] <- NA

  # Apply spatial mask across all dimensions
  for (ssp_idx in 1:dim(data)[3]) {
    for (time_idx in 1:dim(data)[4]) {
      slice <- data[,,ssp_idx,time_idx]
      slice[!sa_mask] <- NA
      data[,,ssp_idx,time_idx] <- slice
    }
  }
  return(data)
}

ens_data <- lapply(vars, extract_data_constrained)

# Get data into usable format ----
ens_df <- lapply(1:length(vars), function(var_idx) {
  variable <- vars[var_idx]
  data_array <- ens_data[[var_idx]]
  lon <- round(ncdf4::ncvar_get(nc_file, "lon"), 5)
  lat <- round(ncdf4::ncvar_get(nc_file, "lat"), 5)
  time <- ncdf4::ncvar_get(nc_file, "time")
  time_units <- ncdf4::ncatt_get(nc_file, "time", "units")$value
  time_origin <- as.POSIXct("1900-01-01 00:00:00", tz = "UTC")
  date <- time_origin + time * 3600  # Convert hours to seconds
  ssp <- ncdf4::ncvar_get(nc_file, "lev")

  bind_rows(lapply(1:dim(data_array)[3], function(ssp_idx) {
    bind_rows(lapply(1:dim(data_array)[4], function(time_idx) {
      data_matrix <- data_array[,,ssp_idx,time_idx]
      df <- as.data.frame(as.table(data_matrix))
      colnames(df) <- c("lon_idx", "lat_idx", variable)
      df$lon <- lon[df$lon_idx]
      df$lat <- lat[df$lat_idx]
      df$date <- date[time_idx]
      df$ssp <- ssp[ssp_idx]
      return(df)
    }))
  })) %>%
    select(ssp, date, lon, lat, all_of(variable)) %>%
    filter(!is.na(get(variable))) %>%
    mutate(date = as.Date(date))
}) %>% 
  reduce(left_join, by = c("lon", "lat", "ssp", "date")) %>%
  mutate(
    ssp = as.character(ssp),
    ssp = case_when(
      ssp == "126" ~ "1-2.6",
      ssp == "245" ~ "2-4.5",
      ssp == "370" ~ "3-7.0",
      ssp == "585" ~ "5-8.5"),
    period = cut(
      date,
      breaks = as.Date(c("2020-01-01", "2039-12-31", "2059-12-31", "2079-12-31", "2099-12-31")),
      labels = c("P1", "P2", "P3", "P4"),
      right = TRUE),
    year = as.integer(substr(date, 1, 4)),
    month = as.integer(substr(date, 6, 7)),
    season = case_when(
      month %in% c(1,2,3) ~ "W",
      month %in% c(4,5,6) ~ "Sp",
      month %in% c(7,8,9) ~ "Su",
      month %in% c(10,11,12) ~ "F"
    )
  ) %>%
  pivot_wider(names_from = season, values_from = mldavg, names_glue = {"{.value}_{season}"}) %>%
  left_join(dep_df, by = c("lon", "lat")) %>%
  mutate(
    mldavg = coalesce(mldavg_W,mldavg_F,mldavg_Su,mldavg_Sp),
    pressure = gsw::gsw_p_from_z(-dep, lat),
    abs_sal = gsw::gsw_SA_from_SP(sobavg, pressure, lon, lat),
    con_temp = gsw::gsw_CT_from_t(abs_sal, tobavg, pressure),
    sw_dens = gsw::gsw_rho(abs_sal, con_temp, pressure),
    bstress = 3.5 * 10^-3 * sw_dens * wobavg^2
  ) %>%
  select(-c(pressure, abs_sal, con_temp, sw_dens, dep, date, year, month))

saveRDS(ens_df, "data/processed/cmip_ens_proj_df.rds")
