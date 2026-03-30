
# Load historical CMIP data from 1993-2014
library(tidyverse)
library(ncdf4)
cmip_path <- R.utils::readWindowsShortcut("data/raw/NAFO_CC_Projections_Data.lnk")$pathname
nc_file <- nc_open(file.path(cmip_path, "nafo_fishingfoot_1993_2014_cmip22_sorall_identity.nc"))

# Read in glorys dataset for depth (will use to calculate bstress)
nc_dep <- nc_open(file.path(cmip_path, "nafo_fishingfoot_1993_2014_glorys.nc"))
dep <- ncvar_get(nc_dep, "depavg")
dep_df <- as.data.frame.table(dep)
colnames(dep_df) <- c("lon_idx", "lat_idx", "dep")
dep_df$lon <- round(ncvar_get(nc_dep, "lon")[dep_df$lon_idx], 5)
dep_df$lat <- round(ncvar_get(nc_dep, "lat")[dep_df$lat_idx], 5)
dep_df <- select(dep_df, lon, lat, dep)

vars <- names(nc_file$var)[!(names(nc_file$var) == "areavg")]

nc_data <- lapply(vars, function(var) {
  data <- ncvar_get(nc_file, var)
  data[data == -9999] <- NA

  # Summarise over CMIP models
  data <- rowMeans(data, dims = 3, na.rm = TRUE)

  # Convert to dataframe
  dimnames(data) <- list(
    lon = round(ncvar_get(nc_file, "lon"), 5),
    lat  = round(ncvar_get(nc_file, "lat"), 5),
    time = ncvar_get(nc_file, "time")
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
  select(-c(pressure, abs_sal, con_temp, sw_dens, dep))

write_csv(nc_data, "data/cleaned/cmip_ens_1993_2014_df.csv")

