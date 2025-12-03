
# Functions for extracting data

# Extract all data and apply spatial mask to all dimensions
extract_data_constrained <- function(variable) {
  data <- ncvar_get(nc_file, variable)
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

dep <- ncvar_get(nc_dep, "depavg")

# Convert to lat/lon dataframe
dep_df <- as.data.frame(as.table(dep))
colnames(dep_df) <- c("lon_idx", "lat_idx", "dep")
dep_df$lon <- lon[dep_df$lon_idx]
dep_df$lat <- lat[dep_df$lat_idx]
dep_df <- dep_df %>%
  select(lon, lat, dep)

# Extract all data from NC file
ens_data <- lapply(vars, extract_data_constrained)

# Convert each 4D array in the list to single long-form dataframe
ens_df <- lapply(1:length(vars), function(var_idx) {
  variable <- vars[var_idx]
  data_array <- ens_data[[var_idx]]
  
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
  reduce(left_join, by = c("ssp", "date", "lon", "lat")) %>%
  mutate(ssp = as.character(ssp),
         ssp = case_when(
           ssp == "126" ~ "1-2.6",
           ssp == "245" ~ "2-4.5",
           ssp == "370" ~ "3-7.0",
           ssp == "585" ~ "5-8.5"
         ),
         period = cut(date,
                      breaks = as.Date(c("2020-01-01", "2039-12-31", "2059-12-31", "2079-12-31", "2099-12-31")),
                      labels = c("P1", "P2", "P3", "P4"),
                      right = TRUE),
         year = as.integer(substr(date, 1, 4)),
         month = as.integer(substr(date, 6, 7)),
         season = case_when(month %in% c(1,2,3) ~ "W",
                            month %in% c(4,5,6) ~ "Sp",
                            month %in% c(7,8,9) ~ "Su",
                            month %in% c(10,11,12) ~ "F")) %>%
  left_join(dep_df, by = c("lon", "lat")) %>%
  mutate(pressure = gsw::gsw_p_from_z(-dep, lat),
         abs_sal = gsw::gsw_SA_from_SP(sobavg, pressure, lon, lat),
         con_temp = gsw::gsw_CT_from_t(abs_sal, tobavg, pressure),
         sw_dens = gsw::gsw_rho(abs_sal, con_temp, pressure),
         bstress = 3.5 * 10^-3 * sw_dens * wobavg^2) %>%
  select(lon, lat, date, ssp, period, year, month, season, dep, all_of(vars), bstress)



saveRDS(ens_df, "data/ens_df.rds")
rm(ens_data,dep,dep_df)



