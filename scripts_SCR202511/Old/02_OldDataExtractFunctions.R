
# Former functions for loading one variable in a list

# Function to extract data for a specific variable, time, ssp, and ensemble
extract_data <- function(variable, ssp_idx = 1, time_idx = 1, ens_idx = 1) {
  data <- ncvar_get(nc_file, variable, 
                    start = c(1, 1, ssp_idx, time_idx, ens_idx),
                    count = c(-1, -1, 1, 1, 1))  # -1 = all (for lat, lon)
  
  # Replace missing values with NA
  data[data == -9999] <- NA
  return(data)
}

# Modified extraction function with spatial constraint
extract_data_constrained <- function(variable, ssp_idx = 1, time_idx = 1, ens_idx = 1) {
  # Extract data as before
  data <- ncvar_get(nc_file, variable, 
                    start = c(1, 1, ssp_idx, time_idx, ens_idx),
                    count = c(-1, -1, 1, 1, 1))
  
  # Replace missing values with NA
  data[data == -9999] <- NA
  
  # Apply spatial mask
  data[!sa_mask] <- NA
  
  return(data)
}


# ens_tosavg <- extract_data_constrained("tosavg")
# 
# # Convert array to long-form dataframe
# ens_tosavg_df <- bind_rows(lapply(1:dim(ens_tosavg)[3], function(ssp_idx) {
#   bind_rows(lapply(1:dim(ens_tosavg)[4], function(time_idx) {
#     data_matrix <- ens_tosavg[,,ssp_idx,time_idx]
#     df <- as.data.frame(as.table(data_matrix))
#     colnames(df) <- c("lon_idx", "lat_idx", "tosavg")
#     df$lon <- lon[df$lon_idx]
#     df$lat <- lat[df$lat_idx]
#     df$date <- date[time_idx]
#     df$ssp <- ssp[ssp_idx]
#     return(df)
#   }))
# }), .id = "ssp") %>%
#   select(ssp, date, lon, lat, tosavg) %>%
#   filter(!is.na(tosavg))
# 


# Function to calculate ensemble statistics
calculate_ensemble_stats_constrained <- function(variable, ssp_idx = 1, time_idx = 1) {
  # Extract data for all ensemble members
  all_data <- array(NA, dim = c(length(lon), length(lat), length(ens)))
  
  for (i in 1:length(ens)) {
    data <- extract_data_constrained(variable, ssp_idx, time_idx, i)
    all_data[, , i] <- data
  }
  
  # Calculate statistics across ensemble dimension
  ensemble_mean <- apply(all_data, c(1, 2), mean, na.rm = FALSE)
  
  return(ensemble_mean)
}


ens_tosavg_ssp_date <- lapply(1:length(ssp), function(ssp_idx) {
  result <- lapply(1:length(time), function(time_idx) {
    calculate_ensemble_stats_constrained("tosavg", ssp_idx = ssp_idx, time_idx = time_idx)
  })
  names(result) <- date
  return(result)
})
names(ens_tosavg_ssp_date) <- ssp

# Convert from a list to a long-form dataframe
ens_tosavg_ssp_date_df <- bind_rows(lapply(names(ens_tosavg_ssp_date), function(ssp_name) {
  bind_rows(lapply(names(ens_tosavg_ssp_date[[ssp_name]]), function(date_name) {
    data_matrix <- ens_tosavg_ssp_date[[ssp_name]][[date_name]]
    df <- as.data.frame(as.table(data_matrix))
    colnames(df) <- c("lon_idx", "lat_idx", "tosavg")
    df$lon <- lon[df$lon_idx]
    df$lat <- lat[df$lat_idx]
    df$date <- as.Date(date_name)
    df$ssp <- ssp_name
    return(df)
  }))
}), .id = "ssp") %>%
  select(ssp, date, lon, lat, tosavg) %>%
  filter(!is.na(tosavg))

# saveRDS(ens_tosavg_ssp_date, file = "data/ens_tosavg_ssp_date.rds")


# Data summarising functions

# Define desired 20-year time periods, P1-P4
time_periods <- list(
  P1 = seq(as.POSIXct("2020-01-15", tz = "UTC"), as.POSIXct("2039-12-15", tz = "UTC"), by = "month"),
  P2 = seq(as.POSIXct("2040-01-15", tz = "UTC"), as.POSIXct("2059-12-15", tz = "UTC"), by = "month"),
  P3 = seq(as.POSIXct("2060-01-15", tz = "UTC"), as.POSIXct("2079-12-15", tz = "UTC"), by = "month"),
  P4 = seq(as.POSIXct("2080-01-15", tz = "UTC"), as.POSIXct("2099-12-15", tz = "UTC"), by = "month")
)

date_idx <- which(date %in% c(time_periods$P1, time_periods$P2, time_periods$P3, time_periods$P4))
time_periods_idx <- list(
  P1 = which(date %in% time_periods$P1),
  P2 = which(date %in% time_periods$P2),
  P3 = which(date %in% time_periods$P3),
  P4 = which(date %in% time_periods$P4)
)


# Average tosavg over each period for each ssp ----
ens_tosavg_ssp_period <- lapply(1:length(ssp), function(ssp_idx) {
  result <- lapply(1:length(time_periods), function(tp_idx) {
    period_data <- ens_tosavg_ssp_date[[ssp_idx]][time_periods_idx[[tp_idx]]]
    period_mean <- Reduce("+", period_data) / length(period_data)
    return(period_mean)
  })
  names(result) <- names(time_periods)
  return(result)
})
names(ens_tosavg_ssp_period) <- ssp

# Summarise over the spatial extent (mean, min, max, sd) to have one value per time period and ssp
ens_tosavg_ssp_period_summary <- lapply(1:length(ssp), function(ssp_idx) {
  ssp_name <- ssp[ssp_idx]
  period_summaries <- lapply(1:length(time_periods), function(tp_idx) {
    period_name <- names(time_periods)[tp_idx]
    data <- ens_tosavg_ssp_period[[ssp_idx]][[tp_idx]]
    
    summary <- c(
      mean = mean(data, na.rm = TRUE),
      min = min(data, na.rm = TRUE),
      max = max(data, na.rm = TRUE),
      sd = sd(data, na.rm = TRUE)
      # range = diff(range(data, na.rm = TRUE))
    )
    return(data.frame(SSP = ssp_name, period = period_name, t(summary)))
  })
  do.call(rbind, period_summaries)
})
names(ens_tosavg_ssp_period_summary) <- ssp
ens_tosavg_ssp_period_summary_df <- do.call(rbind, ens_tosavg_ssp_period_summary)
ens_tosavg_ssp_period_summary_df$period <- as.factor(ens_tosavg_ssp_period_summary_df$period)
ens_tosavg_ssp_period_summary_df$SSP <- as.factor(ens_tosavg_ssp_period_summary_df$SSP)




