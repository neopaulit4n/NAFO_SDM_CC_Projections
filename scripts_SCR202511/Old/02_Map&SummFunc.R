
# Calculate each of the variables (min, max, range, mean) for Period 1 
# (P1: 2020-2039), Period 2 (P2: 2040-2059), Period 3 (P3: 2060-2079), and 
# Period 4 (P4: 2080-2099) for each of the 4 SSPs. As you do each one can you 
# provide me with a figure for the report that is of a good quality and when 
# each variable is done create the graphics as is shown in this doc for the SST 
# but for just the spatial extent for the SDM modeling? It would also be good to 
# get a good figure of that spatial extent.


# Functions for extracting data and mapping

# Function to extract data for a specific variable, time, ssp, and ensemble
extract_data <- function(variable, ssp_idx = 1, time_idx = 1, ens_idx = 1) {
  data <- ncvar_get(nc_file, variable, 
                    start = c(1, 1, ssp_idx, time_idx, ens_idx),
                    count = c(-1, -1, 1, 1, 1))  # -1 = all (for lat, lon)
  
  # Replace missing values with NA
  data[data == -9999] <- NA
  return(data)
}

nc_list_tosavg <- lapply(1:length(ssp), function(ssp_idx) {
  lapply(date_idx, function(time_idx) {
    lapply(1:length(ens), function(ens_idx) {
      extract_data("tosavg", ssp_idx, time_idx, ens_idx)
    })
  })
})




nc_array <- ncvar_get(nc_file, "tosavg")
nc_array_ssplist <- lapply(1:length(ssp), 
                           function(x) ncvar_get(nc_file, "tosavg", start = c(1,1,x,1,1), count = c(-1,-1,1,-1,-1)))
# Convert each array in the list to a long-form dataframe
nc_df_list <- lapply(nc_array_ssplist, function(arr) {
  df <- as.data.frame.table(arr, responseName = "tosavg")
  colnames(df) <- c("lon", "lat", "date", "ens", "tosavg")
  df$lon <- lon
  df$lat <- lat
  df$date <- time
  df$ens <- ens
  return(df)
})


df <- expand_grid(
  lon = lon,
  lat = lat,
  date = time,
  ens = ens
)




calc_period_stats <- function(var_data, time_indices) {
  period_data <- var_data[,,,time_indices,]  # lon, lat, ssp, time, ens
  
  list(
    min = apply(period_data, c(1,2,3,5), min, na.rm = TRUE),
    max = apply(period_data, c(1,2,3,5), max, na.rm = TRUE),
    mean = apply(period_data, c(1,2,3,5), mean, na.rm = TRUE),
    range = apply(period_data, c(1,2,3,5), range, na.rm = TRUE)
  )
}

# Apply to each variable and period
results <- list()
for(var_name in vars) {
  var_data <- ncvar_get(nc, var_name)
  results[[var_name]] <- lapply(periods, function(p) calc_period_stats(var_data, p))
}

df_tosavg_p1_stats <- calc_period_stats(df, dates_p1_idx)
period_data <- df[,,,]

# extract_data <- function(variable) {
#   data <- ncvar_get(nc_file, variable, 
#                     start = c(1, 1, 1, 1, 1),
#                     count = c(-1, -1, -1, -1, -1))  # -1 = all (for lat, lon)
#   
#   # Replace missing values with NA
#   data[data == -9999] <- NA
#   return(data)
# }
# 
# df_tosavg <- extract_data("tosavg")
# 
# # Mean across models
# df_tosavg_ensmean <- apply(df_tosavg, c(1,2,3,4), mean, na.rm = TRUE)
# 
# 
# 
# # Crop to study area
# # spatial_mask <- create_spatial_mask(sa)
# # mask_4d <- array(spatial_mask, dim = dim(df_mean))
# # df_mean_crop <- df_mean[!mask_4d] <- NA
# 
# df_tosavg11 <- df_tosavg[,,1,1,1]
# # tranform into 2-column dataframe
# df_tosavg11 <- as.data.frame(as.table(df_tosavg11))
# colnames(df_tosavg11) <- c("lon", "lat", "tosavg_mean")
# df_tosavg11$lon <- lon
# df_tosavg11$lat <- lat
# df_tosavg11 <- st_as_sf(df_tosavg11, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
# 
# plot(df_tosavg11$tosavg_mean)
# ggplot(data = df_tosavg11) +
#   geom_sf(aes(color = tosavg_mean)) +
#   scale_color_viridis_c()

# Create spatial mask for study area
create_spatial_mask <- function(study_area_sf) {
  
  coords <- expand.grid(lon = as.vector(lon), lat = as.vector(lat))
  coords_info <- list(lon = lon, lat = lat, coords = coords)
  
  # Convert coordinates to sf points
  coords_sf <- st_as_sf(coords_info$coords, 
                        coords = c("lon", "lat"), 
                        crs = st_crs(study_area_sf))
  
  # Check which points are within study area
  within_study_area <- st_within(coords_sf, study_area_sf, sparse = FALSE)
  mask <- rowSums(within_study_area) > 0
  
  # Convert back to matrix form matching netCDF dimensions
  mask_matrix <- matrix(mask, 
                        nrow = length(coords_info$lon), 
                        ncol = length(coords_info$lat))
  
  return(mask_matrix)
}

# Modified extraction function with spatial constraint
extract_data_constrained <- function(variable, study_area_sf, 
                                     ssp_idx = 1, time_idx = 1, ens_idx = 1) {
  # Extract data as before
  data <- ncvar_get(nc_file, variable, 
                    start = c(1, 1, ssp_idx, time_idx, ens_idx),
                    count = c(-1, -1, 1, 1, 1))
  
  # Replace missing values with NA
  data[data == -9999] <- NA
  
  # Create and apply spatial mask
  spatial_mask <- create_spatial_mask(study_area_sf)
  data[!spatial_mask] <- NA
  
  return(data)
}


# Function to create a spatial map for a given variable
plot_spatial_map <- function(variable, time_idx = 1, ssp_idx = 1, ens_idx = 1, 
                             add_bathymetry = TRUE, 
                             bathy_contours = c(-50, -100, -200, -500, -1000, -2000, -3000, -4000, -5000)) {
  data <- extract_data(variable, time_idx, ssp_idx, ens_idx)
  
  # Create a data frame for plotting
  df <- expand.grid(lon = lon, lat = lat)
  df$value <- as.vector(data)
  
  # Remove NA values for plotting
  df <- df[!is.na(df$value), ]
  
  # Get model name
  model_name <- model_names[ens_idx]
  
  # Create the plot
  palette_name <- get_cmocean_palette(variable)
  
  p <- ggplot(df, aes(x = lon, y = lat, fill = value)) +
    geom_tile() +
    scale_fill_cmocean(variable, name = palette_name)
  
  # Add bathymetry contours
  if (add_bathymetry) {
    p <- p + geom_contour(data = bathy_noaa, 
                          aes(x = x, y = y, z = z, fill = NULL), 
                          breaks = bathy_contours,
                          color = "white", 
                          linewidth = 0.3, 
                          alpha = 0.4)
  }
  
  p <- p + geom_sf(data = sa, colour = "black", alpha = 0, inherit.aes = FALSE) +
    scale_x_continuous(limits = c(min(lon)-0.09,max(lon)+0.09), expand = c(0, 0)) +
    scale_y_continuous(limits = c(min(lat)-0.09,max(lat)+0.09), expand = c(0, 0)) +
    # facet_wrap(~ ens_idx) +
    theme_minimal() +
    labs(title = paste(variable, "- Model:", model_name),
         subtitle = paste("Date:", format(dates[time_idx], "%Y-%m"), 
                          "SSP:", ssp[ssp_idx]),
         x = "Longitude", y = "Latitude") +
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5),
          axis.ticks = element_line(colour = "black"))
  
  return(p)
}

# Function to calculate ensemble statistics
calculate_ensemble_stats <- function(variable, ssp_idx = 1, time_idx = 1) {
  # Extract data for all ensemble members
  all_data <- array(NA, dim = c(length(lon), length(lat), length(ens)))
  
  for (i in 1:length(ens)) {
    data <- extract_data(variable, ssp_idx, time_idx, i)
    all_data[, , i] <- data
  }
  
  # Calculate statistics across ensemble dimension
  ensemble_mean <- apply(all_data, c(1, 2), mean, na.rm = FALSE)
  ensemble_sd <- apply(all_data, c(1, 2), sd, na.rm = FALSE)
  ensemble_min <- apply(all_data, c(1, 2), min, na.rm = FALSE)
  ensemble_max <- apply(all_data, c(1, 2), max, na.rm = FALSE)
  
  # return(ensemble_mean)
  return(list(mean = ensemble_mean, sd = ensemble_sd,
              min = ensemble_min, max = ensemble_max))
}


ens_all <- lapply(1:5, function (x) calculate_ensemble_stats("tosavg", time_idx = x, ssp_idx = 1))


# Calculate ensemble stats for every cell/time/ssp
ens_p1_tosavg <- lapply(1:length(ssp),
                     function(y) lapply(dates_p1_idx,  #1:length(time),
                                        function(x) calculate_ensemble_stats("tosavg", time_idx = x, ssp_idx = y)))

# ens_ssp585 <- lapply(1:length(time), 
#        function(x) calculate_ensemble_stats("tosavg", time_idx = x, ssp_idx = 4))
# 
annual_means <- lapply(1:85, function(year) {
  start_idx <- (year - 1) * 12 + 1
  end_idx <- year * 12

  # Extract 12 months for this year
  year_data <- lapply(ens_ssp585[start_idx:end_idx], function(month) month[[1]])

  # Convert to array and calculate mean
  year_array <- array(unlist(year_data),
                      dim = c(nrow(year_data[[1]]),
                              ncol(year_data[[1]]),
                              12))

  apply(year_array, c(1,2), mean, na.rm = TRUE)
})

domain_means <- lapply(1:4, function(ssp) {
  sapply(1:1020, function(month) {
    mean(ens_all[[ssp]][[month]][[1]])
  })
})
names(domain_means) <- paste0("SSP", ssp)
domain_means <- bind_cols(domain_means) %>%
  mutate(date = dates) %>%
  pivot_longer(cols = -date, names_to = "SSP", values_to = "SST")

ggplot(data = domain_means, aes(x = date, y = SST, colour = SSP)) +
  theme_bw() +
  geom_line() +
  stat_smooth(method = "lm")

annual_means <- lapply(1:4, function(ssp) {
  lapply(1:85, function(year) {
    start_idx <- (year - 1) * 12 + 1
    end_idx <- year * 12
    
    # Extract 12 months for this year
    year_data <- lapply(ens_all[[ssp]][start_idx:end_idx], function(month) month[[1]])
    
    # Convert to array and calculate mean
    year_array <- array(unlist(year_data), 
                        dim = c(nrow(year_data[[1]]), 
                                ncol(year_data[[1]]), 
                                12))
    
    apply(year_array, c(1,2), mean, na.rm = TRUE)
  })
})

annual_means_domain <- lapply(1:4, function(ssp) {
  sapply(1:85, function(year) mean(annual_means[[ssp]][[year]]))
})
names(annual_means_domain) <- paste0("SSP",ssp)
annual_means_domain <- bind_cols(annual_means_domain) %>%
  mutate(year = unique(year(dates))) %>%
  pivot_longer(cols = starts_with("SSP"), names_to = "SSP", values_to = "SST")

ggplot(data = annual_means_domain, aes(x = year, y = SST, colour = SSP)) +
  theme_bw() +
  geom_line() +
  stat_smooth(method = "lm", linetype = 2, linewidth = 0.5, se = FALSE) +
  xlab("Year")


# Function to plot ensemble statistics
plot_ensemble_stats <- function(variable, time_idx = 1, ssp_idx = 1, stat = "mean",
                                add_bathymetry = TRUE, 
                                bathy_contours = c(-50, -100, -200, -500, -1000, -2000, -3000, -4000, -5000)) {
  stats <- calculate_ensemble_stats(variable, time_idx, ssp_idx)
  
  # Create data frame
  df <- expand.grid(lon = lon, lat = lat)
  df$value <- as.vector(stats[[stat]])
  df <- df[!is.na(df$value), ]
  
  palette_name <- get_cmocean_palette(variable)
  
  p <- ggplot(df, aes(x = lon, y = lat, fill = value)) +
    geom_tile() +
    scale_fill_cmocean(variable, name = palette_name)
  
  # Add bathymetry contours if requested
  if (add_bathymetry) {
    p <- p + geom_contour(data = bathy_noaa, 
                          aes(x = x, y = y, z = z, fill = NULL), 
                          breaks = bathy_contours,
                          color = "white", 
                          linewidth = 0.3, 
                          alpha = 0.4)
  }
  
  p <- p + geom_sf(data = sa, colour = "black", alpha = 0, inherit.aes = FALSE) +
    scale_x_continuous(limits = c(min(lon)-0.09,max(lon)+0.09), expand = c(0, 0)) +
    scale_y_continuous(limits = c(min(lat)-0.09,max(lat)+0.09), expand = c(0, 0)) +
    theme_minimal() +
    labs(title = paste("Ensemble", stringr::str_to_title(stat), "-", variable),
         subtitle = paste("Date:", format(dates[time_idx], "%Y-%m"), 
                          "SSP:", ssp[ssp_idx]),
         x = "Longitude", y = "Latitude") +
    theme(axis.ticks = element_line(colour = "black"))
  
  return(p)
}

# Function to create animated ensemble statistics over time
create_ensemble_animation <- function(variable, ssp_idx = 1, stat = "mean",
                                      time_subset = NULL, 
                                      add_bathymetry = TRUE,
                                      bathy_contours = c(-50, -100, -200, -500, -1000, -2000, -3000, -4000, -5000),
                                      fps = 2, width = 1000, height = 800,
                                      output_file = NULL) {
  
  # Determine time indices to use
  if (is.null(time_subset)) {
    # Use every 12th time step (roughly annual if monthly data)
    time_indices <- seq(1, length(dates), by = 12)
  } else {
    time_indices <- time_subset
  }
  
  message(paste("Creating animation with", length(time_indices), "time steps..."))
  
  # Calculate ensemble stats for all time steps
  all_stats <- list()
  for (i in seq_along(time_indices)) {
    t_idx <- time_indices[i]
    stats <- calculate_ensemble_stats(variable, t_idx, ssp_idx)
    
    # Create data frame for this time step
    df <- expand.grid(lon = lon, lat = lat)
    df$value <- as.vector(stats[[stat]])
    df$time_idx <- t_idx
    df$date <- dates[t_idx]
    df$year_month <- format(dates[t_idx], "%Y-%m")
    
    # Remove NA values
    df <- df[!is.na(df$value), ]
    
    all_stats[[i]] <- df
  }
  
  # Combine all time steps
  animation_data <- do.call(rbind, all_stats)
  
  # Get consistent color scale across all time steps
  value_range <- range(animation_data$value, na.rm = TRUE)
  palette_name <- get_cmocean_palette(variable)
  
  # Create base plot
  p <- ggplot(animation_data, aes(x = lon, y = lat, fill = value)) +
    geom_tile() +
    scale_fill_cmocean(variable, name = palette_name,
                       limits = value_range)
  
  # Add bathymetry contours if requested
  if (add_bathymetry) {
    p <- p + geom_contour(data = bathy_noaa, 
                          aes(x = x, y = y, z = z, fill = NULL), 
                          breaks = bathy_contours,
                          color = "white", 
                          linewidth = 0.3, 
                          alpha = 0.4)
  }
  
  # Complete the plot with animation
  p <- p + geom_sf(data = sa, colour = "black", alpha = 0, inherit.aes = FALSE) +
    scale_x_continuous(limits = c(min(lon)-0.09,max(lon)+0.09), expand = c(0, 0)) +
    scale_y_continuous(limits = c(min(lat)-0.09,max(lat)+0.09), expand = c(0, 0)) +
    theme_minimal() +
    labs(title = paste("Ensemble", stringr::str_to_title(stat), "-", variable),
         subtitle = "Date: {closest_state}",
         x = "Longitude", 
         y = "Latitude",
         caption = paste("SSP:", ssp[ssp_idx])) +
    theme(plot.title = element_text(size = 14, hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5),
          plot.caption = element_text(size = 10),
          axis.ticks = element_line(colour = "black")) +
    transition_states(year_month,
                      transition_length = 1,
                      state_length = 2) +
    ease_aes('linear')
  
  # Set output filename if not provided
  if (is.null(output_file)) {
    output_file <- paste0("output/ensemble_", stat, "_", variable, "_ssp", ssp_idx, "_animation.gif")
  }
  
  # Create and save animation
  message(paste("Rendering animation to:", output_file))
  anim <- animate(p, 
                  width = width, 
                  height = height, 
                  fps = fps,
                  duration = length(time_indices) / fps,
                  renderer = gifski_renderer(output_file))
  
  message("Animation complete!")
  return(anim)
}


# Make animated annual average SST plot ----

# Function to create animated ensemble statistics over time
create_ensemble_animation_annual_mean_tosavg <- function(variable, ssp_idx = 1, stat = "mean",
                                                         time_subset = NULL,
                                                         add_bathymetry = TRUE,
                                                         bathy_contours = c(-50, -100, -200, -500, -1000, -2000, -3000, -4000, -5000),
                                                         fps = 2, width = 1000, height = 800,
                                                         output_file = NULL) {
  
  message(paste("Creating animation..."))
  # Determine time indices to use
  if (is.null(time_subset)) {
    # Use every 12th time step (roughly annual if monthly data)
    time_indices <- seq(1, length(dates), by = 12)
  } else {
    time_indices <- time_subset
  }
  
  # Calculate ensemble stats for all time steps
  animated_list <- lapply(1:length(annual_means), function(year) {
    df <- expand.grid(lon = lon, lat = lat)
    df$value <- as.vector(annual_means[[year]])
    df$year <- as.vector(unique(year(dates))[year])
    return(df)
  })
  
  # Combine all time steps
  animation_data <- do.call(rbind, animated_list)
  
  # Get consistent color scale across all time steps
  value_range <- range(animation_data$value, na.rm = TRUE)
  palette_name <- get_cmocean_palette(variable)
  
  # Create base plot
  p <- ggplot(animation_data, aes(x = lon, y = lat, fill = value)) +
    geom_tile() +
    scale_fill_cmocean(variable, name = palette_name,
                       limits = value_range)
  
  # Add bathymetry contours if requested
  if (add_bathymetry) {
    p <- p + geom_contour(data = bathy_noaa, 
                          aes(x = x, y = y, z = z, fill = NULL), 
                          breaks = bathy_contours,
                          color = "white", 
                          linewidth = 0.3, 
                          alpha = 0.4)
  }
  
  # Complete the plot with animation
  p <- p + geom_sf(data = sa, colour = "black", alpha = 0, inherit.aes = FALSE) +
    scale_x_continuous(limits = c(min(lon)-0.09,max(lon)+0.09), expand = c(0, 0)) +
    scale_y_continuous(limits = c(min(lat)-0.09,max(lat)+0.09), expand = c(0, 0)) +
    theme_minimal() +
    labs(title = paste("Ensemble", stringr::str_to_title(stat), "-", variable),
         subtitle = "Date: {closest_state}",
         x = "Longitude", 
         y = "Latitude",
         caption = paste("SSP:", ssp[ssp_idx])) +
    theme(plot.title = element_text(size = 14, hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5),
          plot.caption = element_text(size = 10),
          axis.ticks = element_line(colour = "black")) +
    transition_states(year,
                      transition_length = 1,
                      state_length = 2) +
    ease_aes('linear')
  
  # Set output filename if not provided
  if (is.null(output_file)) {
    output_file <- paste0("output/ensemble_annual_", stat, "_", variable, "_ssp", ssp[ssp_idx], "_animation.gif")
  }
  
  # Create and save animation
  message(paste("Rendering animation to:", output_file))
  anim <- animate(p, 
                  width = width, 
                  height = height, 
                  fps = fps,
                  duration = length(time_indices) / fps,
                  renderer = gifski_renderer(output_file))
  
  message("Animation complete!")
  return(anim)
}





# Example usage:
# Plot temperature for first model, first time step, first level
plot1 <- plot_spatial_map("tosavg", time_idx = 1, ssp_idx = 1, ens_idx = 1)
print(plot1)

# Create grid plot for each individual model
plot_11x <- lapply(ens, function(x) plot_spatial_map("tosavg", time_idx = 1, ssp_idx = 1, ens_idx = x))
cowplot::plot_grid(plotlist = plot_11x)

# Plot ensemble mean
plot2 <- plot_ensemble_stats("tosavg", time_idx = 1, ssp_idx = 1, stat = "mean")
print(plot2)

# Plot time series for a specific location
# Example: middle of your domain
lon_middle <- round(length(lon)/2)
lat_middle <- round(length(lat)/2)
plot3 <- plot_timeseries("tosavg", lon_middle, lat_middle, ssp_idx = 1, 
                         ensemble_indices = 1:5)  # First 5 models for clarity
print(plot3)

# Animation
create_ensemble_animation("tosavg", ssp_idx = 1, stat = "mean",
                          time_subset = seq(1, length(dates), by = 12))  # by= number of months
create_ensemble_animation_annual_mean_tosavg("tosavg", ssp_idx = 4, stat = "mean",
                                             time_subset = seq(1, length(dates), by = 12)) 

# Print summary information
cat("Dataset Summary:\n")
cat("Variables:", vars, "\n")
cat("Spatial extent: Lon", range(lon), "Lat", range(lat), "\n")
cat("Time range:", format(range(dates), "%Y-%m-%d"), "\n")
cat("SSPs:", ssp, "\n")
cat("Number of ensemble members:", length(ens), "\n")

# Don't forget to close the file
nc_close(nc_file)

