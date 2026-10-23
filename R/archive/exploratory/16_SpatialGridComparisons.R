# Compare new CMIP layers to old layers

# Load layers ----
old_layers <- list.files(
  path = "output/SetResults/01_env_vars_rasters", 
  pattern = "BT_min",
  full.names = TRUE
) |>
  terra::rast()

new_layers <- list.files(
  path = "output/env_vars_rasters", 
  pattern = "BT_min",
  full.names = TRUE
) |>
  terra::rast()

# Calculate difference between layers ----
diff_layers <- new_layers - old_layers

# Visualise ----
terra::plot(old_layers[[1]])
terra::plot(new_layers[[1]])
terra::plot(diff_layers[[1]])