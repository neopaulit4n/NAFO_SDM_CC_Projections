
# Extrapolations with dsmextra

cat("Computing extrapolations...")
# samples = presence data points with associated covariate values (known) - points
# prediction.grid = grid of covariate values across the study area for which presence is to be predicted (unknown) - points

# vmeoi <- "black_corals"

# Convert selected vme layers to points
extrap_grid <- terra::as.data.frame(vme_layers_rast, xy = TRUE) %>%
  drop_na()

vme_pts <- cmip_comb_df %>%
  filter(VME_Group == vmeoi) %>%
  select(x = Start_Long_DD, y = Start_Lat_DD, all_of(selected_vme_vars)) %>%
  as.data.frame()

extrapolation_area <- dsmextra::compute_extrapolation(samples = vme_pts,
  covariate.names = selected_vme_vars,
  prediction.grid = extrap_grid,
  coordinate.system = sp::CRS(SRS_string = "EPSG:4326"))

# Determine which extrapolation types are not null
extrap_types <- c("univariate", "combinatorial", "analogue")[which(sapply(extrapolation_area$data, nrow)[2:4] > 0)]

# Extract extrapolation rasters
extrapolation_rasters <- lapply(c("ExDet", "mic"), function(method) {
  lapply(extrap_types, function(type) {
    terra::rast(extrapolation_area$rasters[[method]][[type]])
  }) %>%
    set_names(extrap_types)
}) %>%
  set_names(c("ExDet", "mic")) %>%
  unlist(recursive = FALSE)

extrapolation_maps <- lapply(1:length(extrapolation_rasters), function(x) {
  ggplot() +
    tidyterra::geom_spatraster(data = extrapolation_rasters[[x]]) +
    scale_fill_viridis_c(na.value = "transparent") +
    labs(title = names(extrapolation_rasters)[[x]]) +
    theme_minimal()  
})

extrapolation_map_grid <- cowplot::plot_grid(plotlist = extrapolation_maps, ncol = length(extrap_types))
ggsave(paste0("output/03_RF_Map_Outputs/",vmeoi,"_extrapolations_",poi,"_",sspoi,".jpg"), 
  plot = extrapolation_map_grid,
  width = 10, height = 5, dpi = 300)


# Load previous SDM2024 extrapolation maps for comparison ----

# Read TIF
# sdm2024_extana <- terra::rast("output/SDM2024_orig/Black Coral/BlackCorals_raster_output_sens_spe/BlackCorals_ext.analogue.tif")
# sdm2024_extcomb <- terra::rast("output/SDM2024_orig/Black Coral/BlackCorals_raster_output_sens_spe/BlackCorals_ext.combinatorial.tif")
# sdm2024_extuni <- terra::rast("output/SDM2024_orig/Black Coral/BlackCorals_raster_output_sens_spe/BlackCorals_ext.univariate.tif")

# terra::plot(sdm2024_extana)
# terra::plot(sdm2024_extcomb)
# terra::plot(sdm2024_extuni)


# Debugging (NULL combinatorial) ----

# Extract predictor values from samples and prediction grid
# sample_preds <- vme_pts[, selected_vme_vars]
# grid_preds   <- extrap_grid[, selected_vme_vars]

# # Calculate covariance matrix from samples
# cov_matrix <- cov(sample_preds, use = "complete.obs")

# # Check if it's invertible before proceeding
# det(cov_matrix)        # should not be zero or near-zero
# kappa(cov_matrix)      # condition number — very large values indicate problems

# # Calculate column means of samples (centre point)
# sample_means <- colMeans(sample_preds, na.rm = TRUE)

# # Invert the covariance matrix
# cov_inv <- solve(cov_matrix)

# # Calculate Mahalanobis distance for each prediction grid point
# maha_distances <- mahalanobis(
#   x      = grid_preds,
#   center = sample_means,
#   cov    = cov_inv,
#   inverted = TRUE   # tells it cov is already inverted
# )

# # Inspect results
# summary(maha_distances)
# hist(maha_distances, breaks = 50)