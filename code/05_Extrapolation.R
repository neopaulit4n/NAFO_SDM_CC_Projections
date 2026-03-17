
# Extrapolations with dsmextra

cat("Computing extrapolations...")
# samples = presence data points with associated covariate values (known) - points
# prediction.grid = grid of covariate values across the study area for which presence is to be predicted (unknown) - points

# Prepare layers ----
# Convert selected vme layers to points
extrap_grid <- terra::as.data.frame(vme_layers_rast, xy = TRUE) %>%
  drop_na()

## Presence and absence (original) ----
vme_pts_pa <- cmip_comb_df %>%
  filter(VME_Group == vmeoi) %>%
  select(x = Start_Long_DD, y = Start_Lat_DD, all_of(selected_vme_vars)) %>%
  as.data.frame()

## Presence only (refugia) ----
vme_pts_pres <- cmip_comb_df %>%
  filter(VME_Group == vmeoi,
         VME_P_A == "Presence") %>%
  select(x = Start_Long_DD, y = Start_Lat_DD, all_of(selected_vme_vars)) %>%
  as.data.frame()

# Compute extrapolations for each dataset ----
# extrapolation_area <- dsmextra::compute_extrapolation(
#   samples = vme_pts,
#   covariate.names = selected_vme_vars,
#   prediction.grid = extrap_grid,
#   coordinate.system = sp::CRS(SRS_string = "EPSG:4326"))

extrapolation_area <- lapply(list(vme_pts_pa, vme_pts_pres), function(dataset) {
  dsmextra::compute_extrapolation(
    samples = dataset,
    covariate.names = selected_vme_vars,
    prediction.grid = extrap_grid,
    coordinate.system = sp::CRS(SRS_string = "EPSG:4326"))
})


# Extract extrapolation rasters ----
extrapolation_rasters <- lapply(extrapolation_area, function(extrap) {
  # Determine which extrapolation types are not null
  extrap_types <- c("univariate", "combinatorial", "analogue")[which(sapply(extrap$data, nrow)[2:4] > 0)]
  lapply(c("ExDet", "mic"), function(method) {
    lapply(extrap_types, function(type) {
      terra::rast(extrap$rasters[[method]][[type]])
    }) %>%
      set_names(extrap_types)
  }) %>%
    set_names(c("ExDet", "mic")) %>%
    unlist(recursive = FALSE)
}) %>%
  set_names(c("PA","PresenceOnly")) %>%
  unlist()

# Generate maps ----
# extrap2plot <- extrapolation_rasters[-grep("mic.analogue", names(extrapolation_rasters))]
extrapolation_maps <- lapply(extrap2plot, function(x) {
  ggplot() +
    tidyterra::geom_spatraster(data = x) +
    scale_fill_viridis_c(na.value = "transparent") +
    labs(title = names(extrap2plot)) +
    theme_minimal()  
})

extrapolation_map_grid <- cowplot::plot_grid(plotlist = extrapolation_maps, nrow = 2)
ggsave(paste0("output/03_RF_Map_Outputs/",vmeoi,"_extrapolations_",poi,"_",sspoi,".jpg"), 
  plot = extrapolation_map_grid,
  width = 10, height = 5, dpi = 300)

# Extrapolation analysis (extract raster percentages) ----
extrap_analysis <- lapply(list(vme_pts_pa, vme_pts_pres), function(dataset) {
  x <- dsmextra::extrapolation_analysis(
    samples = dataset,
    covariate.names = selected_vme_vars,
    prediction.grid = extrap_grid,
    coordinate.system = sp::CRS(SRS_string = "EPSG:4326"),
    nearby.compute = FALSE,
    map.generate = FALSE
  )
  temp_extrap <- x$extrapolation$summary$extrapolation %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), 
      cols_vary = "slowest",
      names_sep = "\\.",
      names_to = c("Type","metric")) %>%
    pivot_wider(names_from = "metric", values_from = "value") %>%
    rename(freq = n, perc = p) %>%
    mutate(Type = str_to_sentence(Type),
      covariate = "Overall")
  temp_mic <- x$extrapolation$summary$mic %>%
    bind_rows()
  temp <- bind_rows(temp_extrap, temp_mic)
  return(temp)
}) %>%
  set_names("PA", "PresenceOnly") %>%
  bind_rows(.id = "InputData")

write_csv(extrap_analysis, paste0("output/02_Modelling_Outputs/",vmeoi,"/",vmeoi,"_",poi,"_",sspoi,"_extrapolation_percentages.csv"))



# Load previous SDM2024 extrapolation maps for comparison ----

# Read TIF
# sdm2024_extana <- terra::rast("output/SDM2024_orig/Black Coral/BlackCorals_raster_output_sens_spe/BlackCorals_ext.analogue.tif")
# sdm2024_extcomb <- terra::rast("output/SDM2024_orig/Black Coral/BlackCorals_raster_output_sens_spe/BlackCorals_ext.combinatorial.tif")
# sdm2024_extuni <- terra::rast("output/SDM2024_orig/Black Coral/BlackCorals_raster_output_sens_spe/BlackCorals_ext.univariate.tif")
# sdm2024_micana <- terra::rast("output/SDM2024_orig/Sea Pens/SeapensVME_mic.analogue.tif")

# terra::plot(sdm2024_extana)
# terra::plot(sdm2024_extcomb)
# terra::plot(sdm2024_extuni)
# terra::plot(sdm2024_micana)

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