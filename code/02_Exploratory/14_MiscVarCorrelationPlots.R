# Correlation plots of bathy vs BS_range and bathy vs BT_range to help explain correlation matrix results
# Requested by Javier for Boltenia 20260728

# Plot correlations between variable pairs over time
baseline_layers <- c(cmip_layers, bathy_layers[vme_terrain_vars_start]) |>
  terra::rast()
# names(baseline_layers) <- c(names(cmip_layers), vme_terrain_vars)

projection_layers <- unlist(cmip_layers_proj, recursive = FALSE) |>
  lapply(X = _, 
    function(layer) {
      rast_result <- c(layer, bathy_layers[vme_terrain_vars_start]) |>
        terra::rast()
      names(rast_result) <- gsub("GEBCO2024_FS005_StudyArea_","", names(rast_result))
      return(rast_result)
    })

all_layers <- list()
all_layers[[1]] <- baseline_layers
names(all_layers) <- "baseline"
all_layers <- c(all_layers, projection_layers)
all_layers_df <- lapply(all_layers, function(x) terra::as.data.frame(x) |> drop_na())

ggplot(data = all_layers_df$baseline, aes(x = FS005, y = BS_range)) +
  geom_point() +
  theme_classic()
ggsave(paste0(output_folder,"/FS005_vs_BS_range.jpg"), dpi = 300)

ggplot(data = all_layers_df$baseline, aes(x = FS005, y = BT_range)) +
  geom_point() +
  theme_classic()
ggsave(paste0(output_folder,"/FS005_vs_BT_range.jpg"), dpi = 300)