# Output environmental variable layers baseline and projected ----
env_layers <- c(bathy_layers, cmip_layers, unlist(cmip_layers_proj))
env_layer_names <- names(env_layers)
# lapply(1:length(env_layers), function(layer) {
#   terra::writeRaster(env_layers[[layer]], filename = paste0("output/env_vars_rasters/",env_layer_names[layer],".tif"))
# })

# Output three side-by-side plots for baseline/P1 SSP 1/P4 SSP 5 for VMEOI sel vars ----
z <- grep(paste0(selected_cmip_vars, collapse = "|"), env_layer_names, value=TRUE)
env_layers_subset <- lapply(selected_cmip_vars, function(var) {
  layers <- env_layers[grep(var, z, value=TRUE)]
  layers <- layers[c(
    grep(paste0("^", var), names(layers)), 
    grep("P1.1|P4.5", names(layers))
  )] |>
    set_names(c("Reference", "P1 SSP 1-2.6", "P4 SSP 5-8.5")) |>
    terra::rast()
}) |>
  set_names(selected_cmip_vars)

get_cmocean_palette <- function(variable) {
  palette_map <- list(
    "SST" = "thermal",    # Sea surface temperature
    "BT" = "thermal",    # Bottom temperature
    "SSS" = "haline",     # Sea surface salinity
    "BS" = "haline",     # Bottom salinity
    "MLD" = "deep",       # Mixed layer depth
    "BCS" = "speed",      # Bottom water velocity
    "BStr" = "amp"        # Bottom stress
  )
  return(palette_map[[variable]])
}

if (!dir.exists(paste0(output_folder,"/SelVarsMaps"))) dir.create(paste0(output_folder,"/SelVarsMaps"))

lapply(selected_cmip_vars, function(var) {
  ggplot() +
    theme_bw() +
    tidyterra::geom_spatraster(data = env_layers_subset[[var]]) +
    facet_wrap(~ lyr) +
    cmocean::scale_fill_cmocean(
      var,
      name = get_cmocean_palette(str_extract(var, "^(.+?)_", group = 1)),
      na.value = "transparent"
    ) +
    labs(title = var)
  ggsave(paste0(output_folder,"/SelVarsMaps/",var,".jpg"), width = 8, height = 4, dpi = 300)
})

