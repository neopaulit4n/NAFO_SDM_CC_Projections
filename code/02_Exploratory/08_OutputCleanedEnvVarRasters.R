# Output environmental variable layers baseline and projected ----
env_layers <- c(bathy_layers, cmip_layers, unlist(cmip_layers_proj))
env_layer_names <- names(env_layers)
lapply(1:length(env_layers), function(layer) {
  terra::writeRaster(env_layers[[layer]], filename = paste0("output/env_vars_rasters/",env_layer_names[layer],".tif"))
})
