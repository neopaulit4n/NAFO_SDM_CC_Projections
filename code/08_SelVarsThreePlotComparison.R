# Output three side-by-side plots for baseline/P1 SSP 1/P4 SSP 5 for VMEOI sel vars ----
env_layers <- c(bathy_layers, cmip_layers, unlist(cmip_layers_proj))
env_layer_names <- names(env_layers)
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
    "BCS" = "speed",      # Bottom current speed
    "BStr" = "amp"        # Bottom stress
  )
  return(palette_map[[variable]])
}

if (!dir.exists(paste0(output_folder,"/SelVarsMaps"))) dir.create(paste0(output_folder,"/SelVarsMaps"))

p_sel_cmip_vars <- lapply(vme_var_order[vme_var_order %in% selected_cmip_vars], function(var) {
  p <- ggplot() +
    theme_bw() +
    tidyterra::geom_spatraster(data = env_layers_subset[[var]]) +
    facet_wrap(~ lyr) +
    cmocean::scale_fill_cmocean(
      "",
      name = get_cmocean_palette(str_extract(var, "^(.+?)_", group = 1)),
      na.value = "transparent"
    ) +
    labs(title = gsub("_"," ",var))
  ggsave(paste0(output_folder,"/SelVarsMaps/",var,".jpg"), width = 8, height = 4, dpi = 300)
  return(p)
})
names(p_sel_cmip_vars) <- vme_var_order[vme_var_order %in% selected_cmip_vars]

# Create combined final formatted plot ----
p <- patchwork::wrap_plots(p_sel_cmip_vars, ncol = 2, axes = "collect")
p_height <- ifelse(length(p_sel_cmip_vars) > 6, 8, 6)
ggsave(
  paste0(output_folder,"/SelVarsMaps/CombinedThreePlotComparisons.jpg"), 
  plot = p, width = 10, height = p_height, dpi = 300, scale = 1.2)