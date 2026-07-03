# Output three side-by-side plots for baseline/P1 SSP 1/P4 SSP 5 for VMEOI sel vars ----
env_layers <- c(bathy_layers, cmip_layers, unlist(cmip_layers_proj))
env_layer_names <- names(env_layers)
z <- grep(paste0(c("BS_max","SSS_range"), collapse = "|"), env_layer_names, value=TRUE)
env_layers_subset <- lapply(selected_cmip_vars, function(var) {
  layers <- env_layers[grep(var, z, value=TRUE)]
  layers <- layers[c(grep(paste0("^", var), names(layers)), grep("2-4\\.5", names(layers)))] |>
    # set_names(c("Reference", "P1 SSP 1-2.6", "P4 SSP 5-8.5")) |>
    terra::rast()
}) |>
  set_names(selected_cmip_vars)

get_cmocean_palette <- function(variable) {
  palette_map <- list(
    "SST" = "thermal",  # Sea surface temperature
    "BT" = "thermal",  # Bottom temperature
    "SSS" = "haline",  # Sea surface salinity
    "BS" = "haline",  # Bottom salinity
    "MLD" = "deep",  # Mixed layer depth
    "BCS" = "speed",  # Bottom current speed
    "BStr" = "amp"  # Bottom stress
  )
  return(palette_map[[variable]])
}

diff1 <- layers[[3]] - layers[[1]]  # difference between P2 and baseline
diff2 <- layers[[5]] - layers[[1]]  # difference between P4 and baseline
diff3 <- layers[[5]] - layers[[3]]
diff_plots <- c(diff1, diff2, diff3)
names(diff_plots) <- c("P2 SSP 2-4.5 - Reference","P4 SSP 2-4.5 - Reference", "P4 SSP 2-4.5 - P2 SSP 2-4.5")


# if (!dir.exists(paste0(output_folder,"/SelVarsMaps"))) dir.create(paste0(output_folder,"/SelVarsMaps"))

# lapply(selected_cmip_vars, function(var) {
#   ggplot() +
#     theme_bw() +
#     tidyterra::geom_spatraster(data = layers) +
#     facet_wrap(~ lyr) +
#     cmocean::scale_fill_cmocean(
#       var,
#       name = get_cmocean_palette(str_extract(var, "^(.+?)_", group = 1)),
#       na.value = "transparent") +
#     labs(title = var)
#   ggsave(paste0(output_folder,"/SelVarsMaps/",var,".jpg"), width = 8, height = 4, dpi = 300)
# })



ggplot() +
  theme_bw() +
  tidyterra::geom_spatraster(data = diff_plots) +
  facet_wrap(~ lyr) +
  scale_fill_gradient2(name = "Difference", low = scales::muted("blue"), high = scales::muted("red"), na.value = "transparent") +
  labs(title = var)
ggsave(paste0(output_folder,"/",vmeoi,"_SSS_range_diff_maps.jpg"), width = 8, height = 4, dpi = 300)
