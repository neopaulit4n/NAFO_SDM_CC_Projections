# Create extrapolation plots for SSP 2-4.5

# Read in relevant rasters ----
extrapolation_rasters_mask <- lapply(
  list.files(path = paste0(output_folder, "/Extrapolations/rasters"), pattern = "baseline.+?tif$|2-4.5.+?tif$", full.names = TRUE),
  terra::rast
) |>
  set_names(
    tools::file_path_sans_ext(
      gsub(paste0(vmeoi,"_extrap_"),"",list.files(path = paste0(output_folder, "/Extrapolations/rasters"), pattern = "baseline.+?tif$|2-4.5.+?tif$")))
)

# presenceonly/PA

p <- lapply(c("PA","PresenceOnly"), function(dataset) {

  lapply(c("baseline", period_all), function(poi) {
    r <- extrapolation_rasters_mask[grep(poi, names(extrapolation_rasters_mask))]
    p1 <- ggplot() +
      theme_classic() +
      labs(title = ifelse(poi == "baseline", "Reference", paste(poi, "SSP 2-4.5"))) +
      tidyterra::geom_spatraster(data = r[[grep(paste(dataset,"ExDet","analogue", sep = "\\."),names(r))]]) +
      scale_fill_distiller(
        name = "Analogue",
        palette = "Greys",
        limits = c(0, 1),
        direction = 1,
        na.value = "transparent"
      ) +
      ggnewscale::new_scale_fill() +
      tidyterra::geom_spatraster(data = r[[grep(paste(dataset,"ExDet","univariate", sep = "\\."),names(r))]]) +
      scale_fill_distiller(
        name = "Univariate",
        palette = "Oranges",
        limits = lim_uni,
        direction = 1,
        na.value = "transparent"
      ) +
      ggnewscale::new_scale_fill() +
      {if (length(grep(paste(dataset,"ExDet","combinatorial", sep = "\\."),names(r))) > 0) {
      tidyterra::geom_spatraster(data = r[[grep(paste(dataset,"ExDet","combinatorial", sep = "\\."),names(r))]])
      }} +
      scale_fill_distiller(
        name = "Combinatorial",
        palette = "Greens",
        limits = lim_comb,
        direction = 1,
        na.value = "transparent"
      ) +
      theme(
        legend.box = "horizontal",
        legend.position = "inside",
        legend.position.inside = c(0.75, 0.2),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(5, "mm")
      )
    
    mic_layer <- terra::merge(
      r[[grep(paste(dataset,"mic","analogue", sep = "\\."),names(r))]],
      r[[grep(paste(dataset,"mic","univariate", sep = "\\."),names(r))]]
    )
    {if (length(grep(paste(dataset,"ExDet","combinatorial", sep = "\\."),names(r))) > 0) { 
      mic_layer <- terra::merge(
        mic_layer, 
        r[[grep(paste(dataset,"mic","combinatorial", sep = "\\."),names(r))]]
      )
    }}

    p2 <- ggplot() +
      theme_classic() +
      # labs(title = ifelse(poi == "baseline", "Reference", paste(poi, "SSP 2-4.5"))) +
      tidyterra::geom_spatraster(data = mic_layer) +
      paletteer::scale_fill_paletteer_d(palette = "colorBlindness::paletteMartin", 
        name = "Covariate", na.value = "transparent", na.translate = FALSE) +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.75, 0.22),
        legend.box.just = c("left", "top"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(5, "mm")
        )
    return(list(p1,p2) |> set_names("ExDet","MIC"))

  }) |>
    set_names(c("baseline", period_all)) |>
    unlist()
  }) |>
  set_names(c("PA", "PresenceOnly"))

cowplot::plot_grid(plotlist = p[[1]], nrow = 2, byrow = FALSE)

library(patchwork)
p$PA[[1]] + p$PA[[3]] + p$PA[[5]] + p$PA[[7]] + p$PA[[9]] +
  p$PA[[2]] + p$PA[[4]] + p$PA[[6]] + p$PA[[8]] + p$PA[[10]] +
  patchwork::plot_layout(axes = "collect", axis_titles = "collect", nrow = 2)
ggsave(paste0(output_folder,"/Extrapolations/",vmeoi,"_extrapolations_SSP2-4.5_PA.jpg"), width = 22, height = 11, dpi = 500)

p$PresenceOnly[[1]] + p$PresenceOnly[[3]] + p$PresenceOnly[[5]] + p$PresenceOnly[[7]] + p$PresenceOnly[[9]] +
  p$PresenceOnly[[2]] + p$PresenceOnly[[4]] + p$PresenceOnly[[6]] + p$PresenceOnly[[8]] + p$PresenceOnly[[10]] +
  patchwork::plot_layout(axes = "collect", axis_titles = "collect", nrow = 2)
ggsave(paste0(output_folder,"/Extrapolations/",vmeoi,"_extrapolations_SSP2-4.5_PresenceOnly.jpg"), width = 22, height = 11, dpi = 500)

