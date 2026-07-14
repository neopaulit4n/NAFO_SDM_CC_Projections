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

# Create consistent colour scheme for variables in MIC plots ----
mic_pal <- setNames(
  paletteer::palettes_d$colorBlindness$paletteMartin[1:length(c("None",selected_vme_vars))],
  c("None", selected_vme_vars)
)

# Extract limits for univariate and combinatorial legends ----
lim_comb <- lapply(
  grep("ExDet\\.combinatorial", names(extrapolation_rasters_mask)),
  \(x) {
    m <- as.matrix(extrapolation_rasters_mask[[x]])
    m_lim <- c(min(m, na.rm = TRUE), max(m, na.rm = TRUE))
    return(m_lim)})
lim_comb <- c(min(unlist(lapply(lim_comb, "[[", 1))), max(unlist(lapply(lim_comb, "[[", 2))))

lim_uni <- lapply(
  grep("ExDet\\.univariate", names(extrapolation_rasters_mask)),
  \(x) {
    m <- as.matrix(extrapolation_rasters_mask[[x]])
    m_lim <- c(min(m, na.rm = TRUE), max(m, na.rm = TRUE))
    return(m_lim)})
lim_uni <- c(min(unlist(lapply(lim_uni, "[[", 1))), max(unlist(lapply(lim_uni, "[[", 2))))

# Create plots ----
p <- lapply(c("PA","PresenceOnly"), function(dataset) {

  lapply(c("baseline", period_all), function(poi) {
    r <- extrapolation_rasters_mask[grep(poi, names(extrapolation_rasters_mask))]
    p1 <- ggplot() +
      theme_classic() +
      labs(title = ifelse(poi == "baseline", "Reference", paste(poi, "SSP 2-4.5"))) +
      {if (length(grep(paste(dataset,"ExDet","analogue", sep = "\\."),names(r))) > 0) {
      tidyterra::geom_spatraster(data = r[[grep(paste(dataset,"ExDet","analogue", sep = "\\."),names(r))]])
      }} +
      scale_fill_distiller(
        name = "An",
        palette = "Greys",
        limits = c(0, 1),
        direction = 1,
        na.value = "transparent"
      ) +
      ggnewscale::new_scale_fill() +
      {if (length(grep(paste(dataset,"ExDet","univariate", sep = "\\."),names(r))) > 0) {
      tidyterra::geom_spatraster(data = r[[grep(paste(dataset,"ExDet","univariate", sep = "\\."),names(r))]])
      }} +
      scale_fill_distiller(
        name = "Un",
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
        name = "Co",
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
    
    # mic_layer <- terra::merge(
    #   r[[grep(paste(dataset,"mic","analogue", sep = "\\."),names(r))]],
    #   r[[grep(paste(dataset,"mic","univariate", sep = "\\."),names(r))]]
    # )
    # {if (length(grep(paste(dataset,"ExDet","combinatorial", sep = "\\."),names(r))) > 0) { 
    #   mic_layer <- terra::merge(
    #     mic_layer, 
    #     r[[grep(paste(dataset,"mic","combinatorial", sep = "\\."),names(r))]]
    #   )
    # }}

    mic_layer_analogue <- {if (length(grep(paste(dataset,"ExDet","analogue", sep = "\\."),names(r))) > 0) { 
      r[[grep(paste(dataset,"mic","analogue", sep = "\\."),names(r))]]
    }}

    mic_layer_univariate <- {if (length(grep(paste(dataset,"ExDet","univariate", sep = "\\."),names(r))) > 0) { 
      r[[grep(paste(dataset,"mic","univariate", sep = "\\."),names(r))]]
    }}

    mic_layer_combinatorial <- {if (length(grep(paste(dataset,"ExDet","combinatorial", sep = "\\."),names(r))) > 0) { 
      r[[grep(paste(dataset,"mic","combinatorial", sep = "\\."),names(r))]]
    }}

    mic_layer_list <- list(mic_layer_analogue, mic_layer_univariate, mic_layer_combinatorial)
    mic_layer_list <- Filter(Negate(is.null), mic_layer_list)

    mic_layer <- if (length(mic_layer_list) == 0) {
      merged <- NULL
    } else if (length(mic_layer_list) == 1) {
      merged <- mic_layer_list[[1]]  # nothing to merge, use as-is
    } else {
      merged <- do.call(terra::merge, mic_layer_list)
    }

    p2 <- ggplot() +
      theme_classic() +
      # labs(title = ifelse(poi == "baseline", "Reference", paste(poi, "SSP 2-4.5"))) +
      tidyterra::geom_spatraster(data = mic_layer) +
      scale_fill_manual(values = mic_pal, name = "Covariate", na.value = "transparent", na.translate = FALSE) +
      # paletteer::scale_fill_paletteer_d(palette = "colorBlindness::paletteMartin", 
      #   name = "Covariate", na.value = "transparent", na.translate = FALSE) +
      guides(
        fill = guide_legend(
          title.position = "left",
          title.theme = element_text(angle = 90, hjust = 0.5))) +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.25),
        legend.box.just = c("left", "top"),
        legend.title = element_text(size = 8, angle = 90, hjust = 0.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(5, "mm")
        )
    return(list(p1,p2) |> set_names("ExDet","MIC"))

  }) |>
    set_names(c("baseline", period_all)) |>
    unlist()
}) |>
  set_names(c("PA", "PresenceOnly"))

# cowplot::plot_grid(plotlist = p[[1]], nrow = 2, byrow = FALSE)

# Combine using patchwork and export ----
library(patchwork)
patch_p <- p$PA[[1]] + p$PA[[3]] + p$PA[[5]] + p$PA[[7]] + p$PA[[9]] +
  p$PA[[2]] + p$PA[[4]] + p$PA[[6]] + p$PA[[8]] + p$PA[[10]] +
  patchwork::plot_layout(axes = "collect", axis_titles = "collect", nrow = 2)
ggsave(
  paste0(output_folder,"/Extrapolations/",vmeoi,"_extrapolations_SSP2-4.5_PA.jpg"), 
  plot = patch_p, width = 22, height = 11, dpi = 500)

patch_p <- p$PresenceOnly[[1]] + p$PresenceOnly[[3]] + p$PresenceOnly[[5]] + p$PresenceOnly[[7]] + p$PresenceOnly[[9]] +
  p$PresenceOnly[[2]] + p$PresenceOnly[[4]] + p$PresenceOnly[[6]] + p$PresenceOnly[[8]] + p$PresenceOnly[[10]] +
  patchwork::plot_layout(axes = "collect", axis_titles = "collect", nrow = 2)
ggsave(
  paste0(output_folder,"/Extrapolations/",vmeoi,"_extrapolations_SSP2-4.5_PresenceOnly.jpg"), 
  plot = patch_p, width = 22, height = 11, dpi = 500)

