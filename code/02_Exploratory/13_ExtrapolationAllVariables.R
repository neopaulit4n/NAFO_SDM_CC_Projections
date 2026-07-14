# Extrapolation using all variables (no selection)

cmip_rast <- c(bathy_layers[vme_terrain_vars], compact(cmip_layers_proj[[1]][[2]]))
cmip_rast <- cmip_rast[-grep("_range",names(cmip_rast))] |>
  terra::rast()
prediction_grid <- cmip_rast
output_name <- "P1.2-4.5"

cmip_vars_norange <- cmip_vars[-grep("_range", cmip_vars)]

# Prepare layers ----
# Convert selected vme layers to points
extrap_grid <- terra::as.data.frame(prediction_grid, xy = TRUE) %>%
  drop_na()

## Presence and absence (original) ----
vme_pts_pa <- cmip_comb_df %>%
  filter(VME_Group == vmeoi) %>%
  select(x = Start_Long_DD, y = Start_Lat_DD, all_of(vme_terrain_vars), all_of(cmip_vars_norange)) %>%
  as.data.frame()

## Presence only (refugia) ----
vme_pts_pres <- cmip_comb_df %>%
  filter(VME_Group == vmeoi,
         VME_P_A == "Presence") %>%
  select(x = Start_Long_DD, y = Start_Lat_DD, all_of(vme_terrain_vars), all_of(cmip_vars_norange)) %>%
  as.data.frame()

# Compute extrapolations for each dataset ----
extrapolation_area <- lapply(list(vme_pts_pa, vme_pts_pres), function(dataset) {
  dsmextra::compute_extrapolation(
    samples = dataset,
    covariate.names = c(vme_terrain_vars, cmip_vars_norange),
    prediction.grid = extrap_grid,
    coordinate.system = sp::CRS(SRS_string = "EPSG:4326"))
}) %>%
  set_names(c("PA","PresenceOnly"))

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

# Extrapolation layer maps ----

## Prepare rasters for mapping ----
if (!dir.exists(paste0(output_folder,"/Extrapolations"))) dir.create(paste0(output_folder,"/Extrapolations"))
if (!dir.exists(paste0(output_folder,"/Extrapolations/rasters"))) dir.create(paste0(output_folder,"/Extrapolations/rasters"))

extrapolation_rasters_mask <- lapply(1:length(extrapolation_rasters), function(x) {
  layer <- terra::mask(extrapolation_rasters[[x]], !is.na(extrapolation_rasters[[x]]))
  if (x %in% grep("mic", names(extrapolation_rasters))) {
    layer <- terra::as.factor(layer)
    levels(layer) <- data.frame(id = 0:length(selected_vme_vars), covariate = c("None", selected_vme_vars))
  }
  terra::writeRaster(
    layer, 
    filename = paste0(output_folder,"/Extrapolations/rasters/",vmeoi,"_extrap_",output_name,"_",names(extrapolation_rasters)[x],".tif"),
    overwrite = TRUE
  )
  return(layer)
}) %>%
  set_names(names(extrapolation_rasters))

# Fixing combinatorial layers for P3 SSP 3-7.0 for black corals ----
if (output_name == "P3.3-7.0" & vmeoi == "black_corals") {
  replace_layer <- terra::merge(extrapolation_rasters_mask[[8]], extrapolation_rasters_mask[[10]])
  comp_layer <- cmip_layers_proj[[1]][[1]][[1]] %>%
    terra::crop(replace_layer)

  missing_cells <- terra::logic(replace_layer, comp_layer, oper = "is.na") %>%
    terra::mask(comp_layer)

  extrapolation_rasters_mask[[6]] <- terra::resample(extrapolation_rasters_mask[[6]], missing_cells) %>%
    terra::mask(missing_cells) %>%
    terra::crop(missing_cells)

  extrapolation_rasters_mask[[9]] <- terra::resample(extrapolation_rasters_mask[[9]], missing_cells) %>%
    terra::mask(missing_cells) %>%
    terra::crop(missing_cells)

  rm(replace_layer, comp_layer, missing_cells)
}

# Fixing combinatorial layers for P2 SSP 1-2.6 for sea pens ----
# if (output_name == "P2.1-2.6" & vmeoi == "sea_pens") {
#   replace_layer <- terra::merge(extrapolation_rasters_mask[[1]], extrapolation_rasters_mask[[3]])
#   comp_layer <- cmip_layers_proj[[1]][[1]][[1]] %>%
#     terra::crop(replace_layer)

#   missing_cells <- terra::logic(replace_layer, comp_layer, oper = "is.na") %>%
#     terra::mask(comp_layer)

#   extrapolation_rasters_mask[[2]] <- terra::resample(extrapolation_rasters_mask[[2]], missing_cells) %>%
#     terra::mask(missing_cells) %>%
#     terra::crop(missing_cells)

#   extrapolation_rasters_mask[[5]] <- terra::resample(extrapolation_rasters_mask[[5]], missing_cells) %>%
#     terra::mask(missing_cells) %>%
#     terra::crop(missing_cells)

#   rm(replace_layer, comp_layer, missing_cells)
# }

# dsmextra::map_extrapolation(map.type = "extrapolation", extrapolation.object = extrapolation_area[[1]])

# next()

## Extract limits for univariate and combinatorial legends ----
if (nrow(extrapolation_area[[1]]$data$univariate) > 0 | nrow(extrapolation_area[[2]]$data$univariate) > 0) {
  lim_uni <- c(
    min(c(extrapolation_area[[1]]$data$univariate$ExDet, extrapolation_area[[2]]$data$univariate$ExDet)),
    max(c(extrapolation_area[[1]]$data$univariate$ExDet, extrapolation_area[[2]]$data$univariate$ExDet))  
  )  
}
if (nrow(extrapolation_area[[1]]$data$combinatorial) > 0 | nrow(extrapolation_area[[2]]$data$combinatorial) > 0) {
  lim_comb <- c(
    min(c(extrapolation_area[[1]]$data$combinatorial$ExDet, extrapolation_area[[2]]$data$combinatorial$ExDet)),
    max(c(extrapolation_area[[1]]$data$combinatorial$ExDet, extrapolation_area[[2]]$data$combinatorial$ExDet))  
  )
}

## Generate ExDet maps ----
library(patchwork)
extrap_exdet_maps <- lapply(c("PA","PresenceOnly"), function(dataset) {
  ggplot() +
    theme_classic() +
    labs(title = ifelse(dataset == "PA", "Presence + Absence", "Presence Only")) +
    {if (length(grep(paste(dataset,"ExDet","analogue", sep = "\\."),names(extrapolation_rasters_mask))) > 0) {
    tidyterra::geom_spatraster(data = extrapolation_rasters_mask[[grep(paste(dataset,"ExDet","analogue", sep = "\\."),names(extrapolation_rasters_mask))]])
    }} +
    scale_fill_distiller(
      name = "Analogue",
      palette = "Greys",
      limits = c(0, 1),
      direction = 1,
      na.value = "transparent"
    ) +
    ggnewscale::new_scale_fill() +
    {if (length(grep(paste(dataset,"ExDet","univariate", sep = "\\."),names(extrapolation_rasters_mask))) > 0) {
    tidyterra::geom_spatraster(data = extrapolation_rasters_mask[[grep(paste(dataset,"ExDet","univariate", sep = "\\."),names(extrapolation_rasters_mask))]])
    }} +
    scale_fill_distiller(
      name = "Univariate",
      palette = "Oranges",
      limits = lim_uni,
      direction = 1,
      na.value = "transparent"
    ) +
    ggnewscale::new_scale_fill() +
    {if (length(grep(paste(dataset,"ExDet","combinatorial", sep = "\\."),names(extrapolation_rasters_mask))) > 0) {
    tidyterra::geom_spatraster(data = extrapolation_rasters_mask[[grep(paste(dataset,"ExDet","combinatorial", sep = "\\."),names(extrapolation_rasters_mask))]])
    }} +
    scale_fill_distiller(
      name = "Combinatorial",
      palette = "Greens",
      limits = lim_comb,
      direction = 1,
      na.value = "transparent"
    ) +
    theme(
      # axis.text.y = element_blank(),
      # axis.ticks.y = element_blank(),
      # axis.title.y = element_blank(),
      legend.box = "horizontal",
      legend.position = "inside",
      legend.position.inside = c(0.75, 0.2),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.key.size = unit(5, "mm")
    )
})

# Use patchwork to create combined plot
# extrap_exdet_maps <- extrap_exdet_maps[[1]] + extrap_exdet_maps[[2]] + 
#   patchwork::plot_layout(axes = "collect")

# Save
# ggsave(paste0("output/03_RF_Map_Outputs/",vmeoi,"_extrapolations_ExDet_",poi,"_",sspoi,".jpg"), 
#   plot = extrap_exdet_maps,
#   width = 10, height = 5, dpi = 300)

## Create consistent colour scheme for variables in MIC plots ----
mic_pal <- setNames(
  paletteer::palettes_d$colorBlindness$paletteMartin[1:length(c("None",selected_vme_vars))],
  c("None", selected_vme_vars)
)

## Generate MIC maps ----
extrap_mic_maps <- lapply(c("PA","PresenceOnly"), function(dataset) {

  mic_layer_analogue <- {if (length(grep(paste(dataset,"ExDet","analogue", sep = "\\."),names(extrapolation_rasters_mask))) > 0) { 
      extrapolation_rasters_mask[[grep(paste(dataset,"mic","analogue", sep = "\\."),names(extrapolation_rasters_mask))]]
  }}

  mic_layer_univariate <- {if (length(grep(paste(dataset,"ExDet","univariate", sep = "\\."),names(extrapolation_rasters_mask))) > 0) { 
      extrapolation_rasters_mask[[grep(paste(dataset,"mic","univariate", sep = "\\."),names(extrapolation_rasters_mask))]]
  }}

  mic_layer_combinatorial <- {if (length(grep(paste(dataset,"ExDet","combinatorial", sep = "\\."),names(extrapolation_rasters_mask))) > 0) { 
      extrapolation_rasters_mask[[grep(paste(dataset,"mic","combinatorial", sep = "\\."),names(extrapolation_rasters_mask))]]
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

  ggplot() +
    theme_classic() +
    labs(title = ifelse(dataset == "PA", "Presence + Absence", "Presence Only")) +
    tidyterra::geom_spatraster(data = mic_layer) +
    # tidyterra::geom_spatraster(data = extrapolation_rasters_mask[[grep(paste(dataset,"mic","analogue", sep = "\\."),names(extrapolation_rasters_mask))]]) +
    # tidyterra::geom_spatraster(data = extrapolation_rasters_mask[[grep(paste(dataset,"mic","univariate", sep = "\\."),names(extrapolation_rasters_mask))]]) +
    # {if (length(grep(paste(dataset,"ExDet","combinatorial", sep = "\\."),names(extrapolation_rasters_mask))) > 0) {     
    #   tidyterra::geom_spatraster(data = extrapolation_rasters_mask[[grep(paste(dataset,"mic","combinatorial", sep = "\\."),names(extrapolation_rasters_mask))]])
    # }} +
    # paletteer::scale_fill_paletteer_d(palette = "colorBlindness::paletteMartin", 
    #   name = "Covariate", na.value = "transparent", na.translate = FALSE) +
    scale_fill_manual(values = mic_pal, name = "Covariate", na.value = "transparent", na.translate = FALSE) +
    guides(
      fill = guide_legend(
        title.position = "left",
        title.theme = element_text(angle = 90, hjust = 0.5))) +
    # if (dataset == "PresenceOnly") {
    theme(
      # axis.text.y = element_blank(),
      # axis.ticks.y = element_blank(),
      # axis.title.y = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.85, 0.25),
      legend.box.just = c("left", "top"),
      legend.title = element_text(size = 8, angle = 90, hjust = 0.5),
      legend.text = element_text(size = 8),
      legend.key.size = unit(5, "mm")
      )
    # } else {
      # theme(legend.position = "none")
    # }
})

# Use patchwork to create combined plot
# extrap_mic_maps <- extrap_mic_maps[[1]] + extrap_mic_maps[[2]] +
#   patchwork::plot_layout(axes = "collect")

# Save
# ggsave(paste0("output/03_RF_Map_Outputs/",vmeoi,"_extrapolations_MIC_",poi,"_",sspoi,".jpg"), 
#   plot = extrap_mic_maps,
#   width = 10, height = 5, dpi = 300)

# Save both together
extrap_exdet_maps[[1]] + extrap_exdet_maps[[2]] + extrap_mic_maps[[1]] + extrap_mic_maps[[2]] +
  patchwork::plot_layout(axes = "collect")
ggsave(paste0(output_folder,"/Extrapolations/",vmeoi,"_extrapolations_",output_name,".jpg"), 
  width = 10, height = 10, dpi = 300)


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

write_csv(extrap_analysis, paste0(output_folder,"/Extrapolations/",vmeoi,"_extrapolation_percentages_",output_name,".csv"))


# Overlay univariate extrapolation layer with MaxClass map layer for both PA and PresenceOnly datasets ----
# extrap_uni <- extrapolation_rasters$PA.ExDet.univariate
# extrap_uni <- terra::ifel(!is.na(extrap_uni), 1, NA) %>%
#   terra::extend(rf_pred_all$MaxClass$`1-2.6_P1`) %>%
#   terra::subst(NA, 0)
# extrap_uni_maxclass <- rf_pred_all$MaxClass$`1-2.6_P1` + extrap_uni * 2
# levels(extrap_uni_maxclass) <- data.frame(
#   value = 0:3,
#   label = c(
#     "Absence (not extrapolated)",
#     "Presence (not extrapolated)",
#     "Absence (extrapolated)",
#     "Presence (extrapolated)"
#   )
# )

# ggplot() +
#     theme_classic() +    
#     tidyterra::geom_spatraster(data = extrap_uni_maxclass, na.rm = TRUE) +
#     scale_fill_manual(
#       values = c(
#         "Absence (not extrapolated)" = "#ffebcd", 
#         "Presence (not extrapolated)" = "#b87333", 
#         "Absence (extrapolated)" = "coral", 
#         "Presence (extrapolated)" = "coral4"
#       ),
#       na.value = "transparent",
#       na.translate = FALSE
#     ) +
#     # geom_contour(data = bathy_noaa, 
#     #   aes(x = x, y = y, z = z, fill = NULL), 
#     #   breaks = seq(from = -50, to = -5000, by = -250),
#     #   color = "darkgrey", 
#     #   linewidth = 0.3, 
#     #   alpha = 0.4) +
#     theme(legend.position = "right",
#           legend.title = element_blank(),
#           axis.title = element_blank()) +
#     scale_x_continuous(expand = c(0,0)) +
#     scale_y_continuous(expand = c(0,0))
# ggsave(paste0("output/03_RF_Map_Outputs/",vmeoi,"_MaxClass_UnivariateExtrapOverlay_",poi,"_",sspoi,".jpg"), 
#   # plot = extrap_mic_maps,
#   width = 5, height = 5, dpi = 300)

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