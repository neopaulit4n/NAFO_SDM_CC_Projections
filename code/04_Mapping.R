
# Load mapping layers ----

# # Load NOAA bathymetry layer for contours
bathy_noaa <- readRDS("data/raw/Mapping_Layers/bathy_noaa.rds")

# # Load Canadian EEZ boundary from shapefile
# # eez <- sf::st_read("data/raw/Mapping_Layers/eez/eez.shp") %>%
# #   sf::st_as_sf(crs = st_crs(4326)) %>%
# #   sf::st_make_valid() %>%
# #   sf::st_union() %>%
# #   sf::st_cast("POLYGON") %>%
# #   sf::st_cast("LINESTRING") %>%
# #   sf::st_crop(xmin = sa_lims[1]-0.09, xmax = sa_lims[2]+0.09, 
# #               ymin = sa_lims[3]-0.09, ymax = sa_lims[4]+0.09)

# # Load small fishing footprint
# footprint <- sf::st_read("data/raw/Mapping_Layers/FootprintProjectedShp/FootprintAreaProjected.shp") %>%
#   sf::st_as_sf(crs = st_crs(4326)) %>%
#   sf::st_make_valid() %>%
#   sf::st_union() %>%
#   sf::st_cast("LINESTRING")


# Map of overall area and zoomed in plot of NAFO boundary ----
# ggplot() +
#   theme_classic() +
#   # geom_sf(data = sa, aes(colour = "NAFO Study Area"), fill = "lightblue", alpha = 0.8, show.legend = "line") +
#   tidyterra::geom_spatraster(data = rf_pred_pa, na.rm = TRUE) +
#   scale_fill_manual(values = c("0" = "#ffebcd", "1" = "#b87333"),
#                     na.value = "transparent",
#                     na.translate = FALSE,  # remove NAs from legend
#                     labels = c("0" = "Absence", "1" = "Presence")) +
#   geom_sf(data = footprint, aes(colour = "NAFO Fishing Footprint"), show.legend = "line") +  
#   # geom_sf(data = eez, aes(colour = "Canadian EEZ"), inherit.aes = FALSE, linewidth = 1, show.legend = "line") +
#   # Adjust colours
#   scale_colour_manual(name = "Boundary", 
#                       values = c("NAFO Study Area" = "black", "NAFO Fishing Footprint" = "blue", "Canadian EEZ" = "red")) +
#   geom_contour(data = bathy_noaa, 
#                aes(x = x, y = y, z = z, fill = NULL), 
#                breaks = seq(from = -50, to = -5000, by = -250),
#                color = "darkgrey", 
#                linewidth = 0.3, 
#                alpha = 0.4) +
#   coord_sf(xlim = terra::ext(rf_pred_pa)[1:2], 
#            ylim = terra::ext(rf_pred_pa)[3:4], expand = FALSE) +
#   # geom_vline(xintercept = seq(from = terra::ext(rf_pred_pa)[1], 
#   #                             to = terra::ext(rf_pred_pa)[2], 
#   #                             by = terra::res(rf_pred_pa)[1]),
#   #            color = "black", linewidth = 0.1, alpha = 0.1) +
#   # geom_hline(yintercept = seq(from = terra::ext(rf_pred_pa)[3], 
#   #                             to = terra::ext(rf_pred_pa)[4], 
#   #                             by = terra::res(rf_pred_pa)[2]),
#   #            color = "black", linewidth = 0.1, alpha = 0.1) +
#   labs(title = paste("Predicted Presence/Absence for", vme_group),
#        fill = "Prediction", x = "Longitude", y = "Latitude")

# Grid of plots with periods as columns and SSPs as rows for each metric ----

metric_names <- c("MaxClass", "MaxClassF", "MaxClassAvgProb", "CombConf", "CVSum")

## Read in rasters ----
rf_pred_future_all <- lapply(metric_names, function(metric) {
  metric_pred_names <- list.files(paste0("output/02_Modelling_Outputs/", vmeoi), pattern = paste0("rf_res_future_", metric, "_"), full.names = TRUE)
  metric_preds <- terra::rast(metric_pred_names)
  names(metric_preds) <- paste0(str_extract(metric_pred_names, "1-2.6|2-4.5|3-7.0|5-8.5"), "_", str_extract(metric_pred_names, "P[1-4]"))
  metric_preds <- metric_preds[[order(names(metric_preds))]]  # reorder layers by period (P1-P4) within each SSP for facetting
  
  # Factorise MaxClass rasters for plotting
  if (metric == "MaxClass") {
    metric_preds <- terra::as.factor(metric_preds)
  }

  return(metric_preds)
}) %>%
  set_names(metric_names)

rf_pred_current_all <- lapply(metric_names, function(metric) {
  metric_pred_names <- list.files(paste0("output/02_Modelling_Outputs/", vmeoi), pattern = paste0("rf_res_current_", metric, "\\.tif"), full.names = TRUE)
  metric_preds <- terra::rast(metric_pred_names)
  # names(metric_preds) <- paste0(str_extract(metric_pred_names, "1-2.6|2-4.5|3-7.0|5-8.5"), "_", str_extract(metric_pred_names, "P[1-4]"))
  # metric_preds <- metric_preds[[order(names(metric_preds))]]  # reorder layers by period (P1-P4) within each SSP for facetting
  
  # Factorise MaxClass rasters for plotting
  if (metric == "MaxClass") {
    metric_preds <- terra::as.factor(metric_preds)
  }

  return(metric_preds)
}) %>%
  set_names(metric_names)

## Calculate area of predicted presence (MaxClass) ----

# r <- rf_pred_all$MaxClass$`1-2.6_P1`

# # Subset raster to presence cells only (mask out absence)
# presence_raster <- terra::ifel(r == 1, r, NA)

# # Calculate area of each cell in km²
# area_raster <- terra::cellSize(presence_raster, unit = "km", mask = TRUE)

# # Sum all presence cell areas to get total area in km²
# total_area_km2 <- terra::global(area_raster, "sum", na.rm = TRUE)
# print(total_area_km2)

compute_presence_areas <- function(metric) {
  rast_stack <- rf_pred_all[[metric]]
  
  layer_names <- names(rast_stack)
  
  areas <- sapply(layer_names, function(lyr_name) {
    r_layer <- rast_stack[[lyr_name]]
    presence <- terra::ifel(r_layer == 1, r_layer, NA)
    area_rast <- terra::cellSize(presence, unit = "km", mask = TRUE)
    terra::global(area_rast, "sum", na.rm = TRUE)$sum
  })
  
  data.frame(
    lyr = layer_names,
    label = paste0(format(round(areas, 0), big.mark = ","), " km²")
  )
}


## Generate plots ----

# Try alternative method using facet_wrap for each metric
rf_pred_maps <- lapply(metric_names, function(metric) {
  
  # Compute area labels for MaxClass only
  # area_label_layer <- if (metric == "MaxClass") {
  #   area_df <- compute_presence_areas(metric)
  #   geom_label(
  #     data = area_df,
  #     aes(label = label),
  #     x = -50, y = 48, 
  #     hjust = 1.05, vjust = -0.5,
  #     size = 3,
  #     fill = alpha("white", 0.7),
  #     label.size = NA,
  #     inherit.aes = FALSE
  #   )
  # } else {
  #   NULL  # ggplot silently ignores NULL layers
  # }

  # Define fill scale based on metric
  ggtheme_metric <- switch(metric,
    "MaxClass" = function() {
      scale_fill_manual(values = c("0" = "#ffebcd", "1" = "#b87333"),
                    na.value = "transparent",
                    na.translate = FALSE,  # remove NAs from legend
                    labels = c("0" = "Absence", "1" = "Presence"))
    },
    "MaxClassF" = function() {
      scale_fill_binned(breaks = c(0.5,0.6,0.8,0.9,1),
        palette = c("#b06500", "#e5aa70", "#96c8a2", "#008b8b"),
        guide = guide_coloursteps(),
        na.value = "transparent") 
    },
    "MaxClassAvgProb" = function() {scale_fill_gradient2(low = "#1164b4", mid = "#ffff99", high = "#e03c31", midpoint = 0.5, na.value = "transparent")},
    "CombConf" = function() {scale_fill_continuous(palette = "YlGn", na.value = "transparent")},
    "CVSum" = function() {scale_fill_continuous(palette = "YlGn", na.value = "transparent")}
  )

  # Create plot
  # p <- ggplot() +
  #   theme_classic() +    
  #   tidyterra::geom_spatraster(data = rf_pred_future_all[[metric]], na.rm = TRUE) +
  #   facet_wrap(~ lyr, ncol = 4) +
  #   ggtheme_metric() +
  #   geom_contour(data = bathy_noaa, 
  #     aes(x = x, y = y, z = z, fill = NULL), 
  #     breaks = seq(from = -50, to = -5000, by = -250),
  #     color = "darkgrey", 
  #     linewidth = 0.3, 
  #     alpha = 0.4) +
  #   # Add area label per facet
  #   # area_label_layer +
  #   theme(legend.position = "right",
  #         legend.title = element_blank(),
  #         axis.title = element_blank()) +
  #   scale_x_continuous(expand = c(0,0)) +
  #   scale_y_continuous(expand = c(0,0))

  # ggsave(paste0("output/03_RF_Map_Outputs/", vmeoi, "_", metric, "_facet.jpg"), p,
  #   width = 10, height = 10, dpi = 300)

  p_current <- ggplot() +
    theme_classic() +    
    tidyterra::geom_spatraster(data = rf_pred_current_all[[metric]], na.rm = TRUE) +
    ggtheme_metric() +
    geom_contour(data = bathy_noaa, 
      aes(x = x, y = y, z = z, fill = NULL), 
      breaks = seq(from = -50, to = -5000, by = -250),
      color = "darkgrey", 
      linewidth = 0.3, 
      alpha = 0.4) +
    # Add area label per facet
    # area_label_layer +
    theme(legend.position = "right",
          legend.title = element_blank(),
          axis.title = element_blank()) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0))

  ggsave(paste0("output/03_RF_Map_Outputs/", vmeoi, "_current_", metric, ".jpg"), p_current,
    width = 10, height = 10, dpi = 300)
})



# Comparisons with previous SDM2024 rasters ----
# ## Load previous tiff files to compare ----
# sdm2024_pred_stack <- c(
#   MaxClass = terra::rast(filter(sdm2024_raster_output_df, VME_group == vmeoi)$MaxClass) %>%
#     terra::as.factor() %>%
#     terra::project("EPSG:4326"),
#   MaxClassF = terra::rast(filter(sdm2024_raster_output_df, VME_group == vmeoi)$MaxClassF) %>%
#     terra::project("EPSG:4326"),
#   # AvgProb = terra::rast(...),
#   MaxClassAvgProb = terra::rast(filter(sdm2024_raster_output_df, VME_group == vmeoi)$MaxClassAvgProb) %>%
#     terra::project("EPSG:4326"),
#   CombConf = terra::rast(filter(sdm2024_raster_output_df, VME_group == vmeoi)$CombConf) %>%
#     terra::project("EPSG:4326"),
#   CVSum = terra::rast(filter(sdm2024_raster_output_df, VME_group == vmeoi)$CVSum) %>%
#     terra::project("EPSG:4326")
#   )

# # Define plot extents
# plot_xlim <- terra::ext(fold_predictions_spatial_reclass[[1]])[1:2]
# plot_ylim <- terra::ext(fold_predictions_spatial_reclass[[1]])[3:4]

# ## Most frequent class (MaxClass) ----
# plot_sdm2024_MaxClass <- ggplot() +
#   theme_classic() +
#   tidyterra::geom_spatraster(data = sdm2024_pred_stack$MaxClass, na.rm = TRUE) +
#   scale_fill_manual(values = c("0" = "#ffebcd", "1" = "#b87333"),
#                     na.value = "transparent",
#                     na.translate = FALSE,  # remove NAs from legend
#                     labels = c("0" = "Absence", "1" = "Presence")) +
#   geom_contour(data = bathy_noaa, 
#                aes(x = x, y = y, z = z, fill = NULL), 
#                breaks = seq(from = -50, to = -5000, by = -250),
#                color = "darkgrey", 
#                linewidth = 0.3, 
#                alpha = 0.4) +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         axis.title = element_blank()) +
#   scale_x_continuous(limits = plot_xlim, expand = c(0,0)) +
#   scale_y_continuous(limits = plot_ylim, expand = c(0,0)) +
#   labs(title = paste0("SDM2024\n",vmeoi,"\nMaxClass"))

# plot_new_MaxClass <- ggplot() +
#   theme_classic() +
#   tidyterra::geom_spatraster(data = rf_pred_all$`black_corals_P1_1-2.6_rf_res_MaxClass.tif`, na.rm = TRUE) +
#   scale_fill_manual(values = c("0" = "#ffebcd", "1" = "#b87333"),
#                     na.value = "transparent",
#                     na.translate = FALSE,  # remove NAs from legend
#                     labels = c("0" = "Absence", "1" = "Presence")) +
#   geom_contour(data = bathy_noaa, 
#                aes(x = x, y = y, z = z, fill = NULL), 
#                breaks = seq(from = -50, to = -5000, by = -250),
#                color = "darkgrey", 
#                linewidth = 0.3, 
#                alpha = 0.4) +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         axis.title = element_blank()) +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0)) +
#   labs(title = paste("New | Period",poi,"| SSP",sspoi,"\n",vmeoi,"\nMaxClass"))

# cowplot::plot_grid(plot_sdm2024_MaxClass, plot_new_MaxClass, 
#                    # labels = c("SDM2024 MaxClass", "New MaxClass"), 
#                    ncol = 2,
#                    align = "hv")
# ggsave("output/01_BlackCorals_RasterMetrics_OldVSNew/MaxClass_comparison.png", 
#                 width = 10, height = 5, dpi = 300)

# ## Frequency of most frequent class (fraction of runs) ----
# plot_sdm2024_MaxClassF <- ggplot() +
#   theme_classic() +
#   tidyterra::geom_spatraster(data = sdm2024_pred_stack$MaxClassF, na.rm = TRUE) +
#   scale_fill_binned(
#     breaks = c(0.5,0.6,0.8,0.9,1),
#     palette = c("#b06500", "#e5aa70", "#96c8a2", "#008b8b"),
#     guide = guide_coloursteps(),
#     na.value = "transparent"
#   ) +
#   geom_contour(data = bathy_noaa, 
#                aes(x = x, y = y, z = z, fill = NULL), 
#                breaks = seq(from = -50, to = -5000, by = -250),
#                color = "darkgrey", 
#                linewidth = 0.3, 
#                alpha = 0.4) +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         axis.title = element_blank()) +
#   scale_x_continuous(limits = plot_xlim, expand = c(0,0)) +
#   scale_y_continuous(limits = plot_ylim, expand = c(0,0)) +
#   labs(title = paste0("SDM2024\n",vmeoi,"\nMaxClassF"))

# plot_new_MaxClassF <- ggplot() +
#   theme_classic() +
#   tidyterra::geom_spatraster(data = rf_pred_comp$MaxClassF, na.rm = TRUE) +
#   scale_fill_binned(
#     breaks = c(0.5,0.6,0.8,0.9,1),
#     palette = c("#b06500", "#e5aa70", "#96c8a2", "#008b8b"),
#     guide = guide_coloursteps(),
#     na.value = "transparent"
#   ) +
#   geom_contour(data = bathy_noaa, 
#                aes(x = x, y = y, z = z, fill = NULL), 
#                breaks = seq(from = -50, to = -5000, by = -250),
#                color = "darkgrey", 
#                linewidth = 0.3, 
#                alpha = 0.4) +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         axis.title = element_blank()) +
#   scale_x_continuous(limits = plot_xlim, expand = c(0,0)) +
#   scale_y_continuous(limits = plot_ylim, expand = c(0,0)) +
#   labs(title = paste("New | Period",poi,"| SSP",sspoi,"\n",subsample_absences,"|",keep_all_cmip_vars,"\n",vmeoi,"\nMaxClassF"))

# cowplot::plot_grid(plot_sdm2024_MaxClassF, plot_new_MaxClassF, 
#                    # labels = c("SDM2024 MaxClassF", "New MaxClassF"), 
#                    ncol = 2,
#                    align = "hv")
# ggsave("output/01_BlackCorals_RasterMetrics_OldVSNew/MaxClassF_comparison.png", 
#        width = 10, height = 5, dpi = 300)

# ## Average probability of classes ----

# # plot_new_AvgProb_Abs <- ggplot() +
# #   theme_classic() +
# #   tidyterra::geom_spatraster(data = rf_pred_comp$AvgProb[[1]], na.rm = TRUE) +
# #   scale_fill_continuous(palette = c("white","blue")) +
# #   theme(legend.position = "bottom",
# #         legend.title = element_blank())
# # plot_new_AvgProb_Pres <- ggplot() +
# #   theme_classic() +
# #   tidyterra::geom_spatraster(data = rf_pred_comp$AvgProb[[2]], na.rm = TRUE) +
# #   scale_fill_continuous(palette = c("white","blue")) +
# #   theme(legend.position = "bottom",
# #         legend.title = element_blank())
# # cowplot::plot_grid(plot_new_AvgProb_Abs, plot_new_AvgProb_Pres, 
# #                    labels = c("New AvgProb Absence", "New AvgProb Presence"), 
# #                    ncol = 2)
# # 
# # cowplot::plot_grid(plot_sdm2024_AvgProb, plot_new_AvgProb, 
# #                    labels = c("SDM2024 AvgProb", "New AvgProb"), 
# #                    ncol = 2)

# ## Average probability of maximum frequency class ----
# plot_sdm2024_MaxClassAvgProb <- ggplot() +
#   theme_classic() +
#   tidyterra::geom_spatraster(data = sdm2024_pred_stack$MaxClassAvgProb, na.rm = TRUE) +
#   scale_fill_gradient2(low = "#1164b4", mid = "#ffff99", high = "#e03c31", midpoint = 0.5, na.value = "transparent") +
#   geom_contour(data = bathy_noaa, 
#                aes(x = x, y = y, z = z, fill = NULL), 
#                breaks = seq(from = -50, to = -5000, by = -250),
#                color = "darkgrey", 
#                linewidth = 0.3, 
#                alpha = 0.4) +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         axis.title = element_blank()) +
#   scale_x_continuous(limits = plot_xlim, expand = c(0,0)) +
#   scale_y_continuous(limits = plot_ylim, expand = c(0,0)) +
#   labs(title = paste0("SDM2024\n",vmeoi,"\nMaxClassAvgProb"))

# plot_new_MaxClassAvgProb <- ggplot() +
#   theme_classic() +
#   tidyterra::geom_spatraster(data = rf_pred_comp$MaxClassAvgProb, na.rm = TRUE) +
#   scale_fill_gradient2(low = "#1164b4", mid = "#ffff99", high = "#e03c31", midpoint = 0.5, na.value = "transparent") +
#   geom_contour(data = bathy_noaa, 
#                aes(x = x, y = y, z = z, fill = NULL), 
#                breaks = seq(from = -50, to = -5000, by = -250),
#                color = "darkgrey", 
#                linewidth = 0.3, 
#                alpha = 0.4) +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         axis.title = element_blank()) +
#   scale_x_continuous(limits = plot_xlim, expand = c(0,0)) +
#   scale_y_continuous(limits = plot_ylim, expand = c(0,0)) +
#   labs(title = paste("New | Period",poi,"| SSP",sspoi,"\n",subsample_absences,"|",keep_all_cmip_vars,"\n",vmeoi,"\nMaxClassAvgProb"))
  
# cowplot::plot_grid(plot_sdm2024_MaxClassAvgProb, plot_new_MaxClassAvgProb, 
#                    # labels = c("SDM2024 MaxClassAvgProb", "New MaxClassAvgProb"), 
#                    align = "hv",
#                    ncol = 2)
# ggsave("output/01_BlackCorals_RasterMetrics_OldVSNew/MaxClassAvgProb_comparison.png", 
#        width = 10, height = 5, dpi = 300)

# ## Combined confidence metric ----
# plot_sdm2024_CombConf <- ggplot() +
#   theme_classic() +
#   tidyterra::geom_spatraster(data = sdm2024_pred_stack$CombConf, na.rm = TRUE) +
#   scale_fill_continuous(palette = "YlGn", na.value = "transparent") +
#   geom_contour(data = bathy_noaa, 
#                aes(x = x, y = y, z = z, fill = NULL), 
#                breaks = seq(from = -50, to = -5000, by = -250),
#                color = "darkgrey", 
#                linewidth = 0.3, 
#                alpha = 0.4) +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         axis.title = element_blank()) +
#   scale_x_continuous(limits = plot_xlim, expand = c(0,0)) +
#   scale_y_continuous(limits = plot_ylim, expand = c(0,0)) +
#   labs(title = paste0("SDM2024\n",vmeoi,"\nCombConf"))

# plot_new_CombConf <- ggplot() +
#   theme_classic() +
#   tidyterra::geom_spatraster(data = rf_pred_comp$CombConf, na.rm = TRUE) +
#   scale_fill_continuous(palette = "YlGn", na.value = "transparent") +
#   geom_contour(data = bathy_noaa, 
#                aes(x = x, y = y, z = z, fill = NULL), 
#                breaks = seq(from = -50, to = -5000, by = -250),
#                color = "darkgrey", 
#                linewidth = 0.3, 
#                alpha = 0.4) +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         axis.title = element_blank()) +
#   scale_x_continuous(limits = plot_xlim, expand = c(0,0)) +
#   scale_y_continuous(limits = plot_ylim, expand = c(0,0)) +
#   labs(title = paste("New | Period",poi,"| SSP",sspoi,"\n",subsample_absences,"|",keep_all_cmip_vars,"\n",vmeoi,"\nCombConf"))

# cowplot::plot_grid(plot_sdm2024_CombConf, plot_new_CombConf, 
#                    # labels = c("SDM2024 CombConf", "New CombConf"), 
#                    align = "hv",
#                    ncol = 2)
# ggsave("output/01_BlackCorals_RasterMetrics_OldVSNew/CombConf_comparison.png", 
#        width = 10, height = 5, dpi = 300)

# ## Number of models predicting presence ----
# # my_palette <- c("darkblue", paletteer::paletteer_d("colorBlindness::Blue2Orange10Steps", 10))

# plot_sdm2024_CVSum <- ggplot() +
#   theme_classic() +
#   tidyterra::geom_spatraster(data = sdm2024_pred_stack$CVSum, na.rm = TRUE) +
#   scale_fill_continuous(palette = "YlGnBu", na.value = "transparent",
#                         breaks = 0:10) +
#   geom_contour(data = bathy_noaa, 
#                aes(x = x, y = y, z = z, fill = NULL), 
#                breaks = seq(from = -50, to = -5000, by = -250),
#                color = "darkgrey", 
#                linewidth = 0.3, 
#                alpha = 0.4) +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         axis.title = element_blank()) +
#   scale_x_continuous(limits = plot_xlim, expand = c(0,0)) +
#   scale_y_continuous(limits = plot_ylim, expand = c(0,0)) +
#   labs(title = paste0("SDM2024\n",vmeoi,"\nCVSum"))

# plot_new_CVSum <- ggplot() +
#   theme_classic() +
#   tidyterra::geom_spatraster(data = rf_pred_comp$CVSum, na.rm = TRUE) +
#   scale_fill_continuous(palette = "YlGnBu", na.value = "transparent",
#                         breaks = 0:10) +
#   geom_contour(data = bathy_noaa, 
#                aes(x = x, y = y, z = z, fill = NULL), 
#                breaks = seq(from = -50, to = -5000, by = -250),
#                color = "darkgrey", 
#                linewidth = 0.3, 
#                alpha = 0.4) +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         axis.title = element_blank()) +
#   scale_x_continuous(limits = plot_xlim, expand = c(0,0)) +
#   scale_y_continuous(limits = plot_ylim, expand = c(0,0)) +
#   labs(title = paste("New | Period",poi,"| SSP",sspoi,"\n",subsample_absences,"|",keep_all_cmip_vars,"\n",vmeoi,"\nCVSum"))

# cowplot::plot_grid(plot_sdm2024_CVSum, plot_new_CVSum, 
#                    # labels = c("SDM2024 CVSum", "New CVSum"), 
#                    align = "hv",
#                    ncol = 2)
# ggsave(paste0("output/03_RF_Map_Outputs/",vmeoi,"_CVSum_comparison.png"), 
#        width = 10, height = 5, dpi = 300)
# # ggsave("output/01_BlackCorals_RasterMetrics_OldVSNew/CVSum_comparison.png", 
# #        width = 10, height = 5, dpi = 300)






# # transform_cmip_dep_to_raster <- function() {
  
# #   # Select a layer of data to raster
# #   df <- dep_df %>%
# #     sf::st_as_sf(coords = c("lon","lat"), crs = 4326)
  
# #   # Transform into terra vector of points
# #   pts <- terra::vect(df)
  
# #   # Determine resolution of data in degrees
# #   res <- round(ens_df$lon[2]-ens_df$lon[1], 5)
  
# #   # Create template raster
# #   rast_template <- terra::rast(
# #     xmin = sa_lims[1], xmax = sa_lims[2], ymin = sa_lims[3], ymax = sa_lims[4],
# #     resolution = res,
# #     crs = "EPSG:4326"
# #   )
  
# #   # Rasterise points to grid
# #   rast_result <- terra::rasterize(pts, rast_template, field = "dep")
  
# # }

# # dep_rast <- transform_cmip_dep_to_raster()
# # dep_rast_ext <- terra::ext(dep_rast)
# # res <- round(ens_df$lon[2]-ens_df$lon[1], 5)

# # dep_grid_plot <- ggplot() +
# #   theme_classic() +
# #   tidyterra::geom_spatraster(data = cmip_layers[[1]], aes(fill = last)) +
# #   cmocean::scale_fill_cmocean(name = "deep") +
# #   geom_sf(data = sa, aes(colour = "NAFO Study Area"), fill = NA, alpha = 0.8) +
# #   scale_colour_manual(name = "Boundary", 
# #                       values = c("NAFO Study Area" = "black")) +
# #   labs(x = "Longitude", y = "Latitude", fill = "Depth (m)") +
# #   # Add vertical and horizontal lines matching grid cell resolution
# #   coord_sf(xlim = c(dep_rast_ext[1], dep_rast_ext[2]), ylim = c(dep_rast_ext[3], dep_rast_ext[4]), expand = FALSE) +
# #   geom_vline(xintercept = seq(from = dep_rast_ext[1], to = dep_rast_ext[2], by = res),
# #              color = "black", size = 0.1, alpha = 0.2) +
# #   geom_hline(yintercept = seq(from = dep_rast_ext[3], to = dep_rast_ext[4], by = res),
# #              color = "black", size = 0.1, alpha = 0.2)


  