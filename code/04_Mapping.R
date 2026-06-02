# Load mapping layers ----

if (!dir.exists(paste0(output_folder,"/ModellingMaps"))) dir.create(paste0(output_folder,"/ModellingMaps"))

# Load NOAA bathymetry layer for contours
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
metric_names <- c("MaxClass", "MaxClassF", "MaxClassAvgProb", "CombConf", "CVSum", "rawPresenceProb", "rawAbsenceProb")

## Read in rasters ----
rf_pred_proj_all <- lapply(metric_names, function(metric) {
  metric_pred_names <- list.files(paste0(output_folder, "/RFModelRasters"), pattern = paste0("rf_res_proj_", metric, "_"), full.names = TRUE)
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

rf_pred_baseline_all <- lapply(metric_names, function(metric) {
  metric_pred_names <- list.files(paste0(output_folder, "/RFModelRasters"), pattern = paste0("rf_res_baseline_", metric, "\\.tif"), full.names = TRUE)
  metric_preds <- terra::rast(metric_pred_names)
  
  # Factorise MaxClass rasters for plotting
  if (metric == "MaxClass") {
    metric_preds <- terra::as.factor(metric_preds)
  }

  return(metric_preds)
}) %>%
  set_names(metric_names)

## Calculate area of predicted presence (MaxClass) ----
# compute_presence_areas <- function(metric) {
#   rast_stack <- rf_pred_all[[metric]]
  
#   layer_names <- names(rast_stack)
  
#   areas <- sapply(layer_names, function(lyr_name) {
#     r_layer <- rast_stack[[lyr_name]]
#     # Subset raster to presence cells only (mask out absence)
#     presence <- terra::ifel(r_layer == 1, r_layer, NA)
#     # Calculate area of each cell in km²
#     area_rast <- terra::cellSize(presence, unit = "km", mask = TRUE)
#     # Sum all presence cell areas to get total area in km²
#     terra::global(area_rast, "sum", na.rm = TRUE)$sum
#   })
  
#   data.frame(
#     lyr = layer_names,
#     label = paste0(format(round(areas, 0), big.mark = ","), " km²")
#   )
# }


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
    "MaxClassAvgProb" = function() {scale_fill_gradient2(low = "#1164b4", mid = "#ffff99", high = "#e03c31", midpoint = 0.5, na.value = "transparent", limits = c(0, 1))},
    "CombConf" = function() {scale_fill_continuous(palette = "YlGn", na.value = "transparent", limits = c(0, 1))},
    "CVSum" = function() {scale_fill_continuous(palette = "YlGn", na.value = "transparent", limits = c(0, 10))},
    "rawPresenceProb" = function() {scale_fill_continuous(palette = "YlGn", na.value = "transparent", limits = c(0, 1))},
    "rawAbsenceProb" = function() {scale_fill_continuous(palette = "YlOrRd", na.value = "transparent", limits = c(0, 1))}
  )

  # Create plot
  p <- ggplot() +
    theme_classic() +    
    tidyterra::geom_spatraster(data = rf_pred_proj_all[[metric]], na.rm = TRUE) +
    facet_wrap(~ lyr, ncol = 4) +
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

  ggsave(paste0(output_folder,"/ModellingMaps/",vmeoi,"_proj_",metric,"_facet.jpg"), p,
    width = 10, height = 10, dpi = 300)

  p_baseline <- ggplot() +
    theme_classic() +    
    tidyterra::geom_spatraster(data = rf_pred_baseline_all[[metric]], na.rm = TRUE) +
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

  ggsave(paste0(output_folder,"/ModellingMaps/",vmeoi,"_baseline_",metric,".jpg"), p_baseline,
    width = 10, height = 10, dpi = 300)
})

  