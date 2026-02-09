
# Load mapping layers ----

# Load NOAA bathymetry layer for contours
bathy_noaa <- readRDS("data/raw/Mapping_Layers/bathy_noaa.rds")

# Load Canadian EEZ boundary from shapefile
eez <- sf::st_read("data/raw/Mapping_Layers/eez/eez.shp") %>%
  sf::st_as_sf(crs = st_crs(4326)) %>%
  sf::st_make_valid() %>%
  sf::st_union() %>%
  sf::st_cast("POLYGON") %>%
  sf::st_cast("LINESTRING") %>%
  # Crop by sa_lims
  sf::st_crop(xmin = sa_lims[1]-0.09, xmax = sa_lims[2]+0.09, 
              ymin = sa_lims[3]-0.09, ymax = sa_lims[4]+0.09)

# Load small fishing footprint
footprint <- sf::st_read("data/raw/Mapping_Layers/FootprintProjectedShp/FootprintAreaProjected.shp") %>%
  sf::st_as_sf(crs = st_crs(4326)) %>%
  sf::st_make_valid() %>%
  sf::st_union() %>%
  sf::st_cast("LINESTRING")


# Map of overall area and zoomed in plot of NAFO boundary ----
ggplot() +
  theme_classic() +
  # geom_sf(data = sa, aes(colour = "NAFO Study Area"), fill = "lightblue", alpha = 0.8, show.legend = "line") +
  tidyterra::geom_spatraster(data = rf_pred_pa, na.rm = TRUE) +
  scale_fill_manual(values = c("0" = "#ffebcd", "1" = "#b87333"),
                    na.value = "transparent",
                    na.translate = FALSE,  # remove NAs from legend
                    labels = c("0" = "Absence", "1" = "Presence")) +
  geom_sf(data = footprint, aes(colour = "NAFO Fishing Footprint"), show.legend = "line") +  
  geom_sf(data = eez, aes(colour = "Canadian EEZ"), inherit.aes = FALSE, linewidth = 1, show.legend = "line") +
  # Adjust colours
  scale_colour_manual(name = "Boundary", 
                      values = c("NAFO Study Area" = "black", "NAFO Fishing Footprint" = "blue", "Canadian EEZ" = "red")) +
  geom_contour(data = bathy_noaa, 
               aes(x = x, y = y, z = z, fill = NULL), 
               breaks = seq(from = -50, to = -5000, by = -250),
               color = "darkgrey", 
               linewidth = 0.3, 
               alpha = 0.4) +
  coord_sf(xlim = terra::ext(rf_pred_pa)[1:2], 
           ylim = terra::ext(rf_pred_pa)[3:4], expand = FALSE) +
  # geom_vline(xintercept = seq(from = terra::ext(rf_pred_pa)[1], 
  #                             to = terra::ext(rf_pred_pa)[2], 
  #                             by = terra::res(rf_pred_pa)[1]),
  #            color = "black", linewidth = 0.1, alpha = 0.1) +
  # geom_hline(yintercept = seq(from = terra::ext(rf_pred_pa)[3], 
  #                             to = terra::ext(rf_pred_pa)[4], 
  #                             by = terra::res(rf_pred_pa)[2]),
  #            color = "black", linewidth = 0.1, alpha = 0.1) +
  labs(title = paste("Predicted Presence/Absence for", vme_group),
       fill = "Prediction", x = "Longitude", y = "Latitude")

transform_cmip_dep_to_raster <- function() {
  
  # Select a layer of data to raster
  df <- dep_df %>%
    sf::st_as_sf(coords = c("lon","lat"), crs = 4326)
  
  # Transform into terra vector of points
  pts <- terra::vect(df)
  
  # Determine resolution of data in degrees
  res <- round(ens_df$lon[2]-ens_df$lon[1], 5)
  
  # Create template raster
  rast_template <- terra::rast(
    xmin = sa_lims[1], xmax = sa_lims[2], ymin = sa_lims[3], ymax = sa_lims[4],
    resolution = res,
    crs = "EPSG:4326"
  )
  
  # Rasterise points to grid
  rast_result <- terra::rasterize(pts, rast_template, field = "dep")
  
}

dep_rast <- transform_cmip_dep_to_raster()
dep_rast_ext <- terra::ext(dep_rast)
res <- round(ens_df$lon[2]-ens_df$lon[1], 5)

dep_grid_plot <- ggplot() +
  theme_classic() +
  tidyterra::geom_spatraster(data = cmip_layers[[1]], aes(fill = last)) +
  cmocean::scale_fill_cmocean(name = "deep") +
  geom_sf(data = sa, aes(colour = "NAFO Study Area"), fill = NA, alpha = 0.8) +
  scale_colour_manual(name = "Boundary", 
                      values = c("NAFO Study Area" = "black")) +
  labs(x = "Longitude", y = "Latitude", fill = "Depth (m)") +
  # Add vertical and horizontal lines matching grid cell resolution
  coord_sf(xlim = c(dep_rast_ext[1], dep_rast_ext[2]), ylim = c(dep_rast_ext[3], dep_rast_ext[4]), expand = FALSE) +
  geom_vline(xintercept = seq(from = dep_rast_ext[1], to = dep_rast_ext[2], by = res),
             color = "black", size = 0.1, alpha = 0.2) +
  geom_hline(yintercept = seq(from = dep_rast_ext[3], to = dep_rast_ext[4], by = res),
             color = "black", size = 0.1, alpha = 0.2)


