
# Explore data

# Load mapping data ----

# Load NAFO boundary
sa <- terra::rast("data/raw/Mapping_Layers/NRA_BNAM_b_tmp_mean.tif") %>%
  terra::as.polygons(.) %>%
  sf::st_as_sf(.) %>%
  sf::st_transform(4326) %>%
  mutate(NRA_BNAM_b_tmp_mean = 1) %>%
  group_by(NRA_BNAM_b_tmp_mean) %>%
  summarise(geometry = sf::st_union(geometry))

# Load NOAA bathymetry
bathy_noaa <- readRDS("data/raw/Mapping_Layers/bathy_noaa.rds")

# Transform response dataframe into sf
resp_sf <- sf::st_as_sf(resp_df, coords = c("Start_Long_DD", "Start_Lat_DD"), crs = 4326)

##
# Map the NAFO boundary and plot response data points ----
boundary_plot <- ggplot() +
  theme_classic() +
  geom_sf(data = sa, aes(colour = "NAFO Study Area"), fill = "lightblue", alpha = 0.8, show.legend = "line") +
  # Adjust colours
  scale_colour_manual(name = "Boundary", 
                      values = c("NAFO Study Area" = "black")) +
  geom_contour(data = bathy_noaa, 
               aes(x = x, y = y, z = z, fill = NULL), 
               breaks = seq(from = -50, to = -5000, by = -100),
               color = "darkgrey", 
               linewidth = 0.3, 
               alpha = 0.4) +
  # Add response points
  geom_sf(data = resp_sf, aes(fill = as.factor(VME_P_A)), 
          shape = 21, size = 2, alpha = 0.3) +
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
boundary_plot
