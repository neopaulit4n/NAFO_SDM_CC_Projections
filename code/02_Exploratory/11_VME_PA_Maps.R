# Mapping all layers with leaflet
library(tidyverse)
# library(leaflet)

# Load all desired layers ----

# NAFO study area
sa <- terra::rast("data/raw/BNAM_Data_From_Cam/BNAM_From_NAFO_SharePoint/NRA_BNAM_b_cur_avg_max.tif") %>%
  terra::as.polygons(.) %>%
  sf::st_as_sf(.) %>%
  sf::st_transform(4326) %>%
  mutate(NRA_BNAM_b_tmp_mean = 1) %>%
  group_by(NRA_BNAM_b_tmp_mean) %>%
  summarise(geometry = sf::st_union(geometry))

# Bathymetry
bathy_noaa <- readRDS("data/raw/Mapping_Layers/bathy_noaa.rds")

# KDE polygons
kde_poly <- sf::read_sf("data/raw/Mapping_Layers/NAFO_2025_VME_Goup_Threshold_KDE_Polygons/KDE_Analyses_VME_2025_Threshold_Polygons_BlackCoral.shp") %>%
    sf::st_transform(4326)

# VME closures
vme_closures <- sf::read_sf("data/raw/Mapping_Layers/NAFO_VME_closures_2022/NAFO_VME_closures_2022.shp")

# VME presence/absence points
pa <- read_csv("data/cleaned/VME_group_PA_df.csv", show_col_types = FALSE) %>%
  filter(VME_Group == vmeoi) %>%
  sf::st_as_sf(coords = c("Start_Long_DD", "Start_Lat_DD")) %>%
  sf::st_set_crs(4326)
pres <- filter(pa, VME_P_A == 1)
abs <- filter(pa, VME_P_A == 0)

ggplot() +
  theme_classic() +    
  geom_contour_filled(
    data = bathy_noaa, 
    aes(x = x, y = y, z = z), 
    breaks = c(0,-100,-200,-500,-1000,-2000,-2500,-3000,-4000,-10000),
    color = "black", 
    linewidth = 0.1, 
    alpha = 0.4) +
  scale_fill_discrete(palette = "Greys") +
  geom_sf(data = sa, colour = "black", fill = "transparent", linewidth = 1.1) +
  geom_sf(data = abs, colour = "black", size = 0.5) +
  geom_sf(data = pres, colour = "black", fill = "red", shape = 21, size  = 3) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title = element_blank())
ggsave(paste0(output_folder,"/",vmeoi,"_PA_Map.jpg"), dpi = 300, width = 10, height = 10)
