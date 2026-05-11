# Mapping all layers with leaflet
library(tidyverse)
library(leaflet)

# Load all desired layers ----

# NAFO study area
sa <- terra::rast("data/raw/BNAM_Data_From_Cam/BNAM_From_NAFO_SharePoint/NRA_BNAM_b_cur_avg_max.tif") %>%
  terra::as.polygons(.) %>%
  sf::st_as_sf(.) %>%
  sf::st_transform(4326) %>%
  mutate(NRA_BNAM_b_tmp_mean = 1) %>%
  group_by(NRA_BNAM_b_tmp_mean) %>%
  summarise(geometry = sf::st_union(geometry))

# KDE polygons
kde_poly <- sf::read_sf("data/raw/Mapping_Layers/NAFO_2025_VME_Goup_Threshold_KDE_Polygons/KDE_Analyses_VME_2025_Threshold_Polygons_BlackCoral.shp") %>%
    sf::st_transform(4326)

# VME closures
vme_closures <- sf::read_sf("data/raw/Mapping_Layers/NAFO_VME_closures_2022/NAFO_VME_closures_2022.shp")

# VME presence/absence points
pa <- read_csv("data/cleaned/VME_group_PA_df.csv", show_col_types = FALSE) %>%
  filter(VME_Group == vmeoi) %>%
  sf::st_as_sf(coords = c("Start_Long_DD", "Start_Lat_DD"))
pres <- filter(pa, VME_P_A == 1)
abs <- filter(pa, VME_P_A == 0)
# cross_icon <- makeIcon(
#   iconUrl = "data:image/svg+xml;charset=utf-8,%3Csvg xmlns='http://www.w3.org/2000/svg' width='12' height='12' viewBox='0 0 12 12'%3E%3Cline x1='0' y1='0' x2='12' y2='12' stroke='black' stroke-width='2'/%3E%3Cline x1='12' y1='0' x2='0' y2='12' stroke='black' stroke-width='2'/%3E%3C/svg%3E",
#   iconWidth  = 12,
#   iconHeight = 12,
#   iconAnchorX = 6,  # Centre the icon on the point
#   iconAnchorY = 6
# )

# Bathymetry/terrain layers
bathy_layers <- list.files(path = "data/raw/BNAM_Data_From_Cam/Bathymetry_Terrain_From_NAFO_SharePoint", 
                           pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(., nm = basename(.) %>% tools::file_path_sans_ext()) %>%
  lapply(terra::rast) %>%
  lapply(terra::project, "EPSG:4326")

names(bathy_layers) <- gsub("GEBCO2024_FS005_StudyArea_","",names(bathy_layers))
names(bathy_layers)[1] <- "FS005"

# CMIP layers -> baseline and by SSP/period
baseline_layers <- as.list(cmip_layers) %>%
  set_names(names(cmip_layers))

projection_layers <- unlist(map(unlist(cmip_layers_future), function(i) {
  as.list(i) %>%
    set_names(cmip_vars)
}))
projection_layers <- unlist(cmip_layers_future)

all_layers <- c(
  bathy_layers, 
  baseline_layers,
  projection_layers
)

leaflet_map <- reduce(
  names(all_layers),
  function(map, layer_name) {
    addRasterImage(
      map,
      x = all_layers[[layer_name]],
      # colors = bathy_pal,
      group = layer_name,
      opacity = 0.8
    )
  },
  .init = leaflet() |> 
    addTiles() |> 
    addLayersControl(
      baseGroups = names(all_layers), 
      options = layersControlOptions(collapsed = FALSE))
) 


leaflet() %>%
  # fitBounds(-63.2,43.5,-61.8,46) %>%
  
  # Basemap
  addProviderTiles(providers$OpenTopoMap) %>%
  
  addMapPane('baselayer', zIndex = 410) %>%
  addMapPane("NAFO_SA", zIndex = 420) %>%  
  addMapPane('absence', zIndex = 430) %>%  
  addMapPane("presence", zIndex = 435) %>%
  addMapPane("polygons", zIndex = 440) %>%
  
  # NAFO SA boundary
  addPolylines(
    data = sa,
    color = "black",
    group = "NAFO SA",
    options = pathOptions(pane = 'NAFO_SA')
  ) %>%

  # VME closures
  addPolygons(
    data = vme_closures,
    popup = "VME closures",
    opacity = 1,
    weight = 1,
    color = 'red',
    fillOpacity = 0.5,
    group = 'VME closures',
    options = pathOptions(pane = 'polygons')
  ) %>%  
  
  # KDE polygons
  addPolygons(
    data = kde_poly,
    popup = "KDE Polygons",
    # popup = ~htmlEscape(GEOLOGY),
    # fillColor = ~pal_geo(GEOLOGY),
    opacity = 1,
    weight = 1,
    color = 'blue',
    fillOpacity = 0.5,
    group = 'KDE polygons',
    options = pathOptions(pane = 'polygons')
  ) %>%
    
  # VME resence/absence points
  # addCircleMarkers(
  #   data = pa,
  #   # popup = ~htmlEscape(filename),
  #   # fillColor = ~pal_subst(VME_P_A),
  #   fillColor = ~case_match(VME_P_A, 0 ~ "grey", 1 ~ "green"),
  #   opacity = 0.3,
  #   weight = 1,
  #   color = 'black',
  #   fillOpacity = 0.8,
  #   group = 'PA'
  # )

  addMarkers(
    data = abs,
    icon = cross_icon,
    group = "VME Absences",
    options = pathOptions(pane = 'absence')
  ) %>%

  addCircleMarkers(
    data = pres,
    fillColor = "green",
    opacity = 0.3,
    weight = 1,
    color = 'black',
    fillOpacity = 0.8,
    group = 'VME Presences',
    options = pathOptions(pane = 'presence')
  ) %>%
  
  # Layers control
  addLayersControl(
    overlayGroups = c("VME closures","KDE polygons","VME Absences","VME Presences"),
    options = layersControlOptions(collapsed = FALSE)
  ) 


  

      
  # Bathymetry raster
  # addRasterImage(bathy,
  #                maxBytes = 27000000,
  #                colors = pal_bathy,
  #                opacity = 0.8,
  #                group = 'Bathymetry') %>%
    

  


    