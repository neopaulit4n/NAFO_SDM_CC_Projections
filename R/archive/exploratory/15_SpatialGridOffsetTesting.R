# Trying to determine why spatial grid is offset from CMIP points

cross_icon <- makeIcon(
  iconUrl = "data:image/svg+xml;charset=utf-8,%3Csvg xmlns='http://www.w3.org/2000/svg' width='12' height='12' viewBox='0 0 12 12'%3E%3Cline x1='0' y1='0' x2='12' y2='12' stroke='black' stroke-width='2'/%3E%3Cline x1='12' y1='0' x2='0' y2='12' stroke='black' stroke-width='2'/%3E%3C/svg%3E",
  iconWidth  = 12,
  iconHeight = 12,
  iconAnchorX = 6,  # Centre the icon on the point
  iconAnchorY = 6
)

rast_template_poly <- terra::as.polygons(rast_template)
rast_result_poly <- terra::as.polygons(rast_result)

leaflet::leaflet() %>%

  leaflet::addProviderTiles(providers$OpenTopoMap) %>%
  
  leaflet::addRasterImage(
    x = rast_result) %>%  # cmip_layers[[1]]

  # addPolygons(data = rast_template_poly,
  #   fillColor = "transparent",
  #   fillOpacity = 0,
  #   color = "black",
  #   weight = 1,
  #   opacity = 0.7) %>%
  
  addPolygons(data = rast_result_poly,
    fillColor = "transparent",
    fillOpacity = 0,
    color = "blue",
    weight = 1,
    opacity = 0.7) %>%
  # leaflet::addRasterImage(
  #   x = bathy_layers[[1]],  # cmip_layers[[1]]
  # ) %>%
  # Add CMIP points
  leaflet::addMarkers(
    data = pts,
    icon = cross_icon,
    group = "CMIP points"
  ) #%>%
  # leaflet::addCircleMarkers(
  #   data = sa_lims,
  #   color = "red"
  # )
