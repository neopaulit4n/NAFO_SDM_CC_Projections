# Trying to determine why spatial grid is offset from CMIP points

cross_icon <- leaflet::makeIcon(
  iconUrl = "data:image/svg+xml;charset=utf-8,%3Csvg xmlns='http://www.w3.org/2000/svg' width='12' height='12' viewBox='0 0 12 12'%3E%3Cline x1='0' y1='0' x2='12' y2='12' stroke='black' stroke-width='2'/%3E%3Cline x1='12' y1='0' x2='0' y2='12' stroke='black' stroke-width='2'/%3E%3C/svg%3E",
  iconWidth  = 12,
  iconHeight = 12,
  iconAnchorX = 6,  # Centre the icon on the point
  iconAnchorY = 6
)

rast_template_poly <- terra::as.polygons(rast_template)

leaflet::leaflet() %>%

  leaflet::addProviderTiles(leaflet::providers$OpenTopoMap) %>%
  
  leaflet::addRasterImage(
    x = rast_result,
    project = FALSE
  ) %>%  

  leaflet::addPolygons(data = rast_template_poly,
    fillColor = "transparent",
    fillOpacity = 0,
    color = "black",
    weight = 1,
    opacity = 0.7) %>%
  
  # Add study area boundary
  leaflet::addPolygons(data = sa,
    fillColor = "transparent",
    fillOpacity = 0,
    color = "red",
    weight = 1,
    opacity = 0.7) %>%

  # Add CMIP points
  leaflet::addMarkers(
    data = pts,
    icon = cross_icon,
    group = "CMIP points"
  ) %>%
  
  leaflet::addCircleMarkers(
    data = pts_template_ext,
    color = "red"
  )

# Note:
# Leaflet reprojects to 3857 and uses bilinear which causes some vertical distortion of the grid cells (purely visual).
# To remedy this, we can use the following:

# rast_3857 <- terra::project(rast_result, "EPSG:3857", method = "near")

# leaflet::leaflet() %>%
#   leaflet::addProviderTiles(leaflet::providers$OpenTopoMap) %>%
#   leaflet::addRasterImage(x = rast_3857, project = FALSE, method = "ngb") %>%
#   leaflet::addPolygons(data = rast_template_poly, fill = FALSE,
#                        color = "black", weight = 1, opacity = 0.7)