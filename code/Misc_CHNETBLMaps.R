# Ellen wanted to see maps from CHNETBL3/5 and compare to DEM

dem <- terra::rast("data/raw/BNAM_Data_From_Cam/Bathymetry_Terrain_From_NAFO_SharePoint/GEBCO2024_FS005.tif")
chnetbl3 <- terra::rast("data/raw/BNAM_Data_From_Cam/Bathymetry_Terrain_From_NAFO_SharePoint/GEBCO2024_FS005_StudyArea_CHNETBL3.tif")
chnetbl5 <- terra::rast("data/raw/BNAM_Data_From_Cam/Bathymetry_Terrain_From_NAFO_SharePoint/GEBCO2024_FS005_StudyArea_CHNETBL5.tif")
chnetd3 <- terra::rast("data/raw/BNAM_Data_From_Cam/Bathymetry_Terrain_From_NAFO_SharePoint/GEBCO2024_FS005_StudyArea_CHNETD3.tif")
chnetd5 <- terra::rast("data/raw/BNAM_Data_From_Cam/Bathymetry_Terrain_From_NAFO_SharePoint/GEBCO2024_FS005_StudyArea_CHNETD5.tif")
lyrs <- c(dem, chnetbl3, chnetbl5, chnetd3, chnetd5)
names(lyrs) <- c("DEM","CHNETBL3","CHNETBL5","CHNETD3","CHNETD5")

ggplot() +
  tidyterra::geom_spatraster(data = lyrs) +
  scale_fill_viridis_c(na.value = "transparent") +
  facet_wrap(~ lyr) +
  theme_classic()
