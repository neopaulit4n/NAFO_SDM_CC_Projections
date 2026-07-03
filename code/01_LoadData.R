
# Editing and transferring SDM2024 code for use with CC projection model data

library(tidyverse)

# Load response ----
cat("Loading response dataframe\n")
resp_df <- read_csv("data/processed/VME_group_PA_df.csv", show_col_types = FALSE)

# Load predictors ----

## Load terrain variables (static) ----
cat("Loading terrain variables\n")
bathy_layers <- list.files(path = "data/raw/Bathy_Layers", 
                           pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(., nm = basename(.) %>% tools::file_path_sans_ext()) %>%
  lapply(terra::rast) %>%
  lapply(terra::project, "EPSG:4326")

names(bathy_layers) <- gsub("GEBCO2024_FS005_StudyArea_","",names(bathy_layers))
names(bathy_layers)[1] <- "FS005"

## Load ensembled CMIP data ----
cat("Loading CMIP data\n")
cmip_df <- readRDS("data/processed/cmip_ens_proj_df.rds") %>%
  rename(
    BS = sobavg,
    SSS = sosavg,
    BT = tobavg,
    SST = tosavg,
    BCS = wobavg,
    BStr = bstress,
    MLD = mldavg,
    MLD_Su = mldavg_Su,
    MLD_F = mldavg_F,
    MLD_W = mldavg_W,
    MLD_Sp = mldavg_Sp
  )
  
# Load selected terrain static variables ----
terrain_topvars <- read_csv("data/processed/VarImp_2024_2025_SCR_02_TopStaticVarsByVME.csv", show_col_types = FALSE)

# Create baseline dataframe that will build the models used for predictions ----
baseline_df <- readRDS("data/processed/cmip_ens_1993_2014_df.rds") %>%
  rename(
    BS = sobavg,
    SSS = sosavg,
    BT = tobavg,
    SST = tosavg,
    BCS = wobavg,
    BStr = bstress,
    MLD = mldavg,
    MLD_Su = mldavg_Su,
    MLD_F = mldavg_F,
    MLD_W = mldavg_W,
    MLD_Sp = mldavg_Sp
  ) %>%
  summarise(
    across(SST:BStr, 
      list(
        mean = ~mean(.x, na.rm = TRUE),
        min = ~min(.x, na.rm = TRUE),
        max = ~max(.x, na.rm = TRUE), 
        range = ~max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE)
      ),
      .names = "{col}_{fn}"),
    .by = c(lon, lat))

# Average predictors for each period and SSP (will use these to predict on) ----
cmip_df_period_ssp <- cmip_df %>%
  filter(!is.na(period)) %>%
  summarise(
    across(BS:BStr, 
      list(
        mean = ~mean(.x, na.rm = TRUE),
        min = ~min(.x, na.rm = TRUE),
        max = ~max(.x, na.rm = TRUE), 
        range = ~max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE)
      ),
      .names = "{col}_{fn}"),
    .by = c(lon, lat, period, ssp))

# All combinations ----
vme_all <- unique(resp_df$VME_Group)
period_all <- c("P1", "P2", "P3", "P4")
ssp_all <- unique(cmip_df_period_ssp$ssp)

# Get study area extent and create spatial mask ----
sa <- terra::rast("data/raw/Bathy_Layers/GEBCO2024_FS005.tif") %>%
  terra::as.polygons(.) %>%
  sf::st_as_sf(.) %>%
  sf::st_transform(4326) %>%
  summarise(geometry = sf::st_union(geometry))

sa_lims <- sf::st_coordinates(sa) %>% as.data.frame %>%
  select(X,Y)
sa_lims <- c(min(sa_lims$X), max(sa_lims$X),min(sa_lims$Y),max(sa_lims$Y))

# Transform CMIP data to raster ----
cat("Transforming CMIP dataframe to raster layers\n")
transform_cmip_to_raster <- function(data, varstat, poi = NULL, sspoi = NULL) {

  # Select a layer of data to raster
  df <- data %>%
    filter(if(!is.null(sspoi)) ssp == sspoi else TRUE) %>%
    filter(if(!is.null(poi)) period == poi else TRUE) %>%
    select(lon, lat, !!sym(varstat)) %>%
    sf::st_as_sf(coords = c("lon","lat"), crs = 4326)
  
  # Transform into terra vector of points
  pts <- terra::vect(df)
  
  # Determine resolution of data in degrees
  res <- round(sort(unique(data$lon))[2] - sort(unique(data$lon))[1], 5)
  
  # Create template raster
  rast_template <- terra::rast(
    xmin = sa_lims[1], xmax = sa_lims[2], ymin = sa_lims[3], ymax = sa_lims[4],
    resolution = res,
    crs = "EPSG:4326"
  )
  
  # Rasterise points to grid
  rast_result <- terra::rasterize(pts, rast_template, field = varstat)
  
  # Resample to match terrain rasters
  rast_resamp <- terra::resample(rast_result, bathy_layers[[1]], method = "bilinear")

  # Mask by SA
  rast_crop <- terra::mask(rast_resamp, bathy_layers[[1]])

  return(rast_crop)
}

cmip_vars <- colnames(select(cmip_df_period_ssp, contains("_")))

cmip_layers <- lapply(cmip_vars, function(var) {
  transform_cmip_to_raster(data = baseline_df, varstat = var)
}) %>%
  set_names(cmip_vars)

# Extract data from raster layers to the VME response points ----
cat("Extracting data from raster layers to VME PA points\n")
suppressMessages(cmip_pred_df <- lapply(c(bathy_layers, cmip_layers), function(layer) {
  terra::extract(layer, select(resp_df, Start_Long_DD, Start_Lat_DD)) %>%
  select(-ID)
}) %>%
  bind_cols() %>%
  set_names(c(names(bathy_layers), names(cmip_layers)))
)

# Combine predictor and response dataframes ----
cmip_comb_df <- bind_cols(resp_df, cmip_pred_df) %>%
  mutate(VME_P_A = factor(VME_P_A, levels = c(0, 1), labels = c("Absence", "Presence"))) %>%
  drop_na()

# Create CMIP raster layers for projected scenarios ----
cat("Creating projected CMIP raster layers\n")
cmip_layers_proj <- lapply(period_all, function(poi) {
  lapply(ssp_all, function(sspoi) {
    lapply(cmip_vars, function(var) {
      transform_cmip_to_raster(data = cmip_df_period_ssp, sspoi = sspoi, poi = poi, varstat = var)
    }) %>%
      set_names(cmip_vars)
  }) %>%
    set_names(ssp_all)  
}) %>%
  set_names(period_all)

cat("Done!\n")