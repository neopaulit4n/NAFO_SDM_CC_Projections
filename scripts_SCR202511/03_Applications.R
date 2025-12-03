
library(tidyverse)
library(cmocean)

# Load data ----

# Load the ensemble data
ens_df <- readRDS("data/ens_df.rds") %>%
  pivot_wider(names_from = season, values_from = mldavg, names_glue = {"{.value}_{season}"}) %>%
  mutate(mldavg = coalesce(mldavg_W,mldavg_F,mldavg_Su,mldavg_Sp))

# Create mini df of vars, var abbrs, and var units for plots
var_df <- data.frame(var = colnames(select(ens_df, sobavg:mldavg)),
                     var_abbr = c("BS","SSS","BT","SST","BCS","BStr","MLD_W","MLD_Sp","MLD_Su","MLD_F","MLD"),
                     var_unit = c("","","(°C)","(°C)","(m/s)","(Pa)", rep("(m)", 5)))  # no units for salinity

# Load NAFO boundary
sa <- terra::rast("data/NRA_BNAM_b_tmp_mean.tif") %>%
  terra::as.polygons(.) %>%
  sf::st_as_sf(.) %>%
  sf::st_transform(4326) %>%
  mutate(NRA_BNAM_b_tmp_mean = 1) %>%
  group_by(NRA_BNAM_b_tmp_mean) %>%
  summarise(geometry = sf::st_union(geometry))

# 20-year time periods and barplots ----
ens_df_period_cell <- ens_df %>%
  mutate(ssp = as.factor(ssp)) %>%
  filter(!is.na(period)) %>%
  group_by(ssp, period, lon, lat) %>%
  # First summarise each var by lon/lat to get mean/min/max/sd var for each grid cell
  summarise(across(sobavg:mldavg, list(mean = ~mean(.x, na.rm = TRUE),
                                       min = ~min(.x, na.rm = TRUE),
                                       max = ~max(.x, na.rm = TRUE)),
                   .names = "{col}_{fn}")) %>%
  ungroup()
# Summarise mean of each stat over the entire sa
ens_df_period <- ens_df_period_cell %>%
  group_by(ssp, period) %>%
  summarise(across(sobavg_mean:mldavg_max, list(mean = ~mean(.x, na.rm = TRUE)),
                   .names = "{col}_{fn}")) %>%
  ungroup()

# Plot barplots
make_ssp_period_barplot <- function(var) {
  
  var_info <- filter(var_df, var == !!var)
  var_abbr <- var_info$var_abbr
  var_unit <- var_info$var_unit
  
  ggplot(ens_df_period, aes(x = ssp,
                            y = !!sym(paste0(var,"_mean_mean")))) +
    theme_bw() +
    geom_pointrange(aes(ymin = !!sym(paste0(var,"_min_mean")), 
                        ymax = !!sym(paste0(var,"_max_mean")),
                        colour = period),
                    position = position_dodge(0.9)) + #, size = 0.5, fatten = 2) +
    # geom_col(position = position_dodge()) +
    # # geom_bar(stat = "identity", position = position_dodge()) +
    # geom_errorbar(aes(ymin = !!sym(paste0(var,"_min_mean")), ymax = !!sym(paste0(var,"_max_mean"))),
    #               position = position_dodge(0.9), width = 0.25) +
    scale_colour_brewer(name = "Period", palette = "Dark2") +
    labs(x = "SSP",
         y = paste(var_abbr, var_unit)) 
  
  ggsave(paste0("output/05_Pointrange_SSP_Period/", var, "_bySSP_Period.jpg"), width = 6, height = 4)
  
}

lapply(var_df$var, make_ssp_period_barplot)


# Create map of mean/min/max for each var for SSP 126 P1/SSP 585 P4 and overlay sa ----
ssp126_p1 <- ens_df %>%
  filter(ssp == "1-2.6", period %in% "P1") %>%
  group_by(lon, lat) %>%
  summarise(across(mldavg:bstress, list(mean = ~mean(.x, na.rm = TRUE),
                                 min = ~min(.x, na.rm = TRUE),
                                 max = ~max(.x, na.rm = TRUE),
                                 range = ~max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE)),
                   .names = "{col}_{fn}")) %>%
  ungroup()

ssp585_p4 <- ens_df %>%
  filter(ssp == "5-8.5", period %in% "P4") %>%
  group_by(lon, lat) %>%
  summarise(across(mldavg:bstress, list(mean = ~mean(.x, na.rm = TRUE),
                                 min = ~min(.x, na.rm = TRUE),
                                 max = ~max(.x, na.rm = TRUE),
                                 range = ~max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE)),
                   .names = "{col}_{fn}")) %>%
  ungroup()


# Plot every variable-stat combination
library(cmocean)
ssp126_p1_plots <- list()
ssp585_p4_plots <- list()
varstats <- colnames(select(ssp126_p1, -lon, -lat))
var_abbr <- rep(c("MLD","BS","SSS","BT","SST","BCS","BStr"), each = 4)
var_units <- rep(c("(m)","","","(°C)","(°C)","(m/s)","(Pa)"), each = 4)  # no units for salinity

for (varstat in varstats) {
  varstat_idx <- which(varstats == varstat)
  var_abbr <- var_abbr[varstat_idx]
  var_unit <- var_units[varstat_idx]
  palette_name <- get_cmocean_palette(str_remove(varstat, "_.*"))
  
  ssp126_p1_plots[[varstat]] <- ggplot() +
    theme_bw() +
    geom_tile(data = ssp126_p1, aes(x = lon, y = lat, fill = !!sym(varstat))) +
    scale_fill_cmocean(paste(var_abbr, var_unit), name = palette_name,
                       limits = c(min(c(ssp126_p1[[varstat]],ssp585_p4[[varstat]])), 
                                  max(c(ssp126_p1[[varstat]],ssp585_p4[[varstat]])))) +
    # Add NAFO boundary
    geom_sf(data = sa, fill = NA, color = "black", linewidth = 0.8) +
    labs(title = paste(varstat, "SSP 1-2.6, P1"),
         x = "Longitude",
         y = "Latitude") +
    theme(legend.position = "none")
  
  ssp585_p4_plots[[varstat]] <- ggplot() +
    theme_bw() +
    geom_tile(data = ssp585_p4, aes(x = lon, y = lat, fill = !!sym(varstat))) +
    scale_fill_cmocean(paste(var_abbr, var_unit), name = palette_name,
                       limits = c(min(c(ssp126_p1[[varstat]],ssp585_p4[[varstat]])),
                                  max(c(ssp126_p1[[varstat]],ssp585_p4[[varstat]])))) +
    # Add NAFO boundary
    geom_sf(data = sa, fill = NA, color = "black", linewidth = 0.8) +
    labs(title = paste(varstat, "SSP 5-8.5, P4"),
         x = "Longitude",
         y = "Latitude") +
    theme(legend.position = "right",
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  # Save plots as jpg
  if (varstat_abbr == "BCS") {
    cowplot::plot_grid(ssp126_p1_plots[[varstat]], ssp585_p4_plots[[varstat]], 
                       ncol = 2, align = "h", rel_widths = c(0.9,1),
                       labels = c("A", "B"), label_x = c(0.14, 0.05), label_y = 0.92)    
  } else {
    cowplot::plot_grid(ssp126_p1_plots[[varstat]], ssp585_p4_plots[[varstat]], 
                       ncol = 2, align = "h", rel_widths = c(0.9325,1),
                       labels = c("A", "B"), label_x = c(0.14, 0.05), label_y = 0.92)        
  }
  
  ggsave(paste0("output/02_PairsSummStatsMaps/Period/", varstat, "_ssp126p1_ssp585p4.jpg"), width = 10, height = 5)
  
}

ssp126_p1_season <- ens_df_ssp_period_season %>%
  filter(ssp == "1-2.6", period %in% "P1") %>%
  select(lon, lat, contains("mldavg"), -contains("sd"))

ssp585_p4_season <- ens_df_ssp_period_season %>%
  filter(ssp == "5-8.5", period %in% "P4") %>%
  select(lon, lat, contains("mldavg"), -contains("sd"))

# Plot every variable-stat-season combination for MLD
ssp126_p1_season_plots <- list()
ssp585_p4_season_plots <- list()
varstats <- colnames(select(ssp126_p1_season, -lon, -lat))

for (varstat in varstats) {
  
  ssp126_p1_season_plots[[varstat]] <- ggplot() +
    theme_bw() +
    geom_tile(data = ssp126_p1_season, aes(x = lon, y = lat, fill = !!sym(varstat))) +
    scale_fill_cmocean("MLD (m)", name = "deep",
                       limits = c(min(c(ssp126_p1_season[[varstat]],ssp585_p4_season[[varstat]])), 
                                  max(c(ssp126_p1_season[[varstat]],ssp585_p4_season[[varstat]])))) +
    # Add NAFO boundary
    geom_sf(data = sa, fill = NA, color = "black", linewidth = 0.8) +
    labs(title = paste(varstat, "SSP 1-2.6, P1"),
         x = "Longitude",
         y = "Latitude") +
    theme(legend.position = "none")
  
  ssp585_p4_season_plots[[varstat]] <- ggplot() +
    theme_bw() +
    geom_tile(data = ssp585_p4_season, aes(x = lon, y = lat, fill = !!sym(varstat))) +
    scale_fill_cmocean("MLD (m)", name = "deep",
                       limits = c(min(c(ssp126_p1_season[[varstat]],ssp585_p4_season[[varstat]])), 
                                  max(c(ssp126_p1_season[[varstat]],ssp585_p4_season[[varstat]])))) +
    # Add NAFO boundary
    geom_sf(data = sa, fill = NA, color = "black", linewidth = 0.8) +
    labs(title = paste(varstat, "SSP 5-8.5, P4"),
         x = "Longitude",
         y = "Latitude") +
    theme(legend.position = "right",
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  cowplot::plot_grid(ssp126_p1_season_plots[[varstat]], ssp585_p4_season_plots[[varstat]], 
                     ncol = 2, align = "h", rel_widths = c(0.9325,1),
                     labels = c("A", "B"), label_x = c(0.14, 0.05), label_y = 0.92)
  
  ggsave(paste0("output/02_PairsSummStatsMaps/Season_MLD/", varstat, "_ssp126p1_ssp585p4.jpg"), width = 10, height = 5)
  
}



# Create pretty table of mean sst by ssp and period with +/- sd ----
# library(knitr)
# library(kableExtra)
# 
# 
# ens_tosavg_ssp_period_summary_df %>%
#   mutate(mean_sd = sprintf("%.2f ± %.2f", mean, sd)) %>%
#   select(SSP, period, mean_sd) %>%
#   pivot_wider(names_from = period, values_from = mean_sd) %>%
#   kable(caption = "Mean SST (tosavg) by SSP and time period (+/- standard deviation). P1: 2020-2039; P2: 2040-2059; P3: 2060-2079; P4: 2080-2099.") %>%
#   kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
#                 full_width = FALSE, position = "left") %>%
#   row_spec(0, bold = TRUE) %>%
#   column_spec(1, bold = TRUE)
make_ssp_period_table <- function (var) {
  temp <- ens_df %>%
    group_by(ssp, period) %>%
    summarise(mean = mean(!!sym(var), na.rm = TRUE),
              sd = sd(!!sym(var), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(mean_sd = sprintf("%.2f ± %.2f", mean, sd)) %>%
    select(ssp, period, mean_sd) %>%
    pivot_wider(names_from = period, values_from = mean_sd) %>%
    arrange(ssp) %>%
    rename(SSP = ssp,
           `2020-2039` = P1,
           `2040-2059` = P2,
           `2060-2079` = P3,
           `2080-2099` = P4) %>%
    select(-`NA`) %>%
    write_csv(paste0("output/06_Table_MeanSD_SSP_Period/Mean_",var,"_bySSP_Period.csv"))
  # return(temp)
}

lapply(var_df$var, make_ssp_period_table)


# Linegraph of annual var means for each SSP ----

make_annualmean_linegraph <- function (var) {
  var_info <- filter(var_df, var == !!var)
  var_abbr <- var_info$var_abbr
  var_unit <- var_info$var_unit
  
  ggplot(data = ens_df, aes(x = year, y = !!sym(var), colour = ssp, fill = ssp)) +
    theme_bw() +
    stat_summary(fun = mean, geom = "line") +
    # Add 95% confidence shading area
    stat_summary(fun.data = function(x) data.frame(y = mean(x, na.rm = TRUE),
                                                   ymin = mean(x, na.rm = TRUE) - 1.96 * sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)])),
                                                   ymax = mean(x, na.rm = TRUE) + 1.96 * sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))),
                 geom = "ribbon", alpha = 0.2, colour = NA) +
    # Add linear model to show trend
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.5) +
    # Apply manual colour scheme for SSPs
    scale_colour_manual("SSP",
                        values = c("1-2.6" = "#05B", 
                                   "2-4.5" = "darkgreen", 
                                   "3-7.0" = "darkorange", 
                                   "5-8.5" = "#000000")) +
    scale_fill_manual("SSP",
                      values = c("1-2.6" = "#05B", 
                                 "2-4.5" = "darkgreen", 
                                 "3-7.0" = "darkorange", 
                                 "5-8.5" = "#000000")) +
    labs(title = paste("Annual Mean", var, "by SSP"),
         x = "Year",
         y = paste(var_abbr, var_unit))
  
  ggsave(paste0("output/04_Linegraph_AnnualMean/", var, "_bySSP2.jpg"), width = 6, height = 4)
  cat(paste("Saved linegraph for", var, "\n"))
}

lapply(var_df$var, make_annualmean_linegraph)



  
# Old figures from list data ----
# ens_tosavg_ssp_annual_avg <- lapply(1:length(ssp), function(ssp_idx) {
#   ssp_name <- ssp[ssp_idx]
#   ssp_data <- ens_tosavg_ssp_date[[ssp_idx]]
#   
#   # Calculate annual means
#   annual_means <- lapply(1:85, function(year) {
#     start_idx <- (year - 1) * 12 + 1
#     end_idx <- year * 12
#     
#     # Extract 12 months for this year
#     year_data <- sapply(ssp_data[start_idx:end_idx], function(month) mean(month, na.rm = TRUE))
#     
#     # Calculate mean for the whole year
#     mean(year_data, na.rm = TRUE)
#   })
#   
#   names(annual_means) <- 2015:2099
#   return(data.frame(Year = 2015:2099, SST = unlist(annual_means), SSP = ssp_name))
# })
# ens_tosavg_ssp_annual_avg <- do.call(rbind, ens_tosavg_ssp_annual_avg)
# ens_tosavg_ssp_annual_avg$SSP <- as.factor(ens_tosavg_ssp_annual_avg$SSP)
# 
# ggplot(ens_tosavg_ssp_annual_avg, aes(x = Year, y = SST, colour = SSP)) +
#   theme_bw() +
#   geom_line() +
#   # Add dashed lm line for each SSP
#   geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.5) +
#   labs(title = "Annual Mean SST (tosavg) by SSP",
#        x = "Year",
#        y = "SST (°C)")



# Map of overall area and zoomed in plot of NAFO boundary ----

# Load Canadian EEZ boundary from shapefile
eez <- sf::st_read("data/eez/eez.shp") %>%
  sf::st_as_sf(crs = st_crs(4326)) %>%
  sf::st_make_valid() %>%
  sf::st_union() %>%
  sf::st_cast("POLYGON") %>%
  sf::st_cast("LINESTRING") %>%
  # Crop by sa_lims
  sf::st_crop(xmin = sa_lims[1]-0.09, xmax = sa_lims[2]+0.09, 
              ymin = sa_lims[3]-0.09, ymax = sa_lims[4]+0.09)

# Load small fishing footprint
footprint <- sf::st_read("data/FootprintProjectedShp/FootprintAreaProjected.shp") %>%
  sf::st_as_sf(crs = st_crs(4326)) %>%
  sf::st_make_valid() %>%
  sf::st_union() %>%
  sf::st_cast("LINESTRING")

boundary_plot <- ggplot() +
  theme_classic() +
  geom_sf(data = sa, aes(colour = "NAFO Study Area"), fill = "lightblue", alpha = 0.8, show.legend = "line") +
  geom_sf(data = footprint, aes(colour = "NAFO Fishing Footprint"), show.legend = "line") +  
  geom_sf(data = eez, aes(colour = "Canadian EEZ"), inherit.aes = FALSE, linewidth = 1, show.legend = "line") +
  # Adjust colours
  scale_colour_manual(name = "Boundary", 
                      values = c("NAFO Study Area" = "black", "NAFO Fishing Footprint" = "blue", "Canadian EEZ" = "red")) +
  geom_contour(data = bathy_noaa, 
               aes(x = x, y = y, z = z, fill = NULL), 
               breaks = seq(from = -50, to = -5000, by = -100),
               color = "darkgrey", 
               linewidth = 0.3, 
               alpha = 0.4) +
  # Add text to bathymetry contours
  metR::geom_text_contour(data = bathy_noaa, 
                    aes(x = x, y = y, z = z), 
                    breaks = c(-50, seq(from =-200, to = -2000, by = -200), 
                               seq(from = -2500, to = -5000, by = -500)),
                    color = "darkgrey", 
                    size = 3,
                    stroke = 0.1,
                    alpha = 0.6) +
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(limits = c(min(lon)-0.09,max(lon)+0.09), expand = c(0, 0)) +
  scale_y_continuous(limits = c(min(lat)-0.09,max(lat)+0.09), expand = c(0, 0))

  # Add top and right axis lines and ticks
  # theme(axis.line.x.top = element_line(color = "black"),
  #       axis.ticks.x.top = element_line(color = "black"),
  #       axis.line.y.right = element_line(color = "black"),
  #       axis.ticks.y.right = element_line(color = "black"),
  #       axis.ticks.length = unit(0.15, "cm"),
  #       axis.title = element_text(size = 14),
  #       axis.text = element_text(size = 12))
  
  # Add north arrow and scale bar
  # ggspatial::annotation_scale(location = "tl", width_hint = 0.3) +
  # ggspatial::annotation_north_arrow(location = "tl", which_north = "true", 
  #                        pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
  #                        style = "north_arrow_nautical") +
  # theme(legend.position = "none")
ggsave("output/03_SA_Map/Map_SA_bathy.jpg", width = 6, height = 6)



# Create gridded depth plots with NAFO SA boundary ----
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
  tidyterra::geom_spatraster(data = dep_rast, aes(fill = last)) +
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

cowplot::plot_grid(boundary_plot + 
                     theme(legend.position = "bottom", legend.direction = "horizontal",
                           legend.text.position = "bottom",
                           legend.title.position = "top",
                           legend.title = element_text(size = 10, hjust = 0.5),
                           legend.text = element_text(size = 8)), 
                   dep_grid_plot +
                     theme(legend.position = "bottom", legend.direction = "horizontal",
                           legend.text.position = "bottom",
                           legend.title.position = "top",
                           legend.title = element_text(size = 10),
                           legend.text = element_text(size = 6)),
                   ncol = 2, labels = c("A","B"), rel_widths = c(1,1), align = "h")

ggsave("output/03_SA_Map/DualMap_NAFOBoundary_DepGrid.jpg", width = 10, height = 5)


# Create a tiff for every ssp/period and varstat ----
varstats <- colnames(select(ens_df_period_cell, sobavg_mean:mldavg_max))
periods <- unique(ens_df_period_cell$period)
ssps <- unique(ens_df_period_cell$ssp)

lapply(periods, function(period) {
  lapply(ssps, function(ssp) {
    lapply(varstats, function(varstat) {
      
      # Select a layer of data to raster
      df <- ens_df_period_cell %>%
        filter(ssp == ssp, period == period) %>%
        select(lon, lat, !!sym(varstat)) %>%
        sf::st_as_sf(coords = c("lon","lat"), crs = 4326)
      
      # Transform into terra vector of points
      pts <- terra::vect(df)
      
      # Determine resolution of data in degrees
      res <- round(ens_df_period_cell$lon[2]-ens_df_period_cell$lon[1], 5)
      
      # Create template raster
      rast_template <- terra::rast(
        xmin = sa_lims[1], xmax = sa_lims[2], ymin = sa_lims[3], ymax = sa_lims[4],
        resolution = res,
        crs = "EPSG:4326"
      )
      
      # Rasterise points to grid
      rast_result <- terra::rasterize(pts, rast_template, field = varstat)
      
      # # Convert layer_df into a raster layer
      # r <- terra::rast(layer_df, type = "xyz", crs = "EPSG:4326")
      
      # Save raster as tiff
      terra::writeRaster(rast_result, paste0("data/Varstat_SSP_Period_tiff/", varstat, "_", ssp, "_", period, ".tiff"), 
                         overwrite = TRUE)
      
      cat(paste("Saved tiff for", varstat, ssp, period, "\n"))
      
    })
  })
})


