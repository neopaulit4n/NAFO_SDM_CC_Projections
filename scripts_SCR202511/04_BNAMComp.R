
library(terra)
library(cmocean)

# Compare CMIP layers to BNAM output

# BNAM is just the one time frame 1990-2023  so you should compare with the 
# P1 for CMIP and I think the S2-4.5 SSP. I thought a surface map showing 
# differences in case there is a spatial bias and then maybe a correlation.

# List all TIF files, read them into a list, and convert to a long-form dataframe
bnam_layers <- list.files("data/BNAM_Data_From_Cam/BNAM_From_NAFO_SharePoint", pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(., nm = basename(.) %>% tools::file_path_sans_ext()) %>%
  lapply(terra::rast) %>%
  lapply(terra::project, "EPSG:4326")

# Transform CMIP data to raster
transform_cmip_to_raster <- function(varstat, period, ssp) {
  
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
  rast_template <- rast(
    xmin = sa_lims[1], xmax = sa_lims[2], ymin = sa_lims[3], ymax = sa_lims[4],
    resolution = res,
    crs = "EPSG:4326"
  )
  
  # Rasterise points to grid
  rast_result <- rasterize(pts, rast_template, field = varstat)
  
  # Resample to match BNAM raster
  raster_resamp <- resample(rast_result, bnam_layers[[1]], method = "bilinear")
  return(raster_resamp)
  
}

vars <- colnames(select(ens_df_period_cell, contains("_")))

cmip_rast <- lapply(vars, transform_cmip_to_raster, period = "P1", ssp = "2-4.5")
names(cmip_rast) <- vars

names(bnam_layers)[1] <- case_when(str_detect(names(bnam_layers)[1]))

match_cmip_names <- function(bnam_name) {
  first_result <- case_when(
    str_detect(bnam_name, "b_cur") ~ str_replace(bnam_name, "NRA_BNAM_b_cur(_avg)*", "wobavg"),
    str_detect(bnam_name, "b_sal") ~ str_replace(bnam_name, "NRA_BNAM_b_sal(_avg)*", "sobavg"),
    str_detect(bnam_name, "b_stress") ~ str_replace(bnam_name, "NRA_BNAM_b_stress(_avg)*", "bstress"),
    str_detect(bnam_name, "b_tmp") ~ str_replace(bnam_name, "NRA_BNAM_b_tmp(_avg)*", "tobavg"),
    str_detect(bnam_name, "MLD_ann") ~ str_replace(bnam_name, "NRA_BNAM_MLD_ann", "mldavg"),
    str_detect(bnam_name, "MLD_fall") ~ str_replace(bnam_name, "NRA_BNAM_MLD_fall_10_12", "mldavg_F"),
    str_detect(bnam_name, "MLD_spr") ~ str_replace(bnam_name, "NRA_BNAM_MLD_spr_04_06", "mldavg_Sp"),
    str_detect(bnam_name, "MLD_sum") ~ str_replace(bnam_name, "NRA_BNAM_MLD_sum_07_09", "mldavg_Su"),
    str_detect(bnam_name, "MLD_win") ~ str_replace(bnam_name, "NRA_BNAM_MLD_win_01_03", "mldavg_W"),
    # str_detect(bnam_name, "s_cur") ~ str_replace(bnam_name, "NRA_BNAM_s_cur(_avg)*", "wobavg"),
    str_detect(bnam_name, "s_sal") ~ str_replace(bnam_name, "NRA_BNAM_s_sal(_avg)*", "sosavg"),
    str_detect(bnam_name, "s_tmp") ~ str_replace(bnam_name, "NRA_BNAM_s_tmp(_avg)*", "tosavg"),
    str_detect(bnam_name, "_ran") ~ NA_character_,
    TRUE ~ NA_character_
  )
  
  second_result <- ifelse(str_detect(first_result, "_ran"), NA_character_, first_result)
  return(second_result)
  
}

names(bnam_layers) <- sapply(names(bnam_layers), match_cmip_names)

# Remove any that didn't match
bnam_layers <- bnam_layers[!is.na(names(bnam_layers))]

common_names <- intersect(names(bnam_layers), names(cmip_rast))


# Calculate difference and correlation for each common variable
diff_list <- list()
cor_list <- list()
for (varstat in common_names) {
  bnam_raster <- bnam_layers[[varstat]]
  cmip_raster <- cmip_rast[[varstat]]
  
  # Ensure both rasters have the same extent and resolution
  if (!compareGeom(bnam_raster, cmip_raster, stopOnError = FALSE)) {
    cmip_raster <- resample(cmip_raster, bnam_raster, method = "bilinear")
  }
  
  # Calculate difference
  diff_raster <- bnam_raster - cmip_raster
  diff_list[[varstat]] <- diff_raster
  
  # # Calculate correlation (using cell values)
  # bnam_values <- values(bnam_raster, na.rm = TRUE)
  # cmip_values <- values(cmip_raster, na.rm = TRUE)
  # 
  # # Ensure both vectors are of the same length after removing NAs
  # valid_idx <- complete.cases(bnam_values, cmip_values)
  # correlation <- cor(bnam_values[valid_idx], cmip_values[valid_idx])
  # cor_list[[varstat]] <- correlation
}

calculate_rast_diff <- function(varstat) {
  bnam_raster <- bnam_layers[[varstat]]
  cmip_raster <- cmip_rast[[varstat]]
  
  # Ensure both rasters have the same extent and resolution
  if (!compareGeom(bnam_raster, cmip_raster, stopOnError = FALSE)) {
    cmip_raster <- resample(cmip_raster, bnam_raster, method = "bilinear")
  }
  
  # Calculate difference
  diff_raster <- bnam_raster - cmip_raster
  
  return(diff_raster)
}
diff_list <- lapply(common_names, calculate_rast_diff)
names(diff_list) <- common_names



# Plot difference rasters
plot_rast_diff <- function(rast, varname) {
  df <- as.data.frame(rast, xy = TRUE)
  colnames(df)[3] <- "difference"
  
  gradient_scale_limit <- max(abs(c(min(df$difference), max(df$difference))))
  
  p <- ggplot(df, aes(x = x, y = y)) +
    theme_bw() +
    geom_raster(aes(fill = difference)) +
    scale_fill_cmocean(name = "balance", na.value = "transparent",
                       limits = c(-gradient_scale_limit, gradient_scale_limit),  # Symmetric limits around zero
                       # Rescale the legend to be centred at zero
                       values = scales::rescale(c(-gradient_scale_limit,0,gradient_scale_limit))) +
                       # values = scales::rescale(c(min(df$difference),0,max(df$difference)))) +
    coord_fixed(xlim = c(sa_lims[1], sa_lims[2]), ylim = c(sa_lims[3], sa_lims[4]), expand = FALSE) +
    labs(title = paste("Difference (BNAM - CMIP) for", varname),
         x = "Longitude", y = "Latitude",
         fill = "Difference")
  
  return(p)
}

lapply(names(diff_list), function(var) {
  p <- plot_rast_diff(diff_list[[var]], var)
  ggsave(filename = paste0("output/07_BNAM_CMIP_Diff_Maps/BNAM_CMIP_Diff_", var, ".jpg"), plot = p, width = 8, height = 6)
})



# Correlation ----

# Convert raster to dataframe for correlation calculation
raster_to_df <- function(rast) {
  df <- as.data.frame(rast, xy = TRUE)
  colnames(df)[3] <- "value"
  return(df)
}

bnam_dfs <- lapply(bnam_layers[common_names], raster_to_df)
cmip_dfs <- lapply(cmip_rast[common_names], raster_to_df)

# Calculate correlation for each variable
cor_list <- list()
scatterplot_list <- list()
for (var in common_names) {
  bnam_df <- bnam_dfs[[var]]
  cmip_df <- cmip_dfs[[var]]
  
  # Merge dataframes on coordinates
  merged_df <- merge(bnam_df, cmip_df, by = c("x", "y"), suffixes = c("_bnam", "_cmip"))
  
  # Remove rows with NA values
  merged_df <- na.omit(merged_df)
  
  # Calculate correlation
  correlation <- cor(merged_df$value_bnam, merged_df$value_cmip)
  cor_list[[var]] <- correlation
  
  # Save scatterplot for each variable
  scatterplot_list[[var]] <- ggplot(merged_df, aes(x = value_bnam, y = value_cmip)) +
    geom_point(alpha = 0.3) +
    theme_bw() +
    labs(title = paste("Scatter plot of BNAM vs CMIP for", var),
         x = "BNAM Values",
         y = "CMIP Values") +
    # Add lines and add to legend
    geom_smooth(method = "lm", aes(col = "blue")) +
    # Add 1-to-1 line
    geom_abline(aes(color = "red", slope = 1, intercept = 0), linetype = "dashed") +  
    # Add correlation coefficient to the plot
    geom_text(x = Inf, y = Inf, 
              label = paste("Correlation:", round(correlation, 2)), 
              hjust = 1.1, vjust = 1.5, 
              size = 5, color = "blue") +
    scale_colour_manual(name = "Lines", values = c("blue" = "blue", "red" = "red"),
                        labels = c("Regression", "1:1 Line"))
  
  # Save plot
  ggsave(filename = paste0("output/08_BNAM_CMIP_Cor_Plots/BNAM_CMIP_Scatter_", var, ".jpg"), plot = scatterplot_list[[var]], width = 6, height = 6)
}

# Create summary plot of variables and correlations, in order of most to least correlated
cor_df <- enframe(unlist(cor_list), name = "Variable", value = "Correlation") %>%
  arrange(desc(Correlation))
cor_df$Variable <- factor(cor_df$Variable, levels = cor_df$Variable)

ggplot(cor_df, aes(x = Variable, y = Correlation)) +
  geom_col(fill = "black", alpha = 0.3) +
  theme_bw() +
  coord_flip() +
  labs(title = "Correlation between BNAM and CMIP variables",
       x = "Variable",
       y = "Correlation coefficient") +
  # Remove x buffer
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, by = 0.2)) +
  # Remove horizontal grid lines
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
ggsave(filename = "output/08_BNAM_CMIP_Cor_Plots/BNAM_CMIP_Cor_Summary.jpg", width = 8, height = 6)




