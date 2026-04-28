
# Creating functional response curves similar to Javier's Paragorgia and Vazella papers' functional response curves

# Extract predicted presence probability values for each fold
# Extract environmental values
# Plot relationship using loess smoother

# Read in raster layers ----
rf_res_presprob_baseline <- terra::rast(paste0(output_folder,"/rasters/",vmeoi,"_rf_res_baseline_rawPresenceProb.tif"))

rf_res_presprob_proj <- lapply(
  list.files(
    path = paste0(output_folder,"/rasters"),
    pattern = "rf_res_proj_rawPresenceProb",
    full.names = TRUE),
  function(x) {
    terra::rast(x)
}) |>
  set_names(str_extract(
    list.files(path = paste0(output_folder,"/rasters"), pattern = "rf_res_proj_rawPresenceProb"),
    pattern = "P\\d_\\d-\\d\\.\\d"
  ))

# Load NAFO division areas ----
nafo_div <- sf::st_read("data/raw/Mapping_Layers/NAFO_Divisions/NAFO_Divisions_SHP/NAFO_Divisions_2021_poly_not_clipped.shp") |>
  sf::st_transform(4326) |>
  select(Division = Label, geometry)

frc_baseline <- bind_cols(
  terra::as.data.frame(vme_layers_baseline),
  terra::as.data.frame(rf_res_presprob_baseline, xy = TRUE)
) |>
  # pivot_longer(-c(mean, x, y), names_to = "Variable", values_to = "Value") #|>
  # Join NAFO divisions labels to points
  sf::st_as_sf(coords = c("x","y"), crs = 4326) |>
  sf::st_join(nafo_div)

ggplot(data = frc_baseline, aes(x = Value, y = mean_pres_prob)) +
  facet_wrap(~ Variable, scales = "free") +
  # geom_point() +
  # geom_smooth(method = "loess", se = FALSE) +  # takes a long time
  geom_smooth() +  # uses GAM by default
  theme_classic() +
  title(main = "Reference")

frc_proj <- lapply(1:length(rf_res_presprob_proj), function(x) {
  pred_df <- unlist(vme_layers_proj, recursive = FALSE)[[1]] |>
    terra::as.data.frame() |>
    drop_na()
  presprob_df <- terra::as.data.frame(rf_res_presprob_proj[[x]])

  comb_df <- bind_cols(pred_df, presprob_df) |>
    pivot_longer(-mean, names_to = "Variable", values_to = "Value")
  return(comb_df)
}) |>
  set_names(names(rf_res_presprob_proj))

frc_list <- c(baseline = list(frc_baseline), frc_proj)

lapply(1:length(frc_list), function(x) {
  ggplot(data = frc_list[[x]], aes(x = Value, y = mean)) +
    facet_wrap(~ Variable, scales = "free") +
    geom_smooth() +  # uses GAM by default
    theme_classic() +
    labs(
      title = names(frc_list)[x],
      y = "Mean predicted presence probability"
    )
  ggsave(paste0(output_folder,"/FunctionalResponseCurves/",vmeoi,"_FunctionalResponseCurveGAM_",names(frc_list)[x],".jpg"))
})

frc_df <- lapply(1:4, function(sspoi) {
  pred_df <- vme_layers_proj[[sspoi]]
  pred_df <- lapply(1:4, function(poi) {
    terra::as.data.frame(pred_df[[poi]]) |>
      drop_na() |>
      bind_cols(
        rf_res_presprob_proj[[grep(
          paste0(period_all[poi],"_",ssp_all[sspoi]), 
          names(rf_res_presprob_proj))
        ]] |>
          terra::as.data.frame(xy = TRUE) |>
          # Join NAFO division labels to points
          sf::st_as_sf(coords = c("x","y"), crs = 4326) |>
          sf::st_join(nafo_div)
      )
  })
  pred_df[[5]] <- frc_baseline
  pred_df <- pred_df |>
    set_names(c("Reference", period_all)) |>
    bind_rows(.id = "Period")
}) |>
  set_names(ssp_all) |>
  bind_rows(.id = "SSP") |>
  select(-geometry) |>
  pivot_longer(-c(mean, SSP, Period, Division), names_to = "Variable", values_to = "Value")

ggplot(frc_df, aes(x = Value, y = mean)) +
  facet_grid(vars(SSP), vars(Variable), scales = "free") +
  geom_point(aes(colour = Division), alpha = 0.2) +
  geom_smooth(alpha = 0.7, colour = "black") +  # uses GAM by default
  # stat_summary(
  #   fun.data = function(x) {
  #     data.frame(
  #       y = mean(x, na.rm = TRUE),
  #       ymin = mean(x, na.rm = TRUE) - 1.96 * sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)])),
  #       ymax = mean(x, na.rm = TRUE) + 1.96 * sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))
  #     )},
  #     geom = "ribbon", alpha = 0.2, colour = NA) +
  # Apply manual colour scheme for periods
  scale_colour_manual(
    "Period",
    values = c(
      "Reference" = "black",
      "P1" = "#05B", 
      "P2" = "darkgreen", 
      "P3" = "darkorange", 
      "P4" = "red")) +
  # scale_fill_manual(
  #   "Period",
  #   values = c(
  #     "Reference" = "black",
  #     "P1" = "#05B", 
  #     "P2" = "darkgreen", 
  #     "P3" = "darkorange", 
  #     "P4" = "red")) +
  theme_bw() +
  labs(y = "Mean predicted presence probability")

ggsave(paste0(output_folder,"/",vmeoi,"_FunctionalResponseCurveGAM.jpg"),
  width = 12, height = 8, dpi = 500)

# Plots with points coloured by NAFO division ----
z <- filter(frc_df, Period == "Reference", SSP == "1-2.6")
ggplot(z, aes(x = Value, y = mean)) +
  facet_wrap(~ Variable, scales = "free") +
  geom_point(aes(colour = Division), alpha = 0.2) +
  geom_smooth(alpha = 0.7, colour = "black") +  # uses GAM by default
  # stat_summary(
  #   fun.data = function(x) {
  #     data.frame(
  #       y = mean(x, na.rm = TRUE),
  #       ymin = mean(x, na.rm = TRUE) - 1.96 * sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)])),
  #       ymax = mean(x, na.rm = TRUE) + 1.96 * sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))
  #     )},
  #     geom = "ribbon", alpha = 0.2, colour = NA) +
  theme_bw() +
  labs(y = "Mean predicted presence probability", title = "Reference")
ggsave(paste0(output_folder,"/",vmeoi,"_FunctionalResponseCurveGAM_ReferencePts.jpg"),
  width = 12, height = 8, dpi = 500)
