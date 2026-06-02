
# Creating functional response curves similar to Javier's Paragorgia and Vazella papers' functional response curves

# Extract predicted presence probability values for each fold
# Extract environmental values
# Plot relationship using loess smoother

# Read in raster layers ----
rf_res_presprob_baseline <- terra::rast(paste0(output_folder,"/RFModelRasters/",vmeoi,"_rf_res_baseline_rawPresenceProb.tif"))

rf_res_presprob_proj <- lapply(
  list.files(
    path = paste0(output_folder,"/RFModelRasters"),
    pattern = "rf_res_proj_rawPresenceProb",
    full.names = TRUE),
  function(x) {
    terra::rast(x)
}) |>
  set_names(str_extract(
    list.files(path = paste0(output_folder,"/RFModelRasters"), pattern = "rf_res_proj_rawPresenceProb"),
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
  # pivot_longer(-c(mean, x, y), names_to = "Variable", values_to = "Value") |>  # sometimes pivot, sometimes not
  # Join NAFO divisions labels to points
  sf::st_as_sf(coords = c("x","y"), crs = 4326) |>
  sf::st_join(nafo_div)

# ggplot(data = frc_baseline, aes(x = Value, y = mean)) +
#   facet_wrap(~ Variable, scales = "free") +
#   geom_point() +
#   # geom_smooth(method = "loess", se = FALSE) +  # takes a long time
#   geom_smooth() +  # uses GAM by default
#   theme_classic()

# frc_proj <- lapply(1:length(rf_res_presprob_proj), function(x) {
#   pred_df <- unlist(vme_layers_proj, recursive = FALSE)[[1]] |>
#     terra::as.data.frame() |>
#     drop_na()
#   presprob_df <- terra::as.data.frame(rf_res_presprob_proj[[x]])

#   comb_df <- bind_cols(pred_df, presprob_df) |>
#     pivot_longer(-mean, names_to = "Variable", values_to = "Value")
#   return(comb_df)
# }) |>
#   set_names(names(rf_res_presprob_proj))

# frc_list <- c(baseline = list(frc_baseline), frc_proj)

# lapply(1:length(frc_list), function(x) {
#   ggplot(data = frc_list[[x]], aes(x = Value, y = mean)) +
#     facet_wrap(~ Variable, scales = "free") +
#     geom_smooth() +  # uses GAM by default
#     theme_classic() +
#     labs(
#       title = names(frc_list)[x],
#       y = "Mean predicted presence probability"
#     )
#   ggsave(paste0(output_folder,"/FunctionalResponseCurves/",vmeoi,"_FunctionalResponseCurveGAM_",names(frc_list)[x],".jpg"))
# })

frc_list <- lapply(1:4, function(poi) {
  pred_df <- vme_layers_proj[[poi]]
  pred_df <- lapply(1:4, function(sspoi) {
    terra::as.data.frame(pred_df[[sspoi]]) |>
      drop_na() |>
      bind_cols(
        rf_res_presprob_proj[[grep(
          paste0(period_all[poi],"_",ssp_all[sspoi]), 
          names(rf_res_presprob_proj))
        ]] |>
          terra::as.data.frame(xy = TRUE) |>
          # Join NAFO division labels to points
          sf::st_as_sf(coords = c("x","y"), crs = 4326) |>
          sf::st_join(nafo_div) |>
          mutate(Period = poi)
      )
  })
  # pred_df[[5]] <- mutate(frc_baseline, Period = "Reference")
  pred_df <- pred_df |>
    set_names(ssp_all) #|>
    # bind_rows(.id = "SSP")
}) |>
  set_names(period_all) #|>
  # bind_rows(.id = "Period") |>
  # select(-geometry) |>
  # pivot_longer(-c(mean, SSP, Period, Division), names_to = "Variable", values_to = "Value")

frc_list[[5]] <- list(frc_baseline, frc_baseline, frc_baseline, frc_baseline) |>
  set_names(ssp_all)
names(frc_list) <- c(period_all, "Reference")

frc_df <- lapply(frc_list, bind_rows, .id = "SSP") |>
  bind_rows(.id = "Period") |>
  select(-geometry) |>
  pivot_longer(-c(mean, SSP, Period, Division), names_to = "Variable", values_to = "Value")

ggplot(frc_df, aes(x = Value, y = mean, colour = Period)) +
  facet_grid(vars(SSP), vars(Variable), scales = "free") +
  # geom_point(aes(colour = Division), alpha = 0.2) +
  geom_smooth(alpha = 0.7) +  # uses GAM by default
  # Apply manual colour scheme for periods
  scale_colour_manual(
    "Period",
    values = c(
      "Reference" = "black",
      "P1" = "#05B", 
      "P2" = "darkgreen", 
      "P3" = "darkorange", 
      "P4" = "red")) +
  theme_bw() +
  labs(y = "Mean predicted presence probability")
ggsave(paste0(output_folder,"/",vmeoi,"_FunctionalResponseCurveGAM.jpg"),
  width = 12, height = 8, dpi = 300)

# Plots with points coloured by NAFO division ----
z <- filter(frc_df, Period == "Reference", SSP == "1-2.6")
hull_df <- z |>
  group_by(SSP, Period, Division, Variable) |>
  slice(chull(Value, mean))
ggplot(z, aes(x = Value, y = mean)) +
  facet_wrap(~ Variable, scales = "free") +
  geom_point(aes(colour = Division), alpha = 0.05) +
  # geom_polygon(data = hull_df, aes(colour = Division), fill = "transparent", show.legend = FALSE) +  
  geom_smooth(alpha = 0.7, colour = "black") +  # uses GAM by default
  theme_bw() +
  labs(y = "Mean predicted presence probability", title = "Reference") +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
ggsave(paste0(output_folder,"/",vmeoi,"_FunctionalResponseCurveGAM_ReferencePts.jpg"),
  width = 12, height = 8, dpi = 300)
