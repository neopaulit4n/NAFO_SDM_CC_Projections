
# Creating functional response curves similar to Javier's Paragorgia and Vazella papers' functional response curves

# Extract predicted presence probability values for each fold
# Extract environmental values
# Plot relationship using loess smoother

frc_baseline <- bind_cols(
  terra::as.data.frame(vme_layers_baseline),
  terra::as.data.frame(rf_res_presprob_baseline)
) %>%
  pivot_longer(-mean, names_to = "Variable", values_to = "Value")

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
          paste0(period_all[poi],"\\.",ssp_all[sspoi]), 
          names(rf_res_presprob_proj))
        ]] |>
          terra::as.data.frame()
      )
  })
  pred_df[[5]] <- frc_baseline
  pred_df <- pred_df |>
    set_names(c("Reference", period_all)) |>
    bind_rows(.id = "Period")
}) |>
  set_names(ssp_all) |>
  bind_rows(.id = "SSP") |>
  pivot_longer(-c(mean, SSP, Period), names_to = "Variable", values_to = "Value")

ggplot(z, aes(x = Value, y = mean, colour = Period)) +
  # facet_wrap(~Variable, scales = "free") +
  facet_grid(vars(SSP), vars(Variable), scales = "free") +
  geom_smooth(alpha = 0.7) +  # uses GAM by default
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
