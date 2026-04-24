
# Creating functional response curves similar to Javier's Paragorgia and Vazella papers' functional response curves

# Extract predicted presence probability values for each fold
# Extract environmental values
# Plot relationship using loess smoother

frc_df <- bind_cols(
  terra::as.data.frame(vme_layers_baseline),
  terra::as.data.frame(rf_res_presprob_baseline)
) %>%
  rename(mean_pres_prob = mean) %>%
  pivot_longer(-mean_pres_prob, names_to = "Variable", values_to = "Value")

ggplot(data = frc_df, aes(x = Value, y = mean_pres_prob)) +
  facet_wrap(~ Variable, scales = "free") +
  # geom_point() +
  # geom_smooth(method = "loess", se = FALSE) +  # takes a very long time
  geom_smooth() +  # uses GAM by default
  theme_classic()

