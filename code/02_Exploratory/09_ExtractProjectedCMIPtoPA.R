# Extract raster values for PA for projected CMIP layers
suppressMessages(cmip_pred_proj_df <- lapply(unlist(cmip_layers_proj), function(layer) {
  terra::extract(layer, select(resp_df, Start_Long_DD, Start_Lat_DD)) %>%
  select(-ID)
}) %>%
  bind_cols() %>%
  set_names(names(unlist(cmip_layers_proj)))
)

# cmip_pred_proj_df_long <- bind_cols(resp_df, cmip_pred_proj_df) %>%
#   filter(VME_Group == vmeoi) %>%
#   pivot_longer(
#     cols = matches("^P\\d"),
#     names_to = c("period", "ssp", ".value"),
#     names_pattern = "^(P\\d)\\.(\\d-\\d+\\.\\d)\\.(.+)$"
#   )

z <- cmip_pred_df |>
  select(all_of(selected_cmip_vars)) |>
  # mutate(period = "baseline")
  rename_with(~ paste0("P0.0-0.0.", .x))
zz <- cmip_pred_proj_df |>
  # select(VME_Group:ssp, all_of(selected_cmip_vars))
  select(all_of(colnames(cmip_pred_proj_df)[grepl(paste(selected_cmip_vars, collapse = '|'), colnames(cmip_pred_proj_df))]))

cmip_pred_df_selvars <- bind_cols(z,zz) |>
  pivot_longer(
    cols = matches("^P\\d"),
    names_to = c("period", "ssp", ".value"),
    names_pattern = "^(P\\d)\\.(\\d-\\d+\\.\\d)\\.(.+)$"
  ) |>
  mutate(
    period = ifelse(period == "P0", "Reference", period),
    ssp = ifelse(ssp == "0-0.0", "Reference", ssp)
  ) |>
  pivot_longer(cols = -c("period","ssp"), names_to = "variable", values_to = "value") |>
  summarise(
    mean = mean(value, na.rm = TRUE), 
    sd = sd(value, na.rm = TRUE), 
    .by = c(period, ssp, variable)
  )


# write_csv(cmip_pred_proj_df, "output/01_Exploratory/ProjectedCMIPtoBlackCoralsPA_df.csv")
