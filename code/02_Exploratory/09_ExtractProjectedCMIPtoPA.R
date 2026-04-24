
# Extract raster values for PA for projected CMIP layers
suppressMessages(cmip_pred_proj_df <- lapply(unlist(cmip_layers_proj), function(layer) {
  terra::extract(layer, select(resp_df, Start_Long_DD, Start_Lat_DD)) %>%
  select(-ID)
}) %>%
  bind_cols() %>%
  set_names(names(unlist(cmip_layers_proj)))
)

cmip_pred_proj_df1 <- bind_cols(resp_df, cmip_pred_proj_df) %>%
  filter(VME_Group == "black_corals") %>%
  pivot_longer(
    cols = matches("^P\\d"),
    names_to = c("period", "ssp", ".value"),
    names_pattern = "^(P\\d)\\.(\\d-\\d+\\.\\d)\\.(.+)$"
  )

write_csv(cmip_pred_proj_df1, "output/01_Exploratory/ProjectedCMIPtoBlackCoralsPA_df.csv")