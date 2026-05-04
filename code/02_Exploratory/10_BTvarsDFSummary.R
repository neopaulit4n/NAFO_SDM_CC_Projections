
df <- filter(cmip_comb_df, VME_Group == vmeoi) |>
  select(VME_P_A, Start_Long_DD, Start_Lat_DD, all_of(vme_vars)) |>
  select(VME_P_A, Start_Long_DD, Start_Lat_DD, starts_with("BT_")) |>
  filter(VME_P_A == "Presence") |>
  select(-c(VME_P_A))#, Start_Long_DD, Start_Lat_DD))
  # sf::st_as_sf(coords = c("Start_Long_DD","Start_Lat_DD"))

write_csv(df, paste0(output_folder,"/BTstats_VMEPresenceOnly_ReferenceOnly.csv"))

