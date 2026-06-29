# Create table of mean variable values +/- SD for cells with overlapping PA and PresenceOnly data

# Presence and absence (original) ----
vme_pts_pa <- cmip_comb_df |>
  filter(VME_Group == vmeoi) |>
  select(all_of(selected_vme_vars)) |>
  as.data.frame() |>
  pivot_longer(cols = everything(), names_to = "var", values_to = "value") |>
  summarise(mean = mean(value), sd = sd(value), .by = var)

# Presence only (refugia) ----
vme_pts_pres <- cmip_comb_df |>
  filter(VME_Group == vmeoi, VME_P_A == "Presence") |>
  select(all_of(selected_vme_vars)) |>
  as.data.frame() |>
  pivot_longer(cols = everything(), names_to = "var", values_to = "value") |>
  summarise(mean = mean(value), sd = sd(value), .by = var)

# Create table ----
var_table <- bind_rows(vme_pts_pa, vme_pts_pres, .id = "Dataset") |>
  mutate(Dataset = replace_values(Dataset, "1" ~ "Presence + Absence" , "2" ~ "Presence Only"))

write_csv(var_table, file = paste0(output_folder,"/",vmeoi,"_SelectedVariables_PA_PresOnly_mean_sd.csv"))