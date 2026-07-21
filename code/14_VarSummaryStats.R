# Create summary stats table of specific variable for each period/SSP

# Presence only points dataframe
vme_pts_pres <- cmip_comb_df %>%
  filter(VME_Group == vmeoi,
         VME_P_A == "Presence") %>%
  select(x = Start_Long_DD, y = Start_Lat_DD) %>%
  as.data.frame() %>%
  terra::vect()

baseline_layers <- c(cmip_layers, bathy_layers[vme_terrain_vars_start]) |>
  terra::rast()

projection_layers <- unlist(cmip_layers_proj, recursive = FALSE) |>
  lapply(X = _, 
    function(layer) {
      rast_result <- c(layer, bathy_layers[vme_terrain_vars_start]) |>
        terra::rast()
      names(rast_result) <- gsub("GEBCO2024_FS005_StudyArea_","", names(rast_result))
      return(rast_result)
    })

all_layers <- list()
all_layers[[1]] <- baseline_layers
names(all_layers) <- "baseline"
all_layers <- c(all_layers, projection_layers)

var_summ <- lapply(1:length(all_layers), \(x) {
  terra::extract(all_layers[[x]], vme_pts_pres) %>%
    select(BT_min) %>%
    summarise(mean = mean(BT_min, na.rm = TRUE), min = min(BT_min, na.rm = TRUE), max = max(BT_min, na.rm = TRUE), sd = sd(BT_min, na.rm = TRUE)) |>
    mutate(var = "BT_min", PerSSP = names(all_layers)[[x]]) |>
    select(var, PerSSP, everything())
}) |>
  bind_rows() |>
  mutate(across(where(is.numeric), round, 2), PerSSP = ifelse(PerSSP == "baseline","Reference",PerSSP))

write_csv(var_summ, file = paste0(output_folder,"/",vmeoi,"_BT_min_PresOnly_SummaryStatsPeriodSSP.csv"))