
# Editing and transferring SDM2024 code for use with CC projection model data

library(tidyverse)

# Load data ----

# Load response
resp_df <- read_csv("data/cleaned/VME_group_PA_df.csv", show_col_types = FALSE)

# Load predictors

# Load terrain variables (unchanging)
bathy_layers <- list.files(path = "data/raw/BNAM_Data_From_Cam/Bathymetry_Terrain_From_NAFO_SharePoint", 
                           pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(., nm = basename(.) %>% tools::file_path_sans_ext()) %>%
  lapply(terra::rast) %>%
  lapply(terra::project, "EPSG:4326")

# Load BNAM layers (will use these to form predictions, decide which variables to select)
bnam_layers <- list.files("data/raw/BNAM_Data_From_Cam/BNAM_From_NAFO_SharePoint", 
                          pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(., nm = basename(.) %>% tools::file_path_sans_ext()) %>%
  lapply(terra::rast) %>%
  lapply(terra::project, "EPSG:4326")

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
    str_detect(bnam_name, "s_sal") ~ str_replace(bnam_name, "NRA_BNAM_s_sal(_avg)*", "sosavg"),
    str_detect(bnam_name, "s_tmp") ~ str_replace(bnam_name, "NRA_BNAM_s_tmp(_avg)*", "tosavg"),
    TRUE ~ NA_character_
  )
  
  second_result <- ifelse(str_detect(first_result, "_ran"), NA_character_, first_result)
  return(second_result)
  
}

names(bnam_layers) <- sapply(names(bnam_layers), match_cmip_names)
bnam_layers <- bnam_layers[!is.na(names(bnam_layers))]


# Extract predictor values at response locations ----
pred_values <- lapply(bnam_layers, function(layer) {
  terra::extract(layer, 
                 select(resp_df, Start_Long_DD, Start_Lat_DD)) %>%
    select(-ID)
}) %>%
  bind_cols()
colnames(pred_values) <- names(bnam_layers)

# Combine predictor and response dataframes ----
comb_df <- bind_cols(resp_df, pred_values) %>%
  mutate(VME_P_A = as.factor(VME_P_A))

# Remove NA
comb_df_compl <- comb_df %>%
  drop_na()

comb_df_miss <- comb_df[which(!complete.cases(comb_df)),]

# Preliminary RF model test ----
rf_test_form <- as.formula(paste("VME_P_A ~", 
                                 paste(colnames(pred_values), collapse = " + ")))
rf_test <- randomForest::randomForest(formula = rf_test_form,
                        data = filter(comb_df_compl, VME_Group == "boltenia"),
                        importance=TRUE)
rf_test <- lapply(unique(comb_df_compl$VME_Group), function(vme_group) {
  rf_model <- randomForest::randomForest(formula = rf_test_form,
                                        data = filter(comb_df_compl, VME_Group == vme_group),
                                        importance=TRUE)
  return(rf_model)
}) %>%
  set_names(unique(comb_df_compl$VME_Group))

rf_test
