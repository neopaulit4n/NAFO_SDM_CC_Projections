
# Looping through periods/SSPs/VME groups

# Load data ----
source("code/01_LoadData.R")

# Create table of each combination of iterations to save variable selection results ----
vme_all <- unique(resp_df$VME_Group)
period_all <- c("P1", "P2", "P3", "P4")
ssp_all <- unique(cmip_df_period_ssp$ssp)
# period_ssp_all <- apply(expand.grid(period_all, ssp_all), 1, paste, sep = "_", collapse = "_")

vmeoi <- "black_corals"

# Initialise fold metrics dataframe to save results from each fold of each iteration ----
fold_metrics_summary_df <- data.frame(
  VME_Group = character(),
  metric = character(),
  mean_value = numeric(),
  sd_value = numeric()
)

# Loop through each combination of VME group, period, SSP, and subsampling option ----
loop_seed <- 412
for (vmeoi in "black_corals") {  # vme_all

  # Create VMEOI directory if it doesn't exist already
  output_folder <- paste0("output/",vmeoi)  
  if (!dir.exists(output_folder)) dir.create(output_folder)
  if (!dir.exists(paste0(output_folder,"/rasters"))) dir.create(paste0(output_folder,"/rasters"))
  
  # Table to save variable selection results
  var_select_df <- expand_grid(vmeoi = vmeoi, poi = period_all, sspoi = ssp_all) %>%
    arrange(vmeoi, poi, sspoi) %>%
    # Add columns for CMIP variables
    cbind(matrix(NA, nrow = nrow(.), ncol = length(cmip_vars))) %>%
    set_names(c("vmeoi", "poi", "sspoi", cmip_vars))

  # Table to save VIF values
  vif_df <- data.frame(vmeoi = character(),
    poi = character(),
    sspoi = character(),
    variable = character(),
    vif = numeric())
  
  # Variable selection for this VME group
  cat("Running variable selection for VME group:", vmeoi, "\n")
  source("code/02_VariableSelection.R")
  
  # Modelling for this VME group
  cat("Running modelling for VME group:", vmeoi, "\n")
  source("code/03_Modelling.R")
  
  # Summary outputs by VME group
  var_select_df <- filter(var_select_df, vmeoi == vmeoi) %>%
    janitor::adorn_totals("row")  
  write_csv(var_select_df, paste0(output_folder,"/",vmeoi,"_summary_variable_selection.csv"))
  write_csv(vif_df, paste0(output_folder,"/",vmeoi,"_summary_vif_values.csv"))
  write_csv(fold_metrics_summary_df, paste0(output_folder,"/",vmeoi,"_summary_fold_metrics.csv"))

  # Extrapolation outputs loop
  for (i in 1:length(unlist(list(baseline = vme_layers_baseline, unlist(vme_layers_proj))))) {
    prediction_grid <- unlist(list(baseline = vme_layers_baseline, unlist(vme_layers_proj)))[[i]]
    output_name <- names(unlist(list(baseline = vme_layers_baseline, unlist(vme_layers_proj))))[i]
    cat(paste("Computing extrapolations for",names(unlist(list(baseline = vme_layers_baseline, unlist(vme_layers_proj))))[i]))
    source("code/05_Extrapolation.R")
  }  
  
  # Maps per VME group
  source("code/04_Mapping.R")

  # loop_seed <- loop_seed + 1     
}


