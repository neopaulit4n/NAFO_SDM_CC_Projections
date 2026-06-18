
# Looping through periods/SSPs/VME groups

# Load data ----
source("code/01_LoadData.R")

vmeoi <- "sea_pens"

# Initialise fold metrics dataframe to save results from each fold of each iteration ----
fold_metrics_summary_df <- data.frame(
  VME_Group = character(),
  metric = character(),
  mean_value = numeric(),
  sd_value = numeric()
)

# Loop through each combination of VME group, period, SSP, and subsampling option ----
loop_seed <- switch(vmeoi,
  "black_corals" = 412,
  "small_gorgonians" = 412,
  "large_sponges" = 413,
  "sea_pens" = 414
)

# for (vmeoi in vme_all) {

  ## Create VMEOI directory if it doesn't exist already ----
  output_folder <- paste0("output/",vmeoi)
  if (!dir.exists(output_folder)) dir.create(output_folder)
  if (!dir.exists(paste0(output_folder,"/RFModelRasters"))) dir.create(paste0(output_folder,"/RFModelRasters"))

  ## Table to save VIF values ----
  vif_df <- data.frame(vmeoi = character(),
    poi = character(),
    sspoi = character(),
    variable = character(),
    vif = numeric())
  
  ## Variable selection for this VME group ----
  cat("Running variable selection for VME group:", vmeoi, "\n")
  source("code/02_VariableSelection.R")
  
  ## Modelling for this VME group ----
  cat("Running modelling for VME group:", vmeoi, "\n")
  source("code/03_Modelling.R")
  
  ## Summary outputs by VME group ----
  # var_select_df <- filter(var_select_df, vmeoi == vmeoi) %>%
  #   janitor::adorn_totals("row")  
  # write_csv(var_select_df, paste0(output_folder,"/",vmeoi,"_summary_variable_selection.csv"))
  write_csv(vif_df, paste0(output_folder,"/",vmeoi,"_summary_vif_values.csv"))
  write_csv(fold_metrics_summary_df, paste0(output_folder,"/",vmeoi,"_summary_fold_metrics.csv"))

  ## Extrapolation outputs loop ----
  for (i in 1:length(unlist(list(baseline = vme_layers_baseline, unlist(vme_layers_proj))))) {
    prediction_grid <- unlist(list(baseline = vme_layers_baseline, unlist(vme_layers_proj)))[[i]]
    output_name <- names(unlist(list(baseline = vme_layers_baseline, unlist(vme_layers_proj))))[i]
    cat(paste("Computing extrapolations for",names(unlist(list(baseline = vme_layers_baseline, unlist(vme_layers_proj))))[i],"\n"))
    source("code/05_Extrapolation.R")
  }  

  ## Extrapolations for P1-4 for SSP 2-4.5 (specific figures within body of text) ----
  cat("Computing specific extrapolations for SSP 2-4.5\n")
  source("code/06_ExtrapolationsSSP245.R")

  ## Extrapolations percentages tables overall and by variable ----
  cat("Creating extrapolation percentages tables\n")
  source("code/12_ExtrapolationPercTables.R")
  
  ## Maps per VME group ----
  cat("Creating final modelling maps\n")
  source("code/04_Mapping.R")

  ## Output functional response curves ----
  cat("Creating functional response curves\n")
  source("code/09_FunctionalResponseCurves.R")

  ## Environmental variable layer maps ----
  cat("Outputting CMIP selected variables three-plot comparisons\n")
  source("code/08_SelVarsThreePlotComparison.R")

  ## Correlation matrix plots + correlation differences ----
  cat("Extracting correlation matrix tables, plots, and differences\n")
  source("code/11_VarCorrTimeSeries.R")

  ## Selected CMIP vars ADF + Tukey table ----
  cat("Outputting CMIP variable ADF & Tukey results table\n")
  source("code/10_CMIPADFTable.R")

  ## Re-run modelling without variable selection (all CMIP, still only top selected terrain variable) ----
  ### Prepare VME group dataframe ----
  vme_terrain_vars <- filter(terrain_topvars, VME_Group == vmeoi) %>%
    pull(variable)
  vme_vars <- c(vme_terrain_vars, names(cmip_layers))
  vme_df <- filter(cmip_comb_df, VME_Group == vmeoi) %>%
    select(all_of(c("VME_P_A", vme_vars)))

  selected_vme_vars <- vme_vars

  ### Run modelling ----
  source("code/07_ModellingNoVarSel.R")
  write_csv(fold_metrics_summary_df_i, paste0(output_folder,"/",vmeoi,"_summary_fold_metrics_novarsel.csv"))
    
# }


