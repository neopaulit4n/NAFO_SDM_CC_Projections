# Model building

cat("Building RF model and extracting performance metrics across folds...\n")
# Create the 10 folds ---*
set.seed(loop_seed)
folds <- caret::createFolds(vme_df$VME_P_A, k = 10, returnTrain = TRUE)

# Initialise storage for fold metrics ----
fold_metrics_list <- list()
fold_predictions_spatial_baseline <- list()
fold_predictions_spatial_baseline_reclass <- list()
fold_predictions_spatial_proj <- list()
fold_predictions_spatial_proj_reclass <- list()
fold_var_imp <- list()
fold_partialdep <- list()
fold_model <- list()

# Begin looping across folds ----
for (i in 1:10) {
  
  # Get fold indices ----
  train_idx <- folds[[i]]
  test_idx <- setdiff(1:nrow(vme_df), train_idx)
  
  # Train model on this fold ----
  cat("Training fold", i, "\n")  
  set.seed(loop_seed + i)
  fold_model[[i]] <- randomForest::randomForest(
    VME_P_A ~ .,
    data = vme_df[train_idx, ],
    ntree = 500,
    mtry = floor(sqrt(ncol(vme_df) - 1)),
    strata = vme_df$VME_P_A[train_idx],
    replace = FALSE,
    importance = TRUE
  )
  
  cat("  Retrieving fold predictions and metrics\n")
  # Get predictions on held-out test data ----
  rf_test_pred_prob <- predict(fold_model[[i]], 
    newdata = vme_df[test_idx, ], 
    type = "prob")
  
  # Extract probability of positive class (assuming second level) ----
  lev <- levels(vme_df$VME_P_A)
  pred_prob <- rf_test_pred_prob[, lev[2]]
  obs_numeric <- ifelse(vme_df$VME_P_A[test_idx] == lev[2], 1, 0)
  
  # Calculate optimal threshold using Sens=Spec method ----
  threshold_df <- data.frame(
    id = 1:length(test_idx),
    PA = ifelse(vme_df$VME_P_A[test_idx] == lev[2], 1, 0),
    predprob = rf_test_pred_prob[, lev[2]]
  )
  
  opttsh <- PresenceAbsence::optimal.thresholds(
    threshold_df,
    opt.methods = "Sens=Spec"  # if wanted to also test prevalence: c("Sens=Spec", "ObsPrev")
  ) %>%
    pull(predprob)
  
  # Apply optimal threshold ----
  optimal_pred <- ifelse(pred_prob >= opttsh, lev[2], lev[1])
  optimal_pred <- factor(optimal_pred, levels = lev)
  obs_factor <- vme_df$VME_P_A[test_idx]
  
  # Calculate confusion matrix ----
  cm <- caret::confusionMatrix(optimal_pred, obs_factor, positive = lev[2])
  
  # Extract and store metrics ----
  metrics <- unlist(c(cm$overall, cm$byClass))
  metrics["TSS"] <- cm$byClass["Sensitivity"] + cm$byClass["Specificity"] - 1
  metrics["OptThreshold"] <- opttsh
  metrics["Fold"] <- i
  
  fold_metrics_list[[i]] <- metrics

  # Fold variable importance ----
  cat("  Extracting variable importance\n")
  fold_var_imp[[i]] <- as.data.frame(randomForest::importance(fold_model[[i]])) %>%
    rownames_to_column(var = "Variable") %>%
    arrange(desc(MeanDecreaseGini))

  # Fold partial dependence data ----
  cat("  Extracting partial dependence data\n")
  fold_partialdep[[i]] <- lapply(selected_vme_vars, function(var) {
    pdp::partial(fold_model[[i]], 
      pred.var = var,
      plot = FALSE) %>%
      mutate(Variable = colnames(.)[1]) %>%
      rename(value = var)
  })
  
  cat("  Generating spatial predictions under baseline conditions\n")
  # Spatial predictions for this fold ----
  fold_predictions_spatial_baseline[[i]] <- terra::predict(
    vme_layers_baseline,
    fold_model[[i]],
    type = 'prob',
    na.rm = TRUE,
    index = 1:2
  )
  
  # Convert predictions to presence/absence using optimal threshold ----
  fold_predictions_spatial_baseline_reclass[[i]] <- terra::classify(
    fold_predictions_spatial_baseline[[i]][[2]], 
    rcl = matrix(c(-Inf, opttsh, 0,
                    opttsh, Inf, 1), 
                    ncol = 3, byrow = TRUE)
  )

  cat("  Generating spatial predictions under projected scenarios\n")
  # Spatial predictions for this fold ----
  fold_predictions_spatial_proj[[i]] <- lapply(period_all, function(poi) {
    lapply(ssp_all, function(sspoi) {
      terra::predict(
        vme_layers_proj[[poi]][[sspoi]],
        fold_model[[i]],
        type = 'prob',
        na.rm = TRUE,
        index = 1:2
      )      
    }) %>%
      set_names(ssp_all)
  }) %>%
    set_names(period_all)
  
  # Convert predictions to presence/absence using optimal threshold ----
  fold_predictions_spatial_proj_reclass[[i]] <- lapply(period_all, function(poi) {
    lapply(ssp_all, function(sspoi) {
      terra::classify(
        fold_predictions_spatial_proj[[i]][[poi]][[sspoi]][[2]], 
        rcl = matrix(c(-Inf, opttsh, 0,
                        opttsh, Inf, 1), 
                        ncol = 3, byrow = TRUE)
      )      
    }) %>% set_names(ssp_all)
  }) %>% set_names(period_all)

}

cat("Modelling complete. Processing results...\n")

# Combine all fold metrics into a dataframe ----
fold_metrics_df_i <- do.call(rbind, lapply(fold_metrics_list, function(x) {
  data.frame(t(x))
})) %>%
  pivot_longer(cols = -"Fold", names_to = "metric", values_to = "value") %>%
  mutate(VME_Group = vmeoi) %>%
  relocate(VME_Group, Fold, metric, value)
# fold_metrics_df <- bind_rows(fold_metrics_df, fold_metrics_df_i)

# Summarise metrics across folds ----
fold_metrics_summary_df_i <- fold_metrics_df_i %>%
  group_by(VME_Group, metric) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE),
            .groups = "drop")
fold_metrics_summary_df <- bind_rows(fold_metrics_summary_df, fold_metrics_summary_df_i)

# Create spatial predictions stack ----
rf_pred_foldstack_baseline <- terra::rast(fold_predictions_spatial_baseline_reclass)
# terra::writeRaster(rf_pred_foldstack_baseline, 
#   filename = paste0(main_output_folder, vmeoi, "/", vmeoi, "_rf_spatial_predictions_baseline.tif"), overwrite = TRUE)

rf_pred_foldstack_proj <- map(period_all, function(poi) {
  map(ssp_all, function(sspoi) {
    fold_layers <- map(fold_predictions_spatial_proj_reclass, ~ .x[[poi]][[sspoi]])
    terra::rast(fold_layers)
  }) %>% set_names(ssp_all)
}) %>% set_names(period_all) %>%
  unlist()


# Extract variable importance for final selected variables ----
cat("Extracting RF model variable importance metrics...\n")

fold_var_imp_df <- lapply(fold_var_imp, function(fold) {
  fold %>%
    filter(Variable %in% selected_vme_vars) %>%
    mutate(var_imp_rank = row_number())
}) %>%
  bind_rows(.id = "Fold") %>%
  mutate(Variable = fct_reorder(Variable, MeanDecreaseGini, .fun = mean)) %>%
  ungroup()
write.csv(fold_var_imp_df, 
  file = paste0(output_folder,"/",vmeoi,"_table_rf_VarImp.csv"), row.names = FALSE)

ggplot(fold_var_imp_df, aes(y = Variable, x = MeanDecreaseGini)) +
  geom_boxplot() +
  theme_bw() +
  labs(y = "Predictor Variable", x = "Mean Decrease in Gini Index")

ggsave(filename = paste0(output_folder,"/",vmeoi,"_plot_rf_VarImp.jpg"), 
  width = 6, height = 4)

# Extract partial dependence plots for each variable ----
fold_partial_df <- lapply(fold_partialdep, function(fold) {
  bind_rows(fold)
}) %>%
  bind_rows(.id = "Fold") %>%
  arrange(Fold, Variable, value) %>%
  mutate(Fold = as.factor(as.numeric(Fold)),
         Variable = factor(Variable, levels = rev(levels(fold_var_imp_df$Variable))))
write.csv(fold_partial_df, file = paste0(output_folder,"/",vmeoi,"_table_rf_PartialDep.csv"), row.names = FALSE)

ggplot(fold_partial_df, aes(x = value, y = yhat, colour = Fold)) +
  geom_line() +
  facet_wrap(~ Variable, scales = "free_x") +
  theme_bw() +
  labs(x = "Predictor Value", y = "Partial Dependence")

ggsave(filename = paste0(output_folder,"/",vmeoi,"_plot_rf_PartialDep.jpg"),
  width = 10, height = 8)

# Output non-reclassed/non-thresholded presence probability rasters ----

## Baseline ----
rf_res_presprob_baseline <- lapply(fold_predictions_spatial_baseline, `[[`, 2) %>%  # extract Presence layer only
  terra::rast(.) %>%
  terra::mean(.)

rf_res_absprob_baseline <- lapply(fold_predictions_spatial_baseline, `[[`, 1) %>%  # extract Absence layer only
  terra::rast(.) %>%
  terra::mean(.)

terra::writeRaster(rf_res_presprob_baseline, 
  paste0(output_folder,"/rasters/",vmeoi,"_rf_res_baseline_rawPresenceProb.tif"), overwrite = TRUE)
terra::writeRaster(rf_res_absprob_baseline, 
  paste0(output_folder,"/rasters/",vmeoi,"_rf_res_baseline_rawAbsenceProb.tif"), overwrite = TRUE)

## Projections ----
rf_res_presprob_proj <- map(period_all, function(poi) {
  map(ssp_all, function(sspoi) {
    
    fold_layers <- map(fold_predictions_spatial_proj, ~ .x[[poi]][[sspoi]])
    fold_layers <- lapply(fold_layers, `[[`, 2) %>%
      terra::rast(.) %>%
      terra::mean(.)

    terra::writeRaster(fold_layers, 
      paste0(output_folder,"/rasters/",vmeoi,"_rf_res_proj_rawPresenceProb_",poi,"_",sspoi,".tif"), overwrite = TRUE)
    
    return(fold_layers)
  }) %>% set_names(ssp_all)
}) %>% set_names(period_all) %>%
  unlist()

rf_res_absprob_proj <- map(period_all, function(poi) {
  map(ssp_all, function(sspoi) {
    
    fold_layers <- map(fold_predictions_spatial_proj, ~ .x[[poi]][[sspoi]])
    fold_layers <- lapply(fold_layers, `[[`, 1) %>%
      terra::rast(.) %>%
      terra::mean(.)

    terra::writeRaster(fold_layers, 
      paste0(output_folder,"/rasters/",vmeoi,"_rf_res_proj_rawAbsenceProb_",poi,"_",sspoi,".tif"), overwrite = TRUE)
    
    return(fold_layers)
  }) %>% set_names(ssp_all)
}) %>% set_names(period_all) %>%
  unlist()


# Calculate spatial metrics across folds ----

## Baseline ----
### Most frequent class (0/1)
rf_res_baseline_MaxClass <- terra::modal(rf_pred_foldstack_baseline, freq = FALSE)

### Frequency of most frequent class (fraction of runs)
# rf_res_MaxClassF <- terra::modal(rf_pred_foldstack, freq = TRUE) / 10  # old method - doesn't work
rf_res_baseline_freq_count <- sum(rf_pred_foldstack_baseline == rf_res_baseline_MaxClass, na.rm = TRUE)
rf_res_baseline_MaxClassF <- rf_res_baseline_freq_count / 10

### Average probability of classes
rf_res_baseline_AvgProb <- Reduce("+", fold_predictions_spatial_baseline) / 10

### Average probability of maximum frequency class
rf_res_baseline_MaxClassAvgProb <- terra::selectRange(rf_res_baseline_AvgProb, rf_res_baseline_MaxClass + 1)

### Combined confidence metric
rf_res_baseline_CombConf <- rf_res_baseline_MaxClassF * rf_res_baseline_MaxClassAvgProb

### Number of models predicting presence
rf_res_baseline_CVSum <- terra::app(rf_pred_foldstack_baseline, sum, na.rm = TRUE)

### Save baseline rasters for each metric ----
rm(rf_res_baseline_freq_count, rf_res_baseline_AvgProb)
lapply(ls(pattern = "rf_res_baseline"), function(res_name) {
  res_raster <- get(res_name)
  terra::writeRaster(res_raster, 
    filename = paste0(output_folder,"/rasters/",vmeoi,"_",res_name,".tif"), overwrite = TRUE)
})

## Projections ----
for (i in 1:length(rf_pred_foldstack_proj)) {
  comb_name <- names(rf_pred_foldstack_proj)[[i]]

  ### Most frequent class (0/1)
  rf_res_proj_MaxClass <- terra::modal(rf_pred_foldstack_proj[[i]], freq = FALSE)
  # rf_res_proj_MaxClass_reclass <- terra::as.bool(rf_res_proj_MaxClass)

  ### Frequency of most frequent class (fraction of runs)
  # rf_res_MaxClassF <- terra::modal(rf_pred_foldstack, freq = TRUE) / 10  # old method - doesn't work
  rf_res_proj_freq_count <- sum(rf_pred_foldstack_proj[[i]] == rf_res_proj_MaxClass, na.rm = TRUE)
  rf_res_proj_MaxClassF <- rf_res_proj_freq_count / 10

  ### Average probability of classes
  # rf_res_proj_AvgProb <- map(period_all, function(poi) {
  #   map(ssp_all, function(sspoi) {
  #     fold_layers <- map(fold_predictions_spatial_proj, ~ .x[[poi]][[sspoi]])
  #     Reduce("+", fold_layers) / 10
  #   }) %>% set_names(ssp_all)
  # }) %>% set_names(period_all) %>%
  #   unlist()

  rf_res_proj_AvgProb <- map(fold_predictions_spatial_proj, 
    ~ .x[[str_extract(comb_name, "^P\\d")]][[str_extract(comb_name,"\\d-\\d\\.\\d$")]])
  rf_res_proj_AvgProb <- Reduce("+", rf_res_proj_AvgProb) / 10
  
  ### Average probability of maximum frequency class
  rf_res_proj_MaxClassAvgProb <- terra::selectRange(rf_res_proj_AvgProb, rf_res_proj_MaxClass + 1)

  ### Combined confidence metric
  rf_res_proj_CombConf <- rf_res_proj_MaxClassF * rf_res_proj_MaxClassAvgProb

  ### Number of models predicting presence
  rf_res_proj_CVSum <- terra::app(rf_pred_foldstack_proj[[i]], sum, na.rm = TRUE)
  
  ### Save projection rasters for each metric ----
  rm(rf_res_proj_freq_count, rf_res_proj_AvgProb)
  lapply(ls(pattern = "rf_res_proj"), function(res_name) {
    res_raster <- get(res_name)
    terra::writeRaster(res_raster, 
      filename = paste0(output_folder,"/rasters/",vmeoi,"_",res_name,"_",comb_name,".tif"), overwrite = TRUE)
  })

}
