
# Starting models

# Run one of models in full with the new variables for Period 1 and one of the SSPs 
# Then compare it to the BNAM results

# Load data ----
cat("Loading data...\n")
# vmeoi <- "black_corals"
# poi <- "P1"
# sspoi <- "1-2.6"

## Transform relevant CMIP variable data to raster layers
cmip_layers <- lapply(cmip_vars, function(var) {
  transform_cmip_to_raster(data = cmip_df_period_ssp, poi = poi, sspoi = sspoi, varstat = var)
})
names(cmip_layers) <- cmip_vars
# terra::plot(cmip_layers[[9]])

## Extract data from raster layers to the VME response points
suppressMessages(cmip_pred_df <- lapply(c(bathy_layers, cmip_layers), function(layer) {
  terra::extract(layer, 
                 select(resp_df, Start_Long_DD, Start_Lat_DD)) %>%
    select(-ID)
}) %>%
  bind_cols() %>%
  set_names(c(names(bathy_layers), names(cmip_layers)))
)

# Combine predictor and response dataframes ----
cmip_comb_df <- bind_cols(resp_df, cmip_pred_df) %>%
  mutate(VME_P_A = factor(VME_P_A, levels = c(0, 1), labels = c("Absence", "Presence"))) %>%
  drop_na()

# Replicating previous SDM variable elimination approach ----
# 1) Run preliminary RF models with all variables
# 2) Determine variable importance
# 3) Remove correlated variables based on importance ranking
# 4) VIF remaining variables


## Prepare VME group dataframe ----
vme_terrain_vars <- filter(terrain_topvars, VME_Group == vmeoi) %>%
  pull(variable)
vme_vars <- c(vme_terrain_vars, names(cmip_layers))
vme_df <- filter(cmip_comb_df, VME_Group == vmeoi) %>%
  select(all_of(c("VME_P_A", vme_vars)))


## Preliminary RF model formula ----  
rf_prelim_form <- as.formula(paste("VME_P_A ~", 
                                  paste(vme_vars, collapse = " + ")))

## Preliminary RF model run (save/load) ----
cat("Running preliminary RF model...\n")
set.seed(loop_seed)
rf_prelim <- randomForest::randomForest(formula = rf_prelim_form,
                                        data = vme_df,
                                        importance = TRUE)

# saveRDS(rf_model, file = paste0("data/processed/rf_prelim_",vmeoi,".rds"))
# rf_prelim <- readRDS(paste0("data/processed/rf_prelim_",vmeoi,".rds"))

## Extract variable importance metrics ----
cat("Extracting preliminary RF model variable importance metrics...\n")
rf_prelim_imp <- as.data.frame(randomForest::importance(rf_prelim)) %>%
  rownames_to_column(var = "Variable") %>%
  arrange(desc(MeanDecreaseGini))

# Show var imp plot for each VME group, order each by descending importance
# plot_rf_prelim_var_imp <- ggplot(rf_prelim_imp, aes(x = reorder(Variable, MeanDecreaseGini), 
#                           y = MeanDecreaseGini)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +
#   theme_bw() +
#   labs(x = "Predictor Variable", y = "Mean Decrease in Gini Index") +
#   theme(legend.position = "none")

# ggsave(plot_rf_prelim_var_imp, 
#   filename = paste0("output/02_Modelling_Outputs/",
#     vmeoi, "_", poi, "_", sspoi, "_plot_rf_prelim_var_imp.png"), 
#   width = 6, height = 4)

## Plot partial dependence plots ----


## Variable correlations ----

# Calculate VIF for selected variables for each VME group
# Check each VME group's VIF values and see if any are > 10; identify VME groups that need to be re-evaluated
# VIF values need to be < 10
# If values are > 10, need to recompute Spearman correlation at lower threshold 
#   (increments of 0.05) and re-run variable selection until all variables achieve VIF < 10.

### Plot correlation matrix prior to variables selection ----
cat("Calculating variable correlations and plotting correlation matrix...\n")
cor_df <- cor(vme_df[, vme_vars, drop = FALSE], 
                  method = "spearman",
                  use = "pairwise.complete.obs") %>%
  as.data.frame %>%
  pivot_longer(everything(), names_to = "var2", values_to = "cor") %>%
  mutate(var1 = rep(vme_vars, each = length(vme_vars))) %>%
  select(var1, var2, cor)
# write.csv(cor_df, file = paste0("output/02_Modelling_Outputs/",
#   paste(vmeoi, poi, sspoi, "table_cor_AllCMIPVars", sep = "_"), ".csv"), row.names = FALSE)

# plot_cor_allvars <- ggplot(data = cor_df, aes(x = var1, y = var2, fill = cor)) +
#   geom_tile(colour = "black") +
#   geom_text(aes(label = ifelse(cor > 0.7 | cor < -0.7, "*", "")), size = 3) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", limit = c(-1,1), midpoint = 0) +
#   scale_x_discrete(expand = c(0,0)) +
#   scale_y_discrete(expand = c(0,0)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.line = element_blank(),
#         axis.title = element_blank()) +
#   labs(title = paste("Correlation Matrix for VME Group:", vmeoi), fill = "Correlation")

# ggsave(plot_cor_allvars,
#   filename = paste0("output/02_Modelling_Outputs/",
#     vmeoi, "_", poi, "_", sspoi, "_plot_cor_AllCMIPVars.png"), 
#   width = 8, height = 6)


### Remove correlated variables based on importance ranking and VIF values ----
cat("Selecting variables...\n")
select_vme_vars <- function(pred_var_df, vmeoi, cor_threshold = 0.7, verbose = FALSE) {
  
  # cat("\n==============================\n")
  # cat("\nSelecting variables for VME group:", vmeoi, "\n")
  
  vif_threshold_met <- FALSE
  cor_threshold_new <- cor_threshold
  
  while (!vif_threshold_met) {
    
    # Calculate Spearman correlations
    remove_cor_variables <- function(data, priority_list, current_threshold = 0.7, verbose = FALSE) {
      
      vars <- colnames(data)
      
      # Calculate correlation
      cor_matrix <- cor(data, method = "spearman")
      
      # Identify highly correlated pairs and remove lower priority variables
      removed_vars <- c()
      
      for (i in 1:length(priority_list)) {
        var_i <- priority_list[i]
        if (verbose) cat("\nEvaluating variable:", var_i, "\n")
        
        # Skip if this variable was already removed
        if (var_i %in% removed_vars) {
          if (verbose) cat("  Skipped (already removed)\n")
          next
        }
        
        for (j in i:length(priority_list)) {
          var_j <- priority_list[j]
          if (verbose) cat("  Comparing with variable:", var_j, "\n")
          
          # Skip if same variable or already removed
          if (i == j || var_j %in% removed_vars) {
            if (verbose) cat("    Skipped (same variable or already removed)\n")
            next
          }
          
          # Check correlation/association
          if (abs(cor_matrix[var_i, var_j]) > current_threshold) {
            # Remove the lower priority variable
            removed_vars <- c(removed_vars, var_j)
            if (verbose) cat("    Removed", var_j, "due to high association with", var_i,
                            "(", round(cor_matrix[var_i, var_j], 2), ")\n")
          }
        }
      }
      
      # Return remaining variables
      return(setdiff(colnames(data), removed_vars))
    }
    
    # cat("\n   Removing correlated variables (threshold:", cor_threshold_new, ")\n")
    
    uncor_vars <- remove_cor_variables(pred_var_df, 
                                        priority_list = rf_prelim_imp %>%
                                        pull(Variable),
                                        current_threshold = cor_threshold_new,
                                        verbose = verbose)
    
    # Calculate VIF values
    # cat("\n   Calculating VIF values for selected variables\n")
    
    corvif <- function(data, uncorrelated = uncor_vars, verbose = FALSE) {
      data <- as.data.frame(data)
      
      form <- formula(paste("fooy ~ ", paste(strsplit(names(data), " "), collapse = " + ")))
      data <- data.frame(fooy = 1 + rnorm(nrow(data)), data)
      lm_mod <- lm(form, data)

      # Removed aliased variables prior to VIF
      aliased_vars <- alias(lm_mod)
      uncor_unaliased_vars <- uncorrelated[!(uncorrelated %in% rownames(aliased_vars$Complete))]

      # New model with uncorrelated variables
      form_new <- formula(paste("fooy ~ ", paste(uncor_unaliased_vars, collapse = " + ")))
      data_new <- data.frame(fooy = 1 + rnorm(nrow(data)), data[, uncor_unaliased_vars, drop = FALSE])
      lm_mod <- lm(form_new, data_new)
      
      # Calculate VIF
      VIF_result <- data.frame(vif = car::vif(lm_mod))
      if (verbose) print(VIF_result)
      return(VIF_result)
    }
    vif_values <- corvif(pred_var_df[, uncor_vars], verbose = verbose)
    
    # Check if all VIF values are < 10
    if (all(vif_values$vif < 10)) {
      # cat("\n     All VIF values are below threshold for VME group:", vmeoi, "\n")
      vif_threshold_met <- TRUE
      return(list(selected_vars = uncor_vars,
                  vif_values = vif_values,
                  final_cor_threshold = cor_threshold_new))
    }
    else {
      # cat("\n     VIF values exceed threshold (10); lowering correlation threshold and re-evaluating.\n")
      cor_threshold_new <- cor_threshold_new - 0.05
    }
  }
}

vme_var_selection <- select_vme_vars(vme_df[,-1], vmeoi = vmeoi, cor_threshold = 0.7, verbose = FALSE)
selected_vme_vars <- vme_var_selection$selected_vars

#### Update VIF table with this iteration's results
vif_df_i <- data.frame(vmeoi = vmeoi, poi = poi, sspoi = sspoi, 
  variable = rownames(vme_var_selection$vif_values), vif = vme_var_selection$vif_values$vif)
vif_df <- rbind(vif_df, vif_df_i)


# Create final dataframe with selected variables ----
vme_df <- vme_df %>%
  select(all_of(c("VME_P_A", selected_vme_vars)))

# Correlation matrix for selected variables
# cor_df <- cor(vme_df[, selected_vme_vars, drop = FALSE], 
#                   method = "spearman",
#                   use = "pairwise.complete.obs") %>%
#   as.data.frame %>%
#   pivot_longer(everything(), names_to = "var2", values_to = "cor") %>%
#   mutate(var1 = rep(selected_vme_vars, each = length(selected_vme_vars))) %>%
#   select(var1, var2, cor)
# write.csv(cor_df, file = paste0("output/02_StepByStepOutputs/",
#   paste(vmeoi, poi, sspoi, "table_cor_SelCMIPVars", sep = "_"), ".csv"), row.names = FALSE)

# plot_cor_selvars <- ggplot(data = cor_df, aes(x = var1, y = var2, fill = cor)) +
#   geom_tile(colour = "black") +
#   geom_text(aes(label = ifelse(cor > 0.7 | cor < -0.7, "*", "")), size = 3) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", limit = c(-1,1), midpoint = 0) +
#   scale_x_discrete(expand = c(0,0)) +
#   scale_y_discrete(expand = c(0,0)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.line = element_blank(),
#         axis.title = element_blank()) +
#   labs(title = paste("Correlation Matrix for VME Group:", vmeoi), fill = "Correlation")

# ggsave(plot_cor_selvars,
#   filename = paste0("output/02_StepByStepOutputs/", 
#     paste(vmeoi, poi, sspoi, "plot_cor_SelCMIPVars", subsamp_title, ".png", sep = "_")), 
#   width = 8, height = 6)

# Create table of variable selection results for this iteration ----
selected_cmip_vars <- selected_vme_vars[selected_vme_vars %in% names(cmip_layers)]
var_select_df[var_select_df$vmeoi == vmeoi & var_select_df$poi == poi & var_select_df$sspoi == sspoi, 
        selected_cmip_vars] <- 1
var_select_df[var_select_df$vmeoi == vmeoi & var_select_df$poi == poi & var_select_df$sspoi == sspoi, 
        setdiff(names(cmip_layers), selected_cmip_vars)] <- 0


# Prepare variable layer rasters for spatial predictions ----
vme_layers <- c(bathy_layers[vme_terrain_vars],
                compact(cmip_layers[selected_vme_vars]))
vme_layers_rast <- terra::rast(vme_layers)


# Model building ----

cat("Building RF model and extracting performance metrics across folds...\n")
# Create the 10 folds
set.seed(loop_seed)
folds <- caret::createFolds(vme_df$VME_P_A, k = 10, returnTrain = TRUE)

# Initialise storage for fold metrics
fold_metrics_list <- list()
fold_predictions_spatial <- list()
fold_predictions_spatial_reclass <- list()
fold_var_imp <- list()
fold_partialdep <- list()

for (i in 1:10) {
  
  # Get fold indices
  rf_train_idx <- folds[[i]]
  rf_test_idx <- setdiff(1:nrow(vme_df), rf_train_idx)
  
  # Train model on this fold
  cat("Training fold", i, "\n")  
  set.seed(loop_seed + i)
  rf_fold_model <- randomForest::randomForest(
    VME_P_A ~ .,
    data = vme_df[rf_train_idx, ],
    ntree = 500,
    mtry = floor(sqrt(ncol(vme_df) - 1)),
    strata = vme_df$VME_P_A[rf_train_idx],
    replace = FALSE,
    importance = TRUE
  )
  
  cat("  Retrieving fold predictions and metrics\n")
  # Get predictions on held-out test data
  rf_test_pred_prob <- predict(rf_fold_model, 
                            newdata = vme_df[rf_test_idx, ], 
                            type = "prob")
  
  # Extract probability of positive class (assuming second level)
  lev <- levels(vme_df$VME_P_A)
  pred_prob <- rf_test_pred_prob[, lev[2]]
  obs_numeric <- ifelse(vme_df$VME_P_A[rf_test_idx] == lev[2], 1, 0)
  
  # Calculate optimal threshold using Sens=Spec method
  threshold_df <- data.frame(
    id = 1:length(rf_test_idx),
    PA = ifelse(vme_df$VME_P_A[rf_test_idx] == lev[2], 1, 0),
    predprob = rf_test_pred_prob[, lev[2]]
  )
  
  opttsh <- PresenceAbsence::optimal.thresholds(
    threshold_df,
    opt.methods = "Sens=Spec"  # if wanted to also test prevalence: c("Sens=Spec", "ObsPrev")
  ) %>%
    pull(predprob)
  
  # Apply optimal threshold
  optimal_pred <- ifelse(pred_prob >= opttsh, lev[2], lev[1])
  optimal_pred <- factor(optimal_pred, levels = lev)
  obs_factor <- vme_df$VME_P_A[rf_test_idx]
  
  # Calculate confusion matrix
  cm <- caret::confusionMatrix(optimal_pred, obs_factor, positive = lev[2])
  
  # Extract and store metrics
  metrics <- unlist(c(cm$overall, cm$byClass))
  metrics["TSS"] <- cm$byClass["Sensitivity"] + cm$byClass["Specificity"] - 1
  metrics["OptThreshold"] <- opttsh
  metrics["Fold"] <- i
  
  fold_metrics_list[[i]] <- metrics

  # Fold variable importance
  cat("  Extracting variable importance\n")
  fold_var_imp[[i]] <- as.data.frame(randomForest::importance(rf_fold_model)) %>%
    rownames_to_column(var = "Variable") %>%
    arrange(desc(MeanDecreaseGini))

  # Fold partial dependence data
  # cat("  Extracting partial dependence data\n")
  # fold_partialdep[[i]] <- lapply(selected_vme_vars, function(var) {
  #   pdp::partial(rf_fold_model, 
  #     pred.var = var,
  #     plot = FALSE) %>%
  #     mutate(Variable = colnames(.)[1]) %>%
  #     rename(value = var)
  # })
  
  cat("  Generating spatial predictions\n")
  # Spatial predictions for this fold
  fold_predictions_spatial[[i]] <- terra::predict(
    vme_layers_rast,
    rf_fold_model,
    type = 'prob',
    na.rm = TRUE,
    index = 1:2
  )
  
  # Convert predictions to presence/absence using optimal threshold
  fold_predictions_spatial_reclass[[i]] <- terra::classify(
    fold_predictions_spatial[[i]][[2]], 
    rcl = matrix(c(-Inf, opttsh, 0,
                    opttsh, Inf, 1), 
                  ncol = 3, byrow = TRUE)
    )  # %>%
  #   terra::as.factor()
}

cat("Modelling complete. Processing results...\n")

# Combine all fold metrics into a dataframe
fold_metrics_df_i <- do.call(rbind, lapply(fold_metrics_list, function(x) {
  data.frame(t(x))
})) %>%
  pivot_longer(cols = -"Fold", names_to = "metric", values_to = "value") %>%
  mutate(VME_Group = vmeoi, Period = poi, SSP = sspoi) %>%
  relocate(VME_Group, Period, SSP, Fold, metric, value)
# fold_metrics_df <- bind_rows(fold_metrics_df, fold_metrics_df_i)

# Summarise metrics across folds
fold_metrics_summary_df_i <- fold_metrics_df_i %>%
  group_by(VME_Group, Period, SSP, metric) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE),
            .groups = "drop")
fold_metrics_summary_df <- bind_rows(fold_metrics_summary_df, fold_metrics_summary_df_i)

# Create spatial predictions stack
rf_pred_foldstack <- terra::rast(fold_predictions_spatial_reclass)
terra::writeRaster(rf_pred_foldstack, filename = paste0("output/02_Modelling_Outputs/", vmeoi, "/",
  paste(vmeoi, poi, sspoi, "rf_spatial_predictions", sep = "_"), ".tif"), overwrite = TRUE)


# Extract variable importance for final selected variables ----
cat("Extracting RF model variable importance metrics...\n")

# fold_var_imp_df <- lapply(fold_var_imp, function(fold) {
#   fold %>%
#     filter(Variable %in% selected_vme_vars) %>%
#     mutate(var_imp_rank = row_number())
# }) %>%
#   bind_rows(.id = "Fold") %>%
#   mutate(Variable = fct_reorder(Variable, MeanDecreaseGini, .fun = mean)) %>%
#   ungroup()
# write.csv(fold_var_imp_df, file = paste0("output/02_Modelling_Outputs/",
#   paste(vmeoi, poi, sspoi, "table_rf_VarImp", sep = "_"), ".csv"), row.names = FALSE)

# ggplot(fold_var_imp_df, aes(y = Variable, x = MeanDecreaseGini)) +
#   geom_boxplot() +
#   theme_bw() +
#   labs(y = "Predictor Variable", x = "Mean Decrease in Gini Index")

# ggsave(filename = paste0("output/02_Modelling_Outputs/",
#   paste(vmeoi, poi, sspoi, "plot_rf_VarImp", sep = "_"), ".png"), 
#   width = 6, height = 4)


# Extract partial dependence plots for each variable ----
# fold_partial_df <- lapply(fold_partialdep, function(fold) {
#   bind_rows(fold)
# }) %>%
#   bind_rows(.id = "Fold") %>%
#   arrange(Fold, Variable, value) %>%
#   mutate(Fold = as.factor(as.numeric(Fold)),
#          Variable = factor(Variable, levels = rev(levels(fold_var_imp_df$Variable))))
# write.csv(fold_partial_df, file = paste0("output/02_Modelling_Outputs/",
#   paste(vmeoi, poi, sspoi, "table_rf_PartialDep", sep = "_"), ".csv"), row.names = FALSE)

# ggplot(fold_partial_df, aes(x = value, y = yhat, colour = Fold)) +
#   geom_line() +
#   facet_wrap(~ Variable, scales = "free_x") +
#   theme_bw() +
#   labs(x = "Predictor Value", y = "Partial Dependence")

# ggsave(filename = paste0("output/02_Modelling_Outputs/",
#    paste(vmeoi, poi, sspoi, "plot_rf_PartialDep", sep = "_"), ".png"),
#   width = 10, height = 8)


# Calculate spatial metrics across folds ----
## Most frequent class (0/1)
rf_res_MaxClass <- terra::modal(rf_pred_foldstack, freq = FALSE)

## Frequency of most frequent class (fraction of runs)
# rf_res_MaxClassF <- terra::modal(rf_pred_foldstack, freq = TRUE) / 10  # old method - doesn't work
rf_res_freq_count <- sum(rf_pred_foldstack == rf_res_MaxClass, na.rm = TRUE)
rf_res_MaxClassF <- rf_res_freq_count / 10

## Average probability of classes
rf_res_AvgProb <- Reduce("+", fold_predictions_spatial) / 10

## Average probability of maximum frequency class
rf_res_MaxClassAvgProb <- terra::selectRange(rf_res_AvgProb, rf_res_MaxClass + 1)

## Combined confidence metric
rf_res_CombConf <- rf_res_MaxClassF * rf_res_MaxClassAvgProb

## Number of models predicting presence
rf_res_CVSum <- terra::app(rf_pred_foldstack, sum, na.rm = TRUE)


# Save rasters for each metric ----
rm(rf_res_freq_count, rf_res_AvgProb)
lapply(ls(pattern = "rf_res"), function(res_name) {
  res_raster <- get(res_name)
  terra::writeRaster(res_raster, filename = paste0("output/02_Modelling_Outputs/", vmeoi, "/",
    paste(vmeoi, poi, sspoi, res_name, sep = "_"), ".tif"), overwrite = TRUE)
})

# ## Create stack of computed raster metrics across folds
# rf_pred_comp <- c(
#   MaxClass = terra::as.factor(rf_res_MaxClass),
#   MaxClassF = rf_res_MaxClassF,
#   AvgProb = rf_res_AvgProb,
#   MaxClassAvgProb = rf_res_MaxClassAvgProb,
#   CombConf = rf_res_CombConf,
#   CVSum = rf_res_CVSum
# )

# rm(list = ls(pattern = "rf_res_"))
