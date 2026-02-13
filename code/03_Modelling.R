
# Starting models

# Run one of models in full with the new variables for Period 1 and one of the SSPs 
# Then compare it to the BNAM results


# Load data ----
source("code/01_LoadData.R")

# Average predictors for each period and SSP
cmip_df_period_ssp <- cmip_df %>%
  filter(!is.na(period)) %>%
  group_by(lon, lat, period, ssp) %>%
  # summarise(across(sobavg:mldavg, \(x) mean(x, na.rm = TRUE)), .groups = "drop")  
  summarise(across(sobavg:mldavg, list(mean = ~mean(.x, na.rm = TRUE),
                                       min = ~min(.x, na.rm = TRUE),
                                       max = ~max(.x, na.rm = TRUE)),
                   .names = "{col}_{fn}")) %>%
  ungroup() %>%
  # Only keep max mldavg
  select(-(starts_with("mldavg") & !ends_with("_max")))


# Get study area extent and create spatial mask ----
sa <- terra::rast("data/raw/BNAM_Data_From_Cam/BNAM_From_NAFO_SharePoint/NRA_BNAM_b_cur_avg_max.tif") %>%
  terra::as.polygons(.) %>%
  sf::st_as_sf(.) %>%
  sf::st_transform(4326) %>%
  mutate(NRA_BNAM_b_tmp_mean = 1) %>%
  group_by(NRA_BNAM_b_tmp_mean) %>%
  summarise(geometry = sf::st_union(geometry))

sa_lims <- sf::st_coordinates(sa) %>% as.data.frame %>%
  select(X,Y)
sa_lims <- c(min(sa_lims$X), max(sa_lims$X),min(sa_lims$Y),max(sa_lims$Y))

# Transform CMIP data to raster
transform_cmip_to_raster <- function(data, varstat, period, ssp) {
  
  # Data must be ensembled dataframe with period and ssp columns
  
  # Select a layer of data to raster
  df <- data %>%
    filter(ssp == ssp, period == period) %>%
    select(lon, lat, !!sym(varstat)) %>%
    sf::st_as_sf(coords = c("lon","lat"), crs = 4326)
  
  # Transform into terra vector of points
  pts <- terra::vect(df)
  
  # Determine resolution of data in degrees
  # res <- round(ens_df_period_cell$lon[2]-ens_df_period_cell$lon[1], 5)
  res <- round(sort(unique(data$lon))[2] - sort(unique(data$lon))[1], 5)
  
  # Create template raster
  rast_template <- terra::rast(
    xmin = sa_lims[1], xmax = sa_lims[2], ymin = sa_lims[3], ymax = sa_lims[4],
    resolution = res,
    crs = "EPSG:4326"
  )
  # rast_template <- bathy_layers[[1]]
  
  # Rasterise points to grid
  rast_result <- terra::rasterize(pts, rast_template, field = varstat)
  
  # Resample to match BNAM raster
  rast_resamp <- terra::resample(rast_result, bnam_layers[[1]], method = "bilinear")
  return(rast_resamp)
  
}

vars <- colnames(select(cmip_df_period_ssp, contains("_")))

cmip_layers <- lapply(vars, transform_cmip_to_raster, 
                      data = cmip_df_period_ssp, 
                      period = "P1", 
                      ssp = "2-4.5")
names(cmip_layers) <- vars
# terra::plot(cmip_layers[[1]])

cmip_pred_df <- lapply(c(bathy_layers, cmip_layers), function(layer) {
  terra::extract(layer, 
                 select(resp_df, Start_Long_DD, Start_Lat_DD)) %>%
    select(-ID)
}) %>%
  bind_cols() %>%
  set_names(c(names(bathy_layers), names(cmip_layers)))


# Combine predictor and response dataframes ----
cmip_comb_df <- bind_cols(resp_df, cmip_pred_df) %>%
  mutate(VME_P_A = factor(VME_P_A, levels = c(0, 1), labels = c("Absence", "Presence"))) %>%
  drop_na()


# Replicating previous SDM variable elimination approach ----
# 1: Run preliminary RF models with all variables
# 2: Determine variable importance
# 3: Remove correlated variables based on importance ranking
# 4: VIF remaining variables


## Preliminary RF model test ----
vme_group <- "black_corals"
vme_terrain_vars <- filter(terrain_topvars, VME_Group == vme_group) %>%
  pull(variable)
vme_vars <- c(vme_terrain_vars, names(cmip_layers))
vme_df <- filter(cmip_comb_df, VME_Group == vme_group) %>%
  select(all_of(c("VME_P_A", vme_vars)))

rf_prelim_form <- as.formula(paste("VME_P_A ~", 
                                   paste(vme_vars, collapse = " + ")))


## Preliminary RF model run (save/load) ----

# Run through all VME groups
# rf_prelim <- lapply(unique(cmip_comb_df$VME_Group), function(vme_group) {
#   set.seed(804)
#   rf_model <- randomForest::randomForest(formula = rf_prelim_form,
#                                          data = filter(cmip_comb_df, VME_Group == vme_group),
#                                          importance=TRUE)
#   return(rf_model)
# }) %>%
#   set_names(unique(cmip_comb_df$VME_Group))

# Single VME group version
# set.seed(804)
# rf_prelim <- randomForest::randomForest(formula = rf_prelim_form,
#                                        data = vme_df,
#                                        importance=TRUE)
# 
# saveRDS(rf_model, file = paste0("data/processed/rf_prelim_",vme_group,".rds"))
rf_prelim <- readRDS(paste0("data/processed/rf_prelim_",vme_group,".rds"))

# sapply(rf_prelim, print)


## Extract variable importance metrics ----
# rf_prelim_imp <- lapply(rf_prelim, function(model) {
#   imp_df <- as.data.frame(randomForest::importance(model)) %>%
#     rownames_to_column(var = "Variable") %>%
#     arrange(desc(MeanDecreaseGini))
#   return(imp_df)
# }) %>%
#   bind_rows(.id = "VME_Group")

rf_prelim_imp <- as.data.frame(randomForest::importance(rf_prelim)) %>%
  rownames_to_column(var = "Variable") %>%
  arrange(desc(MeanDecreaseGini))

# Show var imp plot for each VME group, order each by descending importance
ggplot(rf_prelim_imp, aes(x = reorder(Variable, MeanDecreaseGini), 
                          y = MeanDecreaseGini)) +
                          # fill = VME_Group)) +
  geom_bar(stat = "identity") +
  # facet_wrap(~ VME_Group, scales = "free_y") +
  coord_flip() +
  theme_bw() +
  labs(x = "Predictor Variable", y = "Mean Decrease in Gini Index") +
  theme(legend.position = "none")

# randomForest::varImpPlot(rf_prelim)


## Plot partial dependence plots ----
# randomForest::partialPlot(x = rf_prelim[[1]], 
#                           pred.data = as.data.frame(filter(comb_df_compl, VME_Group == "black_corals")),
#                           x.var = sosavg_max)


## Remove correlated variables ----

# Calculate VIF for selected variables for each VME group
# Check each VME group's VIF values and see if any are > 10; identify VME groups that need to be re-evaluated
# VIF values need to be < 10
# If values are > 10, need to recompute Spearman correlation at lower threshold 
#   (increments of 0.05) and re-run variable selection until all variables achieve VIF < 10.

select_vme_vars <- function(pred_var_df, vme_group, cor_threshold = 0.7, verbose = FALSE) {
  
  cat("\n==============================\n")
  cat("\nSelecting variables for VME group:", vme_group, "\n")
  
  vif_threshold_met <- FALSE
  cor_threshold_new <- cor_threshold
  
  while (!vif_threshold_met) {
    
    # Calculate Spearman correlations
    remove_cor_variables <- function(data, priority_list, current_threshold = 0.7, verbose = FALSE) {
      
      vars <- colnames(data)
      
      # Calculate correlation
      cor_matrix <- matrix(0, nrow = ncol(data), ncol = ncol(data))
      rownames(cor_matrix) <- colnames(cor_matrix) <- colnames(data)
      
      numeric_cor <- cor(data[, vars, drop = FALSE], 
                         method = "spearman",
                         use = "pairwise.complete.obs")
      cor_matrix[vars, vars] <- numeric_cor
      
      # Prioritise variables
      priority_dict <- setNames(seq_along(priority_list), priority_list)
      all_vars <- colnames(data)
      
      # Assign priorities (variables in priority list come first, by their order)
      var_priority <- sapply(all_vars, function(v) {
        if (v %in% priority_list) {
          return(priority_dict[v])
        } else {
          return(length(priority_list) + 1) # Lower priority for non-listed vars
        }
      })
      
      # Sort variables by priority
      sorted_vars <- all_vars[order(var_priority)]
      
      # Identify highly correlated pairs and remove lower priority variables
      removed_vars <- c()
      
      for (i in 1:length(sorted_vars)) {
        var_i <- sorted_vars[i]
        if (verbose) cat("\nEvaluating variable:", var_i, "\n")
        
        # Skip if this variable was already removed
        if (var_i %in% removed_vars) {
          if (verbose) cat("  Skipped (already removed)\n")
          next
        }
        
        for (j in i:length(sorted_vars)) {
          var_j <- sorted_vars[j]
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
    
    cat("\n   Removing correlated variables (threshold:", cor_threshold_new, ")\n")
    
    uncor_vars <- remove_cor_variables(pred_var_df, 
                                       priority_list = rf_prelim_imp %>%
                                         # filter(VME_Group == vme_group) %>%
                                         pull(Variable),
                                       current_threshold = cor_threshold_new,
                                       verbose = verbose)
    
    # Calculate VIF values
    cat("\n   Calculating VIF values for selected variables\n")
    
    corvif <- function(data, verbose = TRUE) {
      data <- as.data.frame(data)
      
      form    <- formula(paste("fooy ~ ", paste(strsplit(names(data), " "), collapse = " + ")))
      data   <- data.frame(fooy = 1 + rnorm(nrow(data)), data)
      lm_mod  <- lm(form, data)
      
      if (verbose) print(data.frame(vif=car::vif(lm_mod)))
      return(data.frame(vif=car::vif(lm_mod)))
    }
    vif_values <- corvif(pred_var_df[, uncor_vars], verbose = verbose)
    
    # Check if all VIF values are < 10
    if (all(vif_values$vif < 10)) {
      cat("\n     All VIF values are below threshold for VME group:", vme_group, "\n")
      vif_threshold_met <- TRUE
      # uncor_vars[[vme_group]] <- selected_vars_new
      # vars_vif[[vme_group]] <- vif_values_new
      return(list(selected_vars = uncor_vars,
                  vif_values = vif_values,
                  final_cor_threshold = cor_threshold_new))
    }
    else {
      cat("\n     VIF values exceed threshold (10); lowering correlation threshold and re-evaluating.\n")
      cor_threshold_new <- cor_threshold_new - 0.05
    }
  }
}

# selected_vme_vars <- lapply(unique(comb_df_compl$VME_Group), select_vme_vars, cor_threshold = 0.7) %>%
#   set_names(unique(comb_df_compl$VME_Group))

selected_vme_vars <- select_vme_vars(vme_df[,-1], vme_group = vme_group, cor_threshold = 0.7, verbose = TRUE)

# Create final dataframe with selected variables
vme_df_sel <- vme_df %>%
  select(all_of(c("VME_P_A", selected_vme_vars$selected_vars)))


# Model building ----

# Create the 10 folds
set.seed(411)
folds <- caret::createFolds(vme_df$VME_P_A, k = 10, returnTrain = TRUE)

# Initialize storage for fold metrics
fold_metrics_list <- list()
fold_predictions_spatial <- list()
fold_predictions_spatial_reclass <- list()

# Prepare variable layer rasters for spatial predictions
vme_layers <- c(bathy_layers[vme_terrain_vars],
                compact(cmip_layers[selected_vme_vars$selected_vars]))
vme_layers_rast <- terra::rast(vme_layers)

for (i in 1:10) {
  cat("Training fold", i, "\n")
  
  # Get fold indices
  rf_train_idx <- folds[[i]]
  rf_test_idx <- setdiff(1:nrow(vme_df_sel), rf_train_idx)
  
  # Train model on this fold
  set.seed(411)
  rf_fold_model <- randomForest::randomForest(
    VME_P_A ~ .,
    data = vme_df_sel[rf_train_idx, ],
    ntree = 500,
    mtry = floor(sqrt(ncol(vme_df_sel) - 1)),
    strata = vme_df_sel$VME_P_A[rf_train_idx],
    replace = FALSE,
    importance = TRUE
  )
  
  # Get predictions on held-out test data
  rf_test_pred_prob <- predict(rf_fold_model, 
                            newdata = vme_df_sel[rf_test_idx, ], 
                            type = "prob")
  
  # Extract probability of positive class (assuming second level)
  lev <- levels(vme_df_sel$VME_P_A)
  pred_prob <- rf_test_pred_prob[, lev[2]]
  obs_numeric <- ifelse(vme_df_sel$VME_P_A[rf_test_idx] == lev[2], 1, 0)
  
  # Calculate optimal threshold using Sens=Spec method
  threshold_df <- data.frame(
    id = 1:length(rf_test_idx),
    PA = ifelse(vme_df_sel$VME_P_A[rf_test_idx] == lev[2], 1, 0),
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
  obs_factor <- vme_df_sel$VME_P_A[rf_test_idx]
  
  # Calculate confusion matrix
  cm <- caret::confusionMatrix(optimal_pred, obs_factor, positive = lev[2])
  
  # Extract and store metrics
  metrics <- unlist(c(cm$overall, cm$byClass))
  metrics["TSS"] <- cm$byClass["Sensitivity"] + cm$byClass["Specificity"] - 1
  metrics["OptThreshold"] <- opttsh
  metrics["Fold"] <- i
  
  fold_metrics_list[[i]] <- metrics
  
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

# Combine all fold metrics into a dataframe
fold_metrics_df <- do.call(rbind, lapply(fold_metrics_list, function(x) {
  data.frame(t(x))
})) %>%
  pivot_longer(cols = -"Fold", names_to = "metric", values_to = "value")

fold_metrics_summary_df <- fold_metrics_df %>%
  group_by(metric) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE),
            .groups = "drop")


# Create spatial predictions stack
rf_pred_stack <- terra::rast(fold_predictions_spatial_reclass)
# ggplot() +
#   theme_classic() +
#   tidyterra::geom_spatraster(data = rf_pred_stack[[1]], na.rm = TRUE) +
#   labs(title = paste("Predicted Presence/Absence for Fold", 1),
#        x = "Longitude", y = "Latitude")

# Calculate spatial metrics across folds ----
## Most frequent class (0/1)
rf_res_MaxClass <- terra::modal(rf_pred_stack, freq = FALSE)

## Frequency of most frequent class (fraction of runs)
# rf_res_MaxClassF <- terra::modal(rf_pred_stack, freq = TRUE) / 10  # old method - doesn't work
rf_res_freq_count <- sum(rf_pred_stack == rf_res_MaxClass, na.rm = TRUE)
rf_res_MaxClassF <- rf_res_freq_count / 10

## Average probability of classes
# rf_res_AvgProb <- terra::app(terra::rast(fold_predictions_spatial), mean, na.rm = TRUE)
rf_res_AvgProb <- Reduce("+", fold_predictions_spatial) / 10

## Average probability of maximum frequency class
rf_res_MaxClassAvgProb <- terra::selectRange(rf_res_AvgProb, rf_res_MaxClass + 1)

## Combined confidence metric
rf_res_CombConf <- rf_res_MaxClassF * rf_res_MaxClassAvgProb

## Number of models predicting presence
rf_res_NumPres <- terra::app(rf_pred_stack, sum, na.rm = TRUE)

## Create stack of computed raster metrics across folds
rf_pred_comp <- c(
  MaxClass = terra::as.factor(rf_res_MaxClass),
  MaxClassF = rf_res_MaxClassF,
  AvgProb = rf_res_AvgProb,
  MaxClassAvgProb = rf_res_MaxClassAvgProb,
  CombConf = rf_res_CombConf,
  NumPres = rf_res_NumPres
)

rm(list = ls(pattern = "rf_res_"))


# # Custom function with ALL metrics
# rf_metrics <- function(data, lev = NULL, model = NULL) {
#   
#   # Extract predicted probabilities and observed values
#   # data$pred is the class prediction, data$obs is the observed
#   # For binary classification, we need the probability of the positive class
#   pred_prob <- data[, lev[2]]  # probability of positive class
#   obs_numeric <- ifelse(data$obs == lev[2], 1, 0)
#   
#   # Create dataframe for optimal.thresholds function
#   # Format: ID, observed (0/1), predicted probability
#   threshold_data <- data.frame(
#     id = 1:nrow(data),
#     PA = obs_numeric,
#     predprob = pred_prob
#   )
#   
#   # Calculate optimal threshold using Sens=Spec method
#   opttsh <- PresenceAbsence::optimal.thresholds(
#     threshold_data[, c("id", "PA", "predprob")],
#     opt.methods = 'Sens=Spec'
#   ) %>%
#     pull(predprob)
#   
#   # Apply optimal threshold to create new predictions
#   optimal_pred <- ifelse(pred_prob >= opttsh, lev[2], lev[1])
#   optimal_pred <- factor(optimal_pred, levels = lev)
#   
#   # Calculate confusion matrix using optimal threshold predictions
#   cm <- caret::confusionMatrix(optimal_pred, data$obs, positive = lev[2])
#   
#   # Extract metrics
#   metrics <- unlist(c(cm$overall, cm$byClass))
#   metrics["TSS"] <- cm$byClass["Sensitivity"] + cm$byClass["Specificity"] - 1
#   metrics["OptThreshold"] <- opttsh
#   return(metrics)
# }
# 
# train_control <- caret::trainControl(
#   method = "cv",
#   number = 10,
#   index = folds,
#   classProbs = TRUE,  # get probability predictions
#   savePredictions = "all",  # save all CV predictions
#   returnResamp = "all",
#   summaryFunction = rf_metrics
# )
# 
# # Train the model with cross-validation
# set.seed(411)
# rf_model_selvars <- caret::train(
#   VME_P_A ~ .,
#   data = vme_df_sel,
#   method = "rf",
#   trControl = train_control,
#   ntree = 500,
#   strata = vme_df$VME_P_A,
#   replace = FALSE,
#   importance = TRUE,
#   # Prevent tuning - use single mtry value - same one as original code
#   tuneGrid = data.frame(mtry = floor(sqrt(ncol(vme_df) - 1)))
# )
# 
# set.seed(411)
# rf_model_allvars <- caret::train(
#   VME_P_A ~ .,
#   data = vme_df,
#   method = "rf",
#   trControl = train_control,
#   ntree = 500,
#   strata = vme_df$VME_P_A,  # ensures presence/absence ratio remains constant
#   replace = FALSE,
#   importance = TRUE,
#   # Prevent tuning - use single mtry value - same one as original code
#   tuneGrid = data.frame(mtry = floor(sqrt(ncol(vme_df) - 1)))
# )
# 
# # Extract metrics
# rf_metrics_results_selvars <- rf_model_selvars$results
# rf_metrics_results_allvars <- rf_model_allvars$results
# 
# # Extract variable importance
# rf_varimp <- randomForest::importance(rf_model$finalModel) %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "Variable") %>%
#   arrange(desc(MeanDecreaseGini))
# 
# 
# # Predict on CMIP + bathy layers for whole area ----
# vme_layers <- c(bathy_layers[vme_terrain_vars],
#                 compact(cmip_layers[selected_vme_vars$selected_vars]))
# vme_layers_rast <- terra::rast(vme_layers)
# 
# rf_pred <- terra::predict(vme_layers_rast, rf_model_selvars$finalModel, 
#                           type = "prob", na.rm = TRUE, index = 1:2)
# 
# # Convert probability to presence/absence using optimal threshold
# 
# ## Using sens=spec threshold
# opt_threshold <- rf_metrics_results_selvars$OptThreshold
# rf_pred_pa_ss <- terra::classify(rf_pred[[2]], 
#                               # from-to-becomes
#                               rcl = matrix(c(-Inf, opt_threshold, 0,
#                                              opt_threshold, Inf, 1), 
#                                             ncol = 3, byrow = TRUE)) %>%
#   # Change to factor
#   terra::as.factor()
# 
# terra::plot(rf_pred_pa_ss)
# 
# ## Using prevalence as threshold
# prev <- sum(vme_df_sel$VME_P_A == "Presence") / nrow(vme_df_sel)
# rf_pred_pa_prev <- terra::classify(rf_pred[[2]],
#                                    rcl = matrix(c(-Inf, prev, 0,
#                                             prev, Inf, 1),
#                                           ncol = 3, byrow = TRUE)) %>%
#   # Change to factor
#   terra::as.factor()
# 
# 
# # Most frequent class
# rf_pred_pa_mfc <- terra::modal(rf_pred_pa_ss, freq = FALSE)
# terra::plot(rf_pred_pa_mfc)
# 
# 
# rf_pred_pa_mfc_f <- terra::modal(rf_pred, freq = TRUE)
# 
# 
# 
# 
# # Average probability of maximum frequency class across folds
# # rf_pred_pa_mfc_prob <- rf_pred[[2]] * (rf_pred_pa_mfc_f[[1]] / 10)
# 
# terra::plot(rf_pred)

