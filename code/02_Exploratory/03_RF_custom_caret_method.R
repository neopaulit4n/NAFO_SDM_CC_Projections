
# Custom function with ALL metrics
rf_metrics <- function(data, lev = NULL, model = NULL) {
  
  # Extract predicted probabilities and observed values
  # data$pred is the class prediction, data$obs is the observed
  # For binary classification, we need the probability of the positive class
  pred_prob <- data[, lev[2]]  # probability of positive class
  obs_numeric <- ifelse(data$obs == lev[2], 1, 0)
  
  # Create dataframe for optimal.thresholds function
  # Format: ID, observed (0/1), predicted probability
  threshold_data <- data.frame(
    id = 1:nrow(data),
    PA = obs_numeric,
    predprob = pred_prob
  )
  
  # Calculate optimal threshold using Sens=Spec method
  opttsh <- PresenceAbsence::optimal.thresholds(
    threshold_data[, c("id", "PA", "predprob")],
    opt.methods = 'Sens=Spec'
  ) %>%
    pull(predprob)
  
  # Apply optimal threshold to create new predictions
  optimal_pred <- ifelse(pred_prob >= opttsh, lev[2], lev[1])
  optimal_pred <- factor(optimal_pred, levels = lev)
  
  # Calculate confusion matrix using optimal threshold predictions
  cm <- caret::confusionMatrix(optimal_pred, data$obs, positive = lev[2])
  
  # Extract metrics
  metrics <- unlist(c(cm$overall, cm$byClass))
  metrics["TSS"] <- cm$byClass["Sensitivity"] + cm$byClass["Specificity"] - 1
  metrics["OptThreshold"] <- opttsh
  return(metrics)
}

train_control <- caret::trainControl(
  method = "cv",
  number = 10,
  index = folds,
  classProbs = TRUE,  # get probability predictions
  savePredictions = "all",  # save all CV predictions
  returnResamp = "all",
  summaryFunction = rf_metrics
)

# Train the model with cross-validation
set.seed(411)
rf_model_selvars <- caret::train(
  VME_P_A ~ .,
  data = vme_df,
  method = "rf",
  trControl = train_control,
  ntree = 500,
  strata = vme_df$VME_P_A,
  replace = FALSE,
  importance = TRUE,
  # Prevent tuning - use single mtry value - same one as original code
  tuneGrid = data.frame(mtry = floor(sqrt(ncol(vme_df) - 1)))
)

set.seed(411)
rf_model_allvars <- caret::train(
  VME_P_A ~ .,
  data = vme_df,
  method = "rf",
  trControl = train_control,
  ntree = 500,
  strata = vme_df$VME_P_A,  # ensures presence/absence ratio remains constant
  replace = FALSE,
  importance = TRUE,
  # Prevent tuning - use single mtry value - same one as original code
  tuneGrid = data.frame(mtry = floor(sqrt(ncol(vme_df) - 1)))
)

# Extract metrics
rf_metrics_results_selvars <- rf_model_selvars$results
rf_metrics_results_allvars <- rf_model_allvars$results

# Extract variable importance
rf_varimp <- randomForest::importance(rf_model$finalModel) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Variable") %>%
  arrange(desc(MeanDecreaseGini))


# Predict on CMIP + bathy layers for whole area ----
vme_layers <- c(bathy_layers[vme_terrain_vars],
                compact(cmip_layers[selected_vme_vars$selected_vars]))
vme_layers_rast <- terra::rast(vme_layers)

rf_pred <- terra::predict(vme_layers_rast, rf_model_selvars$finalModel, 
                          type = "prob", na.rm = TRUE, index = 1:2)

# Convert probability to presence/absence using optimal threshold

## Using sens=spec threshold
opt_threshold <- rf_metrics_results_selvars$OptThreshold
rf_pred_pa_ss <- terra::classify(rf_pred[[2]], 
                              # from-to-becomes
                              rcl = matrix(c(-Inf, opt_threshold, 0,
                                             opt_threshold, Inf, 1), 
                                            ncol = 3, byrow = TRUE)) %>%
  # Change to factor
  terra::as.factor()

terra::plot(rf_pred_pa_ss)

## Using prevalence as threshold
prev <- sum(vme_df$VME_P_A == "Presence") / nrow(vme_df)
rf_pred_pa_prev <- terra::classify(rf_pred[[2]],
                                   rcl = matrix(c(-Inf, prev, 0,
                                            prev, Inf, 1),
                                          ncol = 3, byrow = TRUE)) %>%
  # Change to factor
  terra::as.factor()


# Most frequent class
rf_pred_pa_mfc <- terra::modal(rf_pred_pa_ss, freq = FALSE)
terra::plot(rf_pred_pa_mfc)


rf_pred_pa_mfc_f <- terra::modal(rf_pred, freq = TRUE)




# Average probability of maximum frequency class across folds
# rf_pred_pa_mfc_prob <- rf_pred[[2]] * (rf_pred_pa_mfc_f[[1]] / 10)

terra::plot(rf_pred)

