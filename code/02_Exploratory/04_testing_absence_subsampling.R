
# Testing absence subsampling for preliminary RF model used to determine variable selection priorities
# Should be called from code/03_Modelling.R if testing_subsampling == TRUE

# Subsample absences to match the number of presences
vme_n_pres <- sum(vme_df$VME_P_A == "Presence")
# set.seed(799)
# vme_df_subsamp <- vme_df %>%
#   group_by(VME_P_A) %>%
#   slice_sample(n = vme_n_pres) %>%
#   ungroup()

## Preliminary RF model formula ----
rf_prelim_form <- as.formula(paste("VME_P_A ~", 
                                   paste(vme_vars, collapse = " + ")))

# Loop over 10 random seeds to see how variable importance changes between folds
rf_prelim_list <- lapply(1:10, function(seed) {
  cat("Running preliminary RF model for seed", seed, "\n")
  lapply(c("subsamp","full"), function(df_type) {
    if (df_type == "full") {
      df <- vme_df 
    } else {
      # Subsample absences to match the number of presences
      set.seed(seed)
      df <- vme_df %>% group_by(VME_P_A) %>% slice_sample(n = vme_n_pres) %>% ungroup()
    }
    set.seed(seed)
    randomForest::randomForest(formula = rf_prelim_form,
                               data = df,
                               importance = TRUE)
  }) %>%
    set_names(c("Subsampled", "Full"))
})

# Extract variable importance lists for each seed and dataset type
rf_prelim_imp <- lapply(rf_prelim_list, function(seed) {
  lapply(seed, function(model) {
    as.data.frame(randomForest::importance(model)) %>%
       rownames_to_column(var = "Variable") %>%
      arrange(desc(MeanDecreaseGini))
  }) %>%
     bind_rows(.id = "dataset_type")
}) %>%
  bind_rows(.id = "Seed") %>%
  arrange(Seed, dataset_type, desc(MeanDecreaseGini)) %>%
  group_by(Seed, dataset_type) %>%
  mutate(var_imp_rank = row_number()) %>%
  ungroup()

# Plot summary variable importance rank over all seeds
plot_rf_prelim_subsamptest_var_imp_rank <- ggplot(data = rf_prelim_imp, 
  aes(x = dataset_type, y = var_imp_rank, fill = dataset_type)) +
  facet_wrap(~ Variable) +
  theme_bw() +
  geom_boxplot() +
  ggtitle("Variable importance rank across 10 random seeds for preliminary RF model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none")
ggsave(plot_rf_prelim_subsamptest_var_imp_rank, 
  filename = paste0("output/02_StepByStepOutputs/",
  vmeoi, "_", poi, "_", sspoi,
  "_plot_subsamptest_rf_prelim_var_imp_rank.png"), 
  width = 8, height = 6)

# Extract RF metrics for each seed and dataset type
rf_prelim_metrics <- lapply(rf_prelim_list, function(seed) {
  lapply(seed, function(model) {
    cm <- caret::confusionMatrix(model$predicted, model$y)
    cm_metrics <- c(cm$overall, cm$byClass)
    cm_metrics["TSS"] <- cm$byClass["Sensitivity"] + cm$byClass["Specificity"] - 1 
    cm_metrics <- data.frame(value = cm_metrics) %>%
      rownames_to_column(var = "Metric")
    return(cm_metrics)
  }) %>%
    bind_rows(.id = "dataset_type")
}) %>%
  bind_rows(.id = "Seed")

# Plot summary metrics over all seeds
plot_rf_prelim_subsamptest_metrics <- ggplot(data = rf_prelim_metrics, 
  aes(x = dataset_type, y = value, fill = dataset_type)) +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_bw() +
  geom_boxplot() +
  ggtitle("Model performance metrics across 10 random seeds for preliminary RF model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none",
        strip.text = element_text(size = 7))
ggsave(plot_rf_prelim_subsamptest_metrics, 
  filename = paste0("output/02_StepByStepOutputs/",
  vmeoi, "_", poi, "_", sspoi,
  "_plot_subsamptest_rf_prelim_metrics.png"), 
  width = 8, height = 6)
