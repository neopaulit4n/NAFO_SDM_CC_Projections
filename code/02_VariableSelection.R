# Variable selection

# Replicating previous SDM variable elimination approach
# 1) Run preliminary RF models with all variables
# 2) Determine variable importance
# 3) Remove correlated variables based on importance ranking
# 4) VIF remaining variables

# Prepare VME group dataframe ----
vme_terrain_vars <- filter(terrain_topvars, VME_Group == vmeoi) %>%
  pull(variable)
vme_vars <- c(vme_terrain_vars, names(cmip_layers))
vme_df <- filter(cmip_comb_df, VME_Group == vmeoi) %>%
  select(all_of(c("VME_P_A", vme_vars)))

# Preliminary RF model formula ----  
rf_prelim_form <- as.formula(paste("VME_P_A ~", paste(vme_vars, collapse = " + ")))

# Preliminary RF model run (save/load) ----
cat("Running preliminary RF model...\n")
set.seed(loop_seed)
rf_prelim <- randomForest::randomForest(formula = rf_prelim_form,
                                        data = vme_df,
                                        importance = TRUE)

# Extract variable importance metrics ----
cat("Extracting preliminary RF model variable importance metrics...\n")
rf_prelim_imp <- as.data.frame(randomForest::importance(rf_prelim)) %>%
  rownames_to_column(var = "Variable") %>%
  arrange(desc(MeanDecreaseGini))

# Show var imp plot for each VME group, order each by descending importance
plot_rf_prelim_var_imp <- ggplot(rf_prelim_imp, aes(x = reorder(Variable, MeanDecreaseGini), 
                          y = MeanDecreaseGini)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() +
  labs(x = "Predictor Variable", y = "Mean Decrease in Gini Index") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 6))

ggsave(plot_rf_prelim_var_imp, 
  filename = paste0(output_folder,"/",vmeoi, "_plot_rf_prelim_var_imp.jpg"), 
  width = 6, height = 4)

# Calculate variable correlations ----

# Calculate VIF for selected variables for each VME group
# Check each VME group's VIF values and see if any are > 10; identify VME groups that need to be re-evaluated
# VIF values need to be < 10
# If values are > 10, need to recompute Spearman correlation at lower threshold 
#   (increments of 0.05) and re-run variable selection until all variables achieve VIF < 10.

cat("Calculating variable correlations...\n")
cor_df <- cor(
  vme_df[, vme_vars, drop = FALSE], 
  method = "spearman",
  use = "complete.obs") %>%
  as.data.frame %>%
  pivot_longer(everything(), names_to = "var2", values_to = "cor") %>%
  mutate(var1 = rep(vme_vars, each = length(vme_vars))) %>%
  select(var1, var2, cor)
# write.csv(cor_df, file = paste0(output_folder,"/",vmeoi,"_table_cor_baseline_AllCMIPVars.csv"), row.names = FALSE)

# plot_cor_allvars <- ggplot(data = cor_df, aes(x = var1, y = var2, fill = cor)) +
#   geom_tile(colour = "black") +
#   geom_text(aes(label = ifelse(cor > 0.7 | cor < -0.7, "*", "")), size = 3) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", limit = c(-1,1), midpoint = 0) +
#   scale_x_discrete(expand = c(0,0)) +
#   scale_y_discrete(expand = c(0,0)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.text = element_text(size = 8),
#         axis.line = element_blank(),
#         axis.title = element_blank()) +
#   labs(fill = "Correlation")

# ggsave(plot_cor_allvars,
#   filename = paste0(output_folder,"/",vmeoi,"_plot_cor_AllCMIPVars_ForVarSelec.jpg"), 
#   width = 8, height = 6)

# Remove correlated variables based on importance ranking and VIF values ----
cat("Selecting variables...\n")
select_vme_vars <- function(pred_var_df, vmeoi, cor_threshold = 0.7, verbose = FALSE) {
  
  # cat("\n==============================\n")
  # cat("\nSelecting variables for VME group:", vmeoi, "\n")
  
  vif_threshold_met <- FALSE
  cor_threshold_new <- cor_threshold
  
  while (!vif_threshold_met) {
    
    # Calculate Spearman correlations and remove correlated variables ----
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

    # Determine if there are aliased variables to remove
    find_aliased_vars <- function() {
      data_alias <- as.data.frame(pred_var_df[, uncor_vars])
      
      form <- formula(paste("fooy ~ ", paste(strsplit(names(data_alias), " "), collapse = " + ")))
      data_alias <- data.frame(fooy = 1 + rnorm(nrow(data_alias)), data_alias)
      lm_mod <- lm(form, data_alias)

      # Find aliased variables
      aliased_vars <- alias(lm_mod)$Complete
      return(aliased_vars)
    }
    aliased_vars_res <- find_aliased_vars()

    # Remove aliased vars prior to VIF if aliased vars exist ----
    if (!is.null(aliased_vars_res)) {
      remove_aliased <- function() {

        # Find aliased variables
        aliased_vars <- aliased_vars_res %>% 
          as.data.frame() %>%        
          rownames_to_column() %>%
          rename(aliased_var = rowname) %>%
          pivot_longer(cols = !aliased_var, names_to = "variable", values_to = "alias_value") %>%
          filter(alias_value != 0)
        aliased_vars <- c(unique(aliased_vars$aliased_var), unique(aliased_vars$variable))

        # Of the aliased variables, eliminate the least important one based on prelim RF results for each aliased variable grouping
        rf_prelim_imp_aliased <- filter(rf_prelim_imp, Variable %in% aliased_vars) %>%
          # Extract which variable to group by
          mutate(var_group = str_extract(Variable, "^(\\w+?)_", group = 1)) %>%
          group_by(var_group) %>%
          slice(tail(row_number(), 1))
        aliased_vars_eliminate <- rf_prelim_imp_aliased$Variable
        uncor_unaliased_vars <- uncor_vars[!(uncor_vars %in% aliased_vars_eliminate)]

        return(uncor_unaliased_vars)
      }
    
      uncor_vars <- remove_aliased()  # remove_aliased(pred_var_df[, uncor_vars])    
    }


    corvif <- function(data, uncorrelated = uncor_vars, verbose = FALSE) {
      data <- as.data.frame(data)
      
      form <- formula(paste("fooy ~ ", paste(strsplit(names(data), " "), collapse = " + ")))
      data <- data.frame(fooy = 1 + rnorm(nrow(data)), data)
      lm_mod <- lm(form, data)
      
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
vif_df_i <- data.frame(
  vmeoi = vmeoi, 
  variable = vme_var_selection$selected_vars, 
  vif = vme_var_selection$vif_values$vif,
  final_cor_thresh = vme_var_selection$final_cor_threshold
)

# Re-append VME terrain variable(s) if not selected
if (!vme_terrain_vars %in% selected_vme_vars) {
  selected_vme_vars <- c(vme_terrain_vars, selected_vme_vars)
  vif_df_i <- rbind(
    vif_df_i, 
    c(
      vmeoi = vmeoi,
      variable = vme_terrain_vars,
      vif = NA,
      final_cor_thresh = NA
    )
  ) %>%
    mutate(across(c(vif, final_cor_thresh), as.numeric))
}

vif_df <- rbind(vif_df, vif_df_i)

# Create final dataframe with selected variables ----
vme_df <- vme_df %>%
  select(all_of(c("VME_P_A", selected_vme_vars)))

# Create table of variable selection results for this iteration ----
selected_cmip_vars <- selected_vme_vars[selected_vme_vars %in% names(cmip_layers)]
var_select_df[var_select_df$vmeoi == vmeoi, selected_cmip_vars] <- 1
var_select_df[var_select_df$vmeoi == vmeoi, setdiff(names(cmip_layers), selected_cmip_vars)] <- 0
# selected_vme_vars <- c(vme_terrain_vars, cmip_vars)  # test without variable selection

# Prepare variable layer rasters for spatial predictions ----
cat("Preparing variable raster layers for spatial predictions...\n")
vme_layers_baseline <- c(bathy_layers[vme_terrain_vars], compact(cmip_layers[selected_vme_vars])) %>%
  terra::rast(.)

vme_layers_proj <- lapply(period_all, function(poi) {
  lapply(ssp_all, function(sspoi) {
    c(bathy_layers[vme_terrain_vars], compact(cmip_layers_proj[[poi]][[sspoi]][selected_vme_vars])) %>%
      terra::rast(.)     
  }) %>%
    set_names(ssp_all)
}) %>%
  set_names(period_all)