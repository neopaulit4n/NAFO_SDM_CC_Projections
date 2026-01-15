
# Testing variable elimination methods

# Method 1: Replicating previous SDM variable elimination approach ----
# 1: Run preliminiary RF models with all variables
# 2: Determine variable importance
# 3: Remove correlated variables based on importance ranking
# 4: VIF remaining variables

## Preliminary RF model test ----
rf_prelim_form <- as.formula(paste("VME_P_A ~", 
                                 paste(colnames(pred_values), collapse = " + ")))

rf_prelim <- lapply(unique(comb_df_compl$VME_Group), function(vme_group) {
  set.seed(804)
  rf_model <- randomForest::randomForest(formula = rf_test_form,
                                         data = filter(comb_df_compl, VME_Group == vme_group),
                                         importance=TRUE)
  return(rf_model)
}) %>%
  set_names(unique(comb_df_compl$VME_Group))

# saveRDS(rf_test, file = "data/processed/prelim_rf_test.rds")
rf_prelim <- readRDS("data/processed/prelim_rf_test.rds")

sapply(rf_prelim, print)


## Extract variable importance metrics ----
rf_prelim_imp <- lapply(rf_prelim, function(model) {
  imp_df <- as.data.frame(randomForest::importance(model)) %>%
    rownames_to_column(var = "Variable") %>%
    arrange(desc(MeanDecreaseGini))
  return(imp_df)
}) %>%
  bind_rows(.id = "VME_Group")

randomForest::varImpPlot(rf_prelim[[1]])

## Plot partial dependence plots ----
randomForest::partialPlot(x = rf_prelim[[1]], 
                          pred.data = as.data.frame(filter(comb_df_compl, VME_Group == "black_corals")),
                          x.var = sosavg_max)


## Remove correlated variables ----

# Function to select variables
remove_cor_variables <- function(data, priority_list, cor_threshold = 0.7, verbose = TRUE) {
  
  vars <- colnames(data)
  
  # Calculate correlation
  cor_matrix <- matrix(0, nrow = ncol(data), ncol = ncol(data))
  rownames(cor_matrix) <- colnames(cor_matrix) <- colnames(data)
  
  numeric_cor <- cor(data[, vars, drop = FALSE], 
                     method = "spearman")
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
      
      # Apply appropriate threshold based on variable types
      current_threshold <- cor_threshold
      
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

# Variable selection by VME group
cor_threshold <- 0.7
uncor_vars <- lapply(unique(comb_df_compl$VME_Group), function(vme_group) {
  
  # Set priority for removal in order of RF preliminary variable importance
  priority_list <- rf_prelim_imp %>%
    filter(VME_Group == vme_group) %>%
    pull(Variable)
  
  # Select variables
  selected_vars <- remove_cor_variables(pred_values, priority_list, cor_threshold = 0.7, verbose = FALSE)
  
  # Sort variables by importance
  rf_prelim_imp %>%
    filter(VME_Group == vme_group, Variable %in% selected_vars) %>%
    arrange(desc(MeanDecreaseGini)) %>%
    pull(Variable)
  
}) %>%
  set_names(unique(comb_df_compl$VME_Group))


## VIF ----

# VIF function
corvif <- function(data, verbose = TRUE) {
  data <- as.data.frame(data)
  
  form    <- formula(paste("fooy ~ ", paste(strsplit(names(data), " "), collapse = " + ")))
  data   <- data.frame(fooy = 1 + rnorm(nrow(data)), data)
  lm_mod  <- lm(form, data)
  
  if (verbose) print(data.frame(vif=car::vif(lm_mod)))
}

# Calculate VIF for selected variables for each VME group
vars_vif <- lapply(names(uncor_vars), function(vme_group) {
  cat("\nVIF: ", vme_group, "\n")
  vif_values <- corvif(pred_values[, uncor_vars[[vme_group]]])
  return(vif_values)
}) %>%
  set_names(names(uncor_vars))

# Check each VME group's VIF values and see if any are > 10; identify VME groups that need to be re-evaluated
# VIF values need to be < 10
# If values are > 10, need to recompute Spearman correlation at lower threshold 
#   (increments of 0.05) and re-run variable selection until all variables achieve VIF < 10.
vme_reval <- names(vars_vif)[sapply(vars_vif, function(vif_df) any(vif_df$vif > 10))]

# Re-evaluate VME groups with VIF > 10
select_vme_vars <- function(vme_group, cor_threshold = 0.7) {
  
  vif_threshold_met <- FALSE
  cor_threshold_new <- cor_threshold
  
  while (!vif_threshold_met) {
    
    # Calculate Spearman correlations
    remove_cor_variables <- function(data, priority_list, current_threshold = 0.7, verbose = TRUE) {
      
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
    uncor_vars <- remove_cor_variables(pred_values, 
                                       priority_list = rf_prelim_imp %>%
                                            filter(VME_Group == vme_group) %>%
                                            pull(Variable),
                                       current_threshold = cor_threshold_new,
                                       verbose = FALSE)
    
    # Calculate VIF values
    corvif <- function(data, verbose = TRUE) {
      data <- as.data.frame(data)
      
      form    <- formula(paste("fooy ~ ", paste(strsplit(names(data), " "), collapse = " + ")))
      data   <- data.frame(fooy = 1 + rnorm(nrow(data)), data)
      lm_mod  <- lm(form, data)
      
      if (verbose) print(data.frame(vif=car::vif(lm_mod)))
      return(data.frame(vif=car::vif(lm_mod)))
    }
    vif_values <- corvif(pred_values[, uncor_vars], verbose = FALSE)
    
    # Check if all VIF values are < 10
    if (all(vif_values$vif < 10)) {
      vif_threshold_met <- TRUE
      # uncor_vars[[vme_group]] <- selected_vars_new
      # vars_vif[[vme_group]] <- vif_values_new
      return(list(selected_vars = uncor_vars,
                  vif_values = vif_values,
                  final_cor_threshold = cor_threshold_new))
    }
    else {
      cor_threshold_new <- cor_threshold_new - 0.05
    }
  }
}

result <- select_vme_vars("bryozoan", cor_threshold = 0.7)

selected_vme_vars <- lapply(unique(comb_df_compl$VME_Group), select_vme_vars, cor_threshold = 0.7) %>%
  set_names(unique(comb_df_compl$VME_Group))


# Method 2: RFE ----