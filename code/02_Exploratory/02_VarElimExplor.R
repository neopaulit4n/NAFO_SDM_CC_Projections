
# Testing variable elimination methods

# Load data ----

library(tidyverse)

# Load response ----
resp_df <- read_csv("data/cleaned/VME_group_PA_df.csv", show_col_types = FALSE)

# Load predictors ----

## Load terrain variables (static) ----
bathy_layers <- list.files(path = "data/raw/BNAM_Data_From_Cam/Bathymetry_Terrain_From_NAFO_SharePoint", 
                           pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(., nm = basename(.) %>% tools::file_path_sans_ext()) %>%
  lapply(terra::rast) %>%
  lapply(terra::project, "EPSG:4326")

names(bathy_layers) <- gsub("GEBCO2024_FS005_StudyArea_","",names(bathy_layers))
names(bathy_layers)[1] <- "FS005"

## Load BNAM layers (will use these to form predictions, decide which variables to select) ----
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
pred_df <- lapply(c(bathy_layers, bnam_layers), function(layer) {
  terra::extract(layer, 
                 select(resp_df, Start_Long_DD, Start_Lat_DD)) %>%
    select(-ID)
}) %>%
  bind_cols()
colnames(pred_df) <- c(names(bathy_layers), names(bnam_layers))


# Combine predictor and response dataframes ----
comb_df <- bind_cols(resp_df, pred_df) %>%
  mutate(VME_P_A = factor(VME_P_A, levels = c(0, 1), labels = c("Absence", "Presence")))

# Remove NA
comb_df_compl <- comb_df %>%
  drop_na()
rm(comb_df)

# comb_df_miss <- comb_df[which(!complete.cases(comb_df)),]


# Method 1: Replicating previous SDM variable elimination approach ----
# 1: Run preliminiary RF models with all variables
# 2: Determine variable importance
# 3: Remove correlated variables based on importance ranking
# 4: VIF remaining variables

library(randomForest)

## Preliminary RF model test ----
rf_prelim_form <- as.formula(paste("VME_P_A ~", 
                                 paste(colnames(pred_df), collapse = " + ")))

## Preliminary RF model run (save/load) ----
# rf_prelim <- lapply(unique(comb_df_compl$VME_Group), function(vme_group) {
#   set.seed(804)
#   rf_model <- randomForest::randomForest(formula = rf_prelim_form,
#                                          data = filter(comb_df_compl, VME_Group == vme_group),
#                                          importance=TRUE)
#   return(rf_model)
# }) %>%
#   set_names(unique(comb_df_compl$VME_Group))

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

# Show var imp plot for each VME group, order each by descending importance
ggplot(rf_prelim_imp, aes(x = reorder(Variable, MeanDecreaseGini), 
                           y = MeanDecreaseGini, 
                           fill = VME_Group)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ VME_Group, scales = "free_y") +
  coord_flip() +
  theme_bw() +
  labs(x = "Predictor Variable", y = "Mean Decrease in Gini Index") +
  theme(legend.position = "none")

randomForest::varImpPlot(rf_prelim$black_corals)

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

select_vme_vars_m1 <- function(vme_group, cor_threshold = 0.7) {
  
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
    
    uncor_vars <- remove_cor_variables(pred_df, 
                                       priority_list = rf_prelim_imp %>%
                                            filter(VME_Group == vme_group) %>%
                                            pull(Variable),
                                       current_threshold = cor_threshold_new,
                                       verbose = FALSE)
    
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
    vif_values <- corvif(pred_df[, uncor_vars], verbose = FALSE)
    
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

selected_vme_vars_m1 <- lapply(unique(comb_df_compl$VME_Group), select_vme_vars_m1, cor_threshold = 0.7) %>%
  set_names(unique(comb_df_compl$VME_Group))

saveRDS(selected_vme_vars_m1, file = "data/processed/02_Exploratory/selected_vme_vars_m1.rds")
selected_vme_vars_m1 <- readRDS("data/processed/02_Exploratory/selected_vme_vars_m1.rds")

# Method 2: RFE ----

# Remove correlated variables > 0.7
remove_cor_m2 <- function(x, cor_threshold = 0.7, verbose = TRUE) {
  
  cor_threshold_new <- cor_threshold
    
  # Calculate Spearman correlations
  remove_cor_variables <- function(data, current_threshold = 0.7, verbose = FALSE) {
    
    vars <- colnames(data)
    
    # Calculate correlation
    cor_matrix <- matrix(0, nrow = ncol(data), ncol = ncol(data))
    rownames(cor_matrix) <- colnames(cor_matrix) <- colnames(data)
    
    numeric_cor <- cor(data[, vars, drop = FALSE], 
                       method = "spearman",
                       use = "pairwise.complete.obs")
    cor_matrix[vars, vars] <- numeric_cor
    
    # Prioritise variables
    # priority_dict <- setNames(seq_along(priority_list), priority_list)
    all_vars <- colnames(data)
    
    # Assign priorities (variables in priority list come first, by their order)
    # var_priority <- sapply(all_vars, function(v) {
    #   if (v %in% priority_list) {
    #     return(priority_dict[v])
    #   } else {
    #     return(length(priority_list) + 1) # Lower priority for non-listed vars
    #   }
    # })
    
    # Sort variables by priority
    sorted_vars <- all_vars#[order(var_priority)]
    
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
  
  uncor_vars <- remove_cor_variables(x,
                                     # priority_list = rf_prelim_imp %>%
                                     #   filter(VME_Group == vme_group) %>%
                                     #   pull(Variable),
                                     current_threshold = cor_threshold_new,
                                     verbose = verbose)
}

uncor_vars_m2 <- remove_cor_m2(distinct(pred_df), cor_threshold = 0.7, verbose = TRUE)

# selected_vme_vars_m2 <- lapply(unique(comb_df_compl$VME_Group), select_vme_vars_m2, cor_threshold = 0.7) %>%
#   set_names(unique(comb_df_compl$VME_Group))




# Set control parameters
# rfe_ctrl <- caret::rfeControl(functions = caret::rfFuncs,  # "rfFuncs" are built-in to caret
#                               rerank = TRUE,
#                               method = "repeatedcv",  # cross-validation
#                               number = 2,  # 10-fold
#                               repeats = 2,
#                               saveDetails = FALSE,
#                               allowParallel = FALSE,
#                               verbose = TRUE)

# Recursive feature selection
library(caret)
library(randomForest)

rfe_x <- as.data.frame(comb_df_compl[comb_df_compl$VME_Group == "sponges", uncor_vars_m2])
rfe_y <- as.factor(filter(comb_df_compl, VME_Group == "sponges")$VME_P_A)
levels(rfe_y) <- c("absent", "present")
# rfe_y <- as.numeric(rfe_y)-1
# rfe_res <- caret::rfe(x = rfe_x,
#                       y = rfe_y,
#                       sizes = c(1:length(selected_vme_vars_m2)),  # test all subset sizes
#                       rfeControl = rfe_ctrl,
#                       metric = "Accuracy")

# Create custom RF functions with explicit parameters
customRF <- list(
  summary = defaultSummary,
  fit = function(x, y, first, last, ...) {
    randomForest(x, y, mtry = floor(sqrt(ncol(x))), ntree = 500)
  },
  pred = function(object, x) {
    predict(object, x)
  },
  rank = function(object, x, y) {
    vimp <- importance(object)
    vimp <- vimp[order(vimp[, 1], decreasing = TRUE), , drop = FALSE]
    data.frame(var = rownames(vimp), Overall = vimp[, 1])
  },
  selectSize = pickSizeBest,
  selectVar = pickVars
)

rfe_ctrl <- rfeControl(
  functions = customRF,
  method = "cv",
  number = 5,
  allowParallel = FALSE,
  verbose = TRUE
)

set.seed(123)
rfe_res <- rfe(x = rfe_x,
               y = rfe_y,
               sizes = c(1:16),
               rfeControl = rfe_ctrl)
print(rfe_res)

ggplot(rfe_res, metric = "Accuracy")
ggplot(rfe_res, metric = "Kappa")

# Compare random forest with method 1 vs method 2 ----

# Comparison using cross-validation using identical folds for both methods

library(pROC)

# Custom function with ALL metrics
rf_metrics <- function(data, lev = NULL, model = NULL) {
  cm <- caret::confusionMatrix(data$pred, data$obs, positive = lev[2])
  
  sens <- as.numeric(cm$byClass["Sensitivity"])
  spec <- as.numeric(cm$byClass["Specificity"])
  tss <- sens + spec - 1
  
  roc_obj <- pROC::roc(data$obs, data[, lev[2]], quiet = TRUE)
  
  # Return as named vector with as.numeric to strip attributes
  out <- c(as.numeric(cm$overall["Accuracy"]),
           as.numeric(cm$overall["Kappa"]),
           sens,
           spec,
           tss,
           as.numeric(auc(roc_obj)))
  
  names(out) <- c("Acc", "Kap", "Sens", "Spec", "TSS", "AUC")
  
  return(out)
}


# Set up training control
set.seed(123)
train_index <- caret::createFolds(y = filter(comb_df_compl, VME_Group == "sponges")$VME_P_A, 
                                          k = 10, 
                                          returnTrain = TRUE,
                                  list = TRUE)

set.seed(123)
train_control <- trainControl(index = train_index,
                              method = "cv",
                     number = 10,
                     classProbs = TRUE,
                     summaryFunction = rf_metrics,
                     savePredictions = "final")

# Train both models with identical CV folds
rf_m1 <- caret::train(form = as.formula(paste("VME_P_A ~", 
                                              paste(selected_vme_vars_m1$sponges$selected_vars, collapse = " + "))),
                      data = filter(comb_df_compl, VME_Group == "sponges"),
                      method = "rf",
                      trControl = train_control,
                      metric = "AUC")

rf_m2 <- caret::train(form = as.formula(paste("VME_P_A ~", 
                                              paste(rfe_res$optVariables, collapse = " + "))),
                      data = filter(comb_df_compl, VME_Group == "sponges"),
                      method = "rf",
                      trControl = train_control,
                      metric = "AUC")

# Compare all metrics
rf_mcomp_res <- resamples(list(Method1 = rf_m1, Method2 = rf_m2))
saveRDS(rf_mcomp_res, file = "data/processed/02_Exploratory/VarElimMethodComparisonResults_sponges.rds")

# View summary statistics
summary(rf_mcomp_res)

# Statistical comparison of differences
diff_results <- diff(rf_mcomp_res)
summary(diff_results)

# Visualise comparisons
dotplot(rf_mcomp_res)
bwplot(rf_mcomp_res)




# Method 3: SMOTE to fix data imbalance between presence/absence ----

comb_df_compl_bc <- filter(comb_df_compl, 
                           VME_Group == "black_corals")

# MDS to detect outliers among presence points
set.seed(123)
mds_bc <- cmdscale(d = dist(comb_df_compl_bc[, colnames(pred_df)]),
                   k = 2)

ggplot(as.data.frame(mds_bc), aes(x = V1, y = V2)) +
  geom_point() +
  theme_bw() +
  labs(x = "MDS1", y = "MDS2") +
  ggtitle("MDS of Black Coral Presence Points")

# No obvious outliers

# MDS of all data and colour-code by presence/absence
comb_df_compl_bc_all <- filter(comb_df_compl, 
                               VME_Group == "black_corals") %>%
  select(all_of(colnames(pred_df)), VME_P_A) %>%
  mutate(source = "Original")
set.seed(123)
mds_bc_all <- cmdscale(d = dist(comb_df_compl_bc_all[, colnames(pred_df)]),
                       k = 2)
mds_bc_all_df <- as.data.frame(mds_bc_all) %>%
  mutate(VME_P_A = comb_df_compl_bc_all$VME_P_A)

ggplot(as.data.frame(mds_bc_all_df), aes(x = V1, y = V2, colour = VME_P_A)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  labs(x = "MDS1", y = "MDS2") +
  ggtitle("MDS of Black Coral Presence/Absence Points") +
  scale_colour_manual(values = c("0" = "blue", "1" = "red"),
                      labels = c("0" = "Absence", "1" = "Presence"),
                      name = "VME Presence/Absence")


smote_df <- ROSE::ROSE(formula = rf_prelim_form, 
                       data = comb_df_compl_bc_all, 
                       p = 0.5, 
                       seed = 1)$data %>%
  mutate(source = "SMOTE")

mds_bc_smote <- cmdscale(d = dist(smote_df[, colnames(pred_df)]),
                       k = 2) %>%
  as.data.frame()
mds_bc_smote_df <- mutate(mds_bc_smote, VME_P_A = smote_df$VME_P_A) %>%
  filter(VME_P_A == "1") %>%
  mutate(source = "SMOTE")

smote_orig_df <- bind_rows(comb_df_compl_bc_all, smote_df)

mds_bc_smote_all <- cmdscale(d = dist(smote_orig_df[, colnames(pred_df)]),
                           k = 2)


ggplot(mds_bc_smote_df, aes(x = V1, y = V2, colour = VME_P_A)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  labs(x = "MDS1", y = "MDS2") +
  ggtitle("MDS of Black Coral SMOTE Presence/Absence Points") +
  scale_colour_manual(values = c("0" = "blue", "1" = "red"),
                      labels = c("0" = "Absence", "1" = "Presence"),
                      name = "VME Presence/Absence")

ggplot(mds_bc_smote_all, aes(x = V1, y = V2, colour = smote_orig_df$source)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  labs(x = "MDS1", y = "MDS2") +
  ggtitle("MDS of Black Coral Original vs SMOTE Data") +
  scale_colour_manual(values = c("Original" = "blue", "SMOTE" = "red"),
                      name = "Data Source")


rf_m3 <- randomForest::randomForest(formula = as.formula(paste("VME_P_A ~", 
                                                                 paste(colnames(pred_df), collapse = " + "))),
                                    data = smote_df,
                                    importance=TRUE)
rf_m3

# Try with repeated cross-validation
set.seed(123)

train_control <- caret::trainControl(method = "repeatedcv", 
                              number = 5, 
                              repeats = 3)

train_index <- caret::createDataPartition(y = smote_df$VME_P_A, 
                                          p = 0.8, 
                                          list = FALSE)
training_df <- smote_df[train_index, ]
testing_df <- smote_df[-train_index, ]

rf_m3_cv <- caret::train(form = as.formula(paste("VME_P_A ~", 
                                            paste(colnames(pred_df), collapse = " + "))),
                  data = smote_df,
                  method = "rf",
                  trControl = train_control,
                  importance = TRUE,
                  metric = "Accuracy")
rf_m3_cv

var_imp <- caret::varImp(rf_m3_cv)

# Plot var imp
ggplot(var_imp, aes(x = reorder(rownames(var_imp$importance), Overall), 
                           y = Overall)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_bw() +
  labs(x = "Predictor Variable", y = "Importance") +
  ggtitle("Variable Importance from RF with SMOTE Data")


