
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

sapply(rf_prelim, print)


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

# Model building with selected variables ----

# Create the 10 folds
set.seed(411)
folds <- caret::createFolds(vme_df$VME_P_A, k = 10, returnTrain = TRUE)

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
  # metrics <- c(
  #   "N" = cm[[3]][[1]],
  #   "Accuracy" = cm[[3]][[5]],
  #   "NIR" = cm[[3]][[6]],
  #   "P" = cm[[3]][[2]],
  #   "Sensitivity" = cm[[4]][[1]],
  #   "Specificity" = cm[[4]][[2]],
  #   "BalancedAcc" = cm[[4]][[11]],
  #   "TSS" = cm[[4]][[1]] + cm[[4]][[2]] - 1,
  #   "OptThreshold" = opttsh
  # )
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
rf_model <- caret::train(
  VME_P_A ~ .,
  data = vme_df_sel,
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
rf_metrics_results <- rf_model$results

# Extract variable importance
rf_varimp <- randomForest::importance(rf_model$finalModel) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Variable") %>%
  arrange(desc(MeanDecreaseGini))


# Predict on CMIP + bathy layers for whole area ----

vme_layers_rast <- terra::rast(vme_layers)

rf_pred <- terra::predict(vme_layers_rast, rf_model$finalModel, type = 'prob', na.rm = TRUE)


# Original model construction ----
## Extract predictions and process results

mdata <- vme_df %>%
  rename(resp = VME_P_A)

### Set number of cross-validation runs required
nruns <- 10 

### Set up training data ---

# Set up empty lists for looping through
train.sets <- list()
test.sets <- list()

# Split for 10 random subsets (list of row numbers), selects 90% of rows, keeping balance of classes equal, times = 10 runs
set.seed(123)
trainIndexP <- caret::createDataPartition(mdata$resp, p = 0.9, # Repeated sampling
                                          times = nruns)
trainIndexK <- caret::createFolds(mdata$resp, k = nruns) # K-fold

# Create 10 x separate train and tests sets using K-fold
for (j in 1:nruns){
  train.sets[[j]] <- mdata[unname(unlist(trainIndexK[-j])),]
  test.sets[[j]] <- mdata[trainIndexK[[j]],]
}

# Save the datasets
# save(train.sets,test.sets,file=paste0("data/SDM2024/Train_Test_",rvar,".RData"))


### Drop unnecessary layers from predictors ---
# dr <- names(predictors)
# dr <- dr[!dr %in% predsel]
# predictors <- dropLayer(predictors, dr)
# predictors


### Set up lists and tables for outputs ---

ffs <- list() # Create empty list for forests
imps <- list() # create empty list for importances
res <- list() # create empty list for results
tshs = NULL # create empty object for a list of optimal thresholds
cvpred <- NULL # create empty object for a stack of model class predictions
cvpred.cps <- list() # create empty list for model probability predictions
plotdata <- NULL # create empty object for partial plot data

# Create empty table for collecting all model performance statistics
class.res.all <- data.frame(Name=character(0),
                            Run=character(0),
                            N=character(0),
                            Acc=numeric(0),
                            NIR=numeric(0),
                            P=numeric(0),
                            Kappa=numeric(0),
                            Sensitivity=numeric(0),
                            Specificity=numeric(0),
                            BalancedAcc=numeric(0),
                            TSS=numeric(0),
                            stringsAsFactors =F)


# Below code is a loop that runs x 10
# Before running the whole loop, test the code by running just 1 model (run j=1)
j=1
for (j in 1:10){
  
  train <- train.sets[[j]]
  test <- test.sets[[j]]
  
  set.seed(123)
  ffs[[j]] <- randomForest::randomForest(resp ~ ., data = train,
                           ntree = 500, 
                           strata = resp,  # this makes sure sampling is balanced
                           replace = FALSE,
                           importance = TRUE, 
                           keep.forest = TRUE)
  
  results <- as.data.frame(rownames(test))  # check results
  results$actual <- test[[1]]  # adds column to results - P/A as factor
  results$PA <- as.character(test[[1]])  
  results$PA <- ifelse(results$PA == "Presence", 1, 0)  # change factor to numeric
  
  # Predict class with model j
  results$predicted <- as.data.frame(predict(ffs[j],test))[,1]  # outputs factor
  
  # Predicted probability is of PRESENCE
  results$predprob <- as.data.frame(predict(ffs[j],test,type='prob'))[,2]  # Check second column is presence!
  names(results)[1] <- "id"
  
  # Choose own optimal probability threshold: ID,observed, predicted. Various threshold methods, 
  # but 'Sens=Spec' returns equal amounts true and false positive classifications
  
  # require(PresenceAbsence)
  
  opttsh <- results %>%
    select(id, PA, predprob) %>%
    mutate(predprob2 = runif(nrow(.), 0, 1)) %>%
    PresenceAbsence::optimal.thresholds(opt.methods = 'Sens=Spec') %>%
    pull(predprob)
  tshs <- c(tshs, opttsh)
  
  # Presence by threshold, adds column for optimal thresholded class
  results <- results %>%
    mutate(optimal = as.factor(case_when(predprob >= opttsh ~ '1',
                                         TRUE ~ '0')))
  
  # Calculate confusion matrix for predictions by model i
  results.matrix <- caret::confusionMatrix(results$optimal, as.factor(results$PA), positive = '1')
  results.matrix
  
  
  # Get overall accuracy measures for model validation run i
  class.res.all[j,2] <- j
  class.res.all[j,3] <- nrow(test)
  class.res.all[j,4] <- results.matrix[[3]][[1]]
  class.res.all[j,5] <- results.matrix[[3]][[5]]
  class.res.all[j,6] <- results.matrix[[3]][[6]]
  class.res.all[j,7] <- results.matrix[[3]][[2]]
  class.res.all[j,8] <- results.matrix[[4]][[1]]
  class.res.all[j,9] <- results.matrix[[4]][[2]]
  class.res.all[j,10] <- results.matrix[[4]][[11]]
  class.res.all[j,11] <- results.matrix[[4]][[1]] + results.matrix[[4]][[2]] - 1 
  
  class.res.all$Name <- vme_group
  class.res.all
  
  imps[[j]] <- list(round(randomForest::importance(ffs[[j]]), 2))
  # require(pdp)
  
  p=1
  pcl = 1 # Presence class
  # for (p in 1:length(vme_vars)) {
  #   
  #   predname <- vme_vars[p]    
  #   pdata <- pdp::partial(ffs[[j]], pred.var = predname, which.class = pcl,
  #                    plot = FALSE, train = mdata, grid.resolution = 100, prob = TRUE)
  #   temp <- data.frame(Name = vme_group,
  #                      run = j,
  #                      predvar = predname,
  #                      class = pcl,
  #                      x = pdata[[1]],
  #                      y = pdata[[2]])
  #   plotdata <- rbind(plotdata,temp)
  #   
  # }
  
  vme_layers <- c(bathy_layers[vme_terrain_vars],
                            compact(cmip_layers[selected_vme_vars$selected_vars]))
  
  # Predict rasters
  rnn <-  paste0('Run',j) # Set layer name
  # Check if there is already a raster stack - if not create one
  if (is.null(cvpred)){
    # Probabilities for each class
    cvpred.cps[[rnn]] <- terra::predict(vme_vars,ffs[[j]],type='prob',index=1:numclass)
    # Presence/Absence raster from applying to the threshold to presence probability
    cvpred  <- stack(cut(cvpred.cps[[rnn]]$layer.2,breaks=c(-1,tshs[j],1)))
    cvpred <- cvpred - 1
    names(cvpred) <- rnn
  } else {
    # Probabilities for each class
    cvpred.cps[[rnn]] <- predict(predictors,ffs[[j]],type='prob',index=1:numclass)
    # Presence/Absence raster from applying to the threshold to presence probability
    tmpl <- cut(cvpred.cps[[rnn]]$layer.2,breaks=c(-1,tshs[j],1))-1
    names(tmpl) <- rnn
    cvpred <- addLayer(cvpred,tmpl)
  }
  
}


# Save models and validation results, plot data and importances
save(ffs,plotdata,class.res.all,imps,file = paste0('data/processed/SDM2024_rerun/',rvar,'/RF_Results_',rvar,'.RData'))


#### Look at the Validation statistics ----
require(matrixStats)
# Calculate averages and standard deviations for validation statistics
callavevalsB <- colMeans(class.res.all[,4:11])
callsdvalsB <- colSds(as.matrix(class.res.all[,4:11]))

# Combine values in a table
BCallvalsT <- data.frame(Accmean=round(callavevalsB[1],2),
                         Accsd=round(callsdvalsB[1],2),
                         Pmean=round(callavevalsB[3],2),
                         Psd=round(callsdvalsB[3],2),
                         Kmean=round(callavevalsB[4],2),
                         Ksd=round(callsdvalsB[4],2),
                         Sensmean=round(callavevalsB[5],2),
                         Senssd=round(callsdvalsB[5],2),
                         Specmean=round(callavevalsB[6],2),
                         Specsd=round(callsdvalsB[6],2),
                         BAmean=round(callavevalsB[7],2),
                         BAsd=round(callsdvalsB[7],2),
                         TSSmean=round(callavevalsB[8],2),
                         TSSsd=round(callsdvalsB[8],2))

# Rename columns
names(BCallvalsT) <- c("Accmean","Accsd","Pmean","Psd","Kmean","Ksd",
                       "Sensmean","Senssd","Specmean","Specsd",
                       "BAmean","BAsd","TSSmean","TSSsd")


# Print table
BCallvalsT

asg.perf <- data.table(N = nrow(train),
                       'Sensitivity'= paste(BCallvalsT$Sensmean, '/u00B1',BCallvalsT$Senssd),
                       'Specificity' =  paste(BCallvalsT$Specmean, '/u00B1',BCallvalsT$Specsd),
                       'Kappa' = paste(BCallvalsT$Kmean, '/u00B1',BCallvalsT$Ksd) ,
                       'Balanced Accuracy'= paste(BCallvalsT$BAmean, '/u00B1',BCallvalsT$BAsd),
                       'TSS'=paste(BCallvalsT$TSSmean, '/u00B1',BCallvalsT$TSSsd))

asg.perf[, data.table(t(.SD), keep.rownames=TRUE),] %>%
  kbl('html',digits = 2,escape = FALSE, col.names = c('Statistic','Mean /u00B1 SD'),
      caption='Performance statistics') %>%
  kable_classic(full_width = F, position = "left",fixed_thead = T) %>%
  row_spec(0, bold = T)  %>%
  column_spec(1:2, width = "3cm") 


