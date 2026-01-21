
# RF model building following variable elimination

# BNAM SDM model building ----

# We are building the model inside the 10 loops and predict with it and get the validation at the 
# same time.

### Set up data and constants ----

# Set name of response variable to be used in results tables
tax = 'BlackCorals' #'SeaPens'
# Set the name of the positive class
pcl = '1'
# Predictor variable constants
preds <-   predsel # for backwards compatibility  
facvars <- NULL # gear code when running full model
predselnf <- predsel # this code doesn't do anything here, is backwards compatibility

# Rename the response column (just to fit with existing code)
setnames(mdata,1,'resp')

### Set number of cross-validation runs required ---
nruns <- 10 

### Set up training data ---

# Set up empty lists for looping through
train.sets <- list()
test.sets <- list()

# Split for 10 random subsets (list of row numbers), selects 90% of rows, keeping balance of classes equal, times = 10 runs
trainIndexP <- createDataPartition(mdata$resp, p = .90, # Repeated sampling
                                   times = nruns)
trainIndexK <- createFolds(mdata$resp,k=nruns) # K-fold

# Create 10 x separate train and tests sets using K-fold
for (j in 1:nruns){
  
  train.sets[[j]] <- mdata[unname(unlist(trainIndexK[-j])),]
  test.sets[[j]] <- mdata[trainIndexK[[j]],]
  
  next}

# Save the datasets
save(train.sets,test.sets,file=paste0("data/SDM2024/Train_Test_",rvar,".RData"))


### Drop unnecessary layers from predictors ---
dr <- names(predictors)
dr <- dr[!dr %in% predsel]
predictors <- dropLayer(predictors, dr)
predictors


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

for (j in 1:10){
  
  train <- train.sets[[j]]
  test <- test.sets[[j]]
  
  ffs[[j]] <- randomForest(resp ~.,data=train,
                           ntree=500, 
                           strata=resp,
                           replace=FALSE,
                           importance=T, 
                           keep.forest= T)
  
  
  results <- as.data.frame(rownames(test)) #check results
  results$actual <- test[[1]] #adds column to results - P/A as factor
  results$PA <- as.numeric(as.character(test[[1]])) # changes factor to numeric
  
  # Predict class with model j
  results$predicted <- as.data.frame(predict(ffs[j],test))[,1] # outputs factor
  
  # Predicted probability is of PRESENCE
  results$predprob <- as.data.frame(predict(ffs[j],test,type='prob'))[,2] # Check second column is presence!
  names(results)[1] <- "id"
  
  # Choose own optimal probability threshold: ID,observed, predicted. Various threshold methods, 
  # but 'Sens=Spec' returns equal amounts true and false positive classifications
  
  require(PresenceAbsence)
  
  opttsh <- results %>%
    dplyr::select(id,PA,predprob) %>%
    optimal.thresholds(opt.methods = 'Sens=Spec') %>%
    pull(predprob)
  tshs <- c(tshs, opttsh)
  
  # Presence by threshold, adds column for optimal thresholded class
  results <- results %>%
    mutate(optimal=as.factor(case_when(predprob>=opttsh ~ '1',
                                       TRUE~'0')))
  
  # Calculate confusion matrix for predictions by model i
  results.matrix <- confusionMatrix(results$optimal, results$actual,positive = '1',)
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
  
  class.res.all$Name <- tax
  class.res.all
  
  
  imps[[j]] <- list(round(randomForest::importance(ffs[[j]]), 2))
  
  require(pdp)
  
  for (p in 1:length(predselnf)) {
    
    pdata <- partial(ffs[[j]],pred.var = predselnf[p],which.class = pcl,
                     plot = FALSE,train=mdata,grid.resolution=100,prob = TRUE)
    predname <- predselnf[p]
    temp <- data.frame(Name=tax,run=j,predvar=predselnf[p],class=pcl,x=pdata[[1]],y=pdata[[2]])
    plotdata <- rbind(plotdata,temp)
    
    
  }
  
  ## Predict rasters
  rnn <-  paste0('Run',j) # Set layer name
  # Check if there is already a raster stack - if not create one
  if (is.null(cvpred)){
    # Probabilities for each class
    cvpred.cps[[rnn]] <- predict(predictors,ffs[[j]],type='prob',index=1:numclass) 
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
  
  next
}



# ESIAOI RF model building ----

# Cross-validation
set.seed(411)
cv_folds <- createFolds(div, k = 10, list = TRUE)

# Function to evaluate regression model performance
evaluate_rf_regression <- function(data, response_var, cv_folds) {
  # Initialize performance metrics
  rmse_values <- numeric(length(cv_folds))
  r_squared_values <- numeric(length(cv_folds))
  variable_importance <- list()
  
  for (i in seq_along(cv_folds)) {
    # Split data
    train_data <- data[-cv_folds[[i]], ]
    test_data <- data[cv_folds[[i]], ]
    
    # Train model
    rf_model <- randomForest(
      formula = as.formula(paste(response_var, "~ .")),
      data = train_data,
      ntree = 500,
      importance = TRUE
    )
    
    # Store variable importance
    variable_importance[[i]] <- importance(rf_model)
    
    # Make predictions
    predictions <- predict(rf_model, test_data)
    
    # Calculate metrics
    rmse_values[i] <- sqrt(mean((predictions - test_data[[response_var]])^2))
    sst <- sum((test_data[[response_var]] - mean(test_data[[response_var]]))^2)
    sse <- sum((predictions - test_data[[response_var]])^2)
    r_squared_values[i] <- 1 - (sse / sst)
  }
  
  # Average variable importance across folds
  avg_importance <- Reduce(`+`, variable_importance) / length(variable_importance)
  
  return(list(
    # rf_model = rf_model,
    list(
      mean_rmse = mean(rmse_values),
      sd_rmse = sd(rmse_values),
      se_rmse = sd(rmse_values) / sqrt(length(cv_folds)),
      mean_r_squared = mean(r_squared_values),
      sd_r_squared = sd(r_squared_values),
      se_r_squared = sd(r_squared_values) / sqrt(length(cv_folds))
    ),
    importance = avg_importance
  ))
}
