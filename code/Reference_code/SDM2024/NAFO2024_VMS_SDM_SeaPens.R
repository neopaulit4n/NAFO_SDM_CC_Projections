###############################################################################
######          NAFO Sea pen distribution models 2024              ############
###############################################################################

#Load packages and library files
library(raster)
# library(maptools)            # no longer available
library(randomForest)
library(ranger)
library(tidyverse)
library(vtable)
library(pdp)
library(sf)
library(data.table)
library(ggcorrplot)
library(patchwork)
library(caret)
# library(dsmextra)            # no longer available on CRAN


## Colour palette
cpl <- c('#d4ebe7','#cbbcbb','#f5f1f1','#172957','#66afad')
names(cpl) <- c('lt','dbe','lbe','dbl','dt')

OneB_theme <-
  ggplot2::theme(axis.title.y = element_text(vjust=4,  size=12,colour="black"),
                 axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
                 axis.text.x  = element_text(vjust=0.5, size=12,colour="black"),
                 axis.title.x  = element_text(vjust=-4, size=12,colour="black"),
                 strip.background = element_rect(fill=cpl['lbe']),
                 strip.text.x = element_text(size=12, face="bold"),
                 panel.grid.major = element_line(colour=cpl['lbe']),
                 panel.grid.minor = element_line(colour=cpl['lbe']),
                 panel.background = element_rect(fill="white"),
                 plot.margin = ggplot2::margin(0.5, 0.5, 1, 0.5, "cm"))



######################################################################
###### STEP 1: IMPORT THE RESPONSE AND ENVIRONMENTAL DATA ############
######################################################################

# Set working directory
# wdir = ("C:/Users/AD06/OneDrive - CEFAS/VME/NAFO2024")
# setwd(wdir)

# Set response variable title for model outputs
rvar = "BlackCorals"  #'SeaPens'


### Environmental data -----

# Directory containing environmental rasters
# rasterdir = "/ENVDATA/FINAL"  
rasterdir <- "data/raw/SDM2024/FINAL Environmental_Variables"

# List of raster files
# predictorfiles = list.files(path = paste(wdir, rasterdir, sep=""), pattern = "\\.tif$", full.names=T)
predictorfiles = list.files(path = rasterdir, 
                            pattern = "tif", 
                            full.names = TRUE,
                            recursive = TRUE)
# Now read the raster data (create a raster stack)
predictors = raster::stack(predictorfiles,RAT=F)
# Confirm raster stack with all raster layers present
predictors
names(predictors)
names(predictors)[72:73] <- c("NRA_fishing_effort_1km","NRA_fishing_effort_5km")
# Plot raster files
plot(predictors)

# coordinate system for rasters
rprj = st_crs(predictors) # or if is not included in data set manually e.g. CRS('+proj=longlat +datum=WGS84 +no_defs')

#######################

### Response data ---

# Directory containing response data 
respdir = ("data/raw/SDM2024/FINAL Response Variables")


# Read and investigate csv
responsefile = read.csv(paste(respdir,"Black corals/black_corals.csv", sep="/"), header=TRUE)
head(responsefile)
dim(responsefile)

# Response variable name
respvar = "VME_P_A"
# Coordinate variable names
xyvars = c("Start_Long_DD","Start_Lat_DD")

# Select response and coordinate columns
responsedata = data.frame(pa=responsefile[,respvar],
                          x=responsefile[,xyvars[1]], 
                          y=responsefile[,xyvars[2]])
head(responsedata)
dim(responsedata)

# Check response column is a factor
if (!is.factor(responsedata$pa)) {
  responsedata$pa <- as.factor(responsedata$pa)
}
# Check factor levels 
levels(responsedata$pa)

# Final data removing NAs
response = responsedata[complete.cases(responsedata),]

## Convert to spatial --
# Define coordinate system
pprj = CRS('+proj=longlat +datum=WGS84 +no_defs')
# Convert to sf
response_sp = sf::st_as_sf(response,coords=c('x','y'),crs=sf::st_crs(4326))
response_sp

# Check response and environmental variables are in the same coordinate system
if (pprj != rprj) {
  response_sp <- st_transform(response_sp, rprj)
}

st_write(response_sp,dsn = paste0('data/SDM2024/processed/',rvar,'_Data/',rvar,'_Data.shp'),append = FALSE)


#########################################################################
#### STEP 2. EXTRACT THE VALUES OF PREDICTORS AT RESPONSE LOCATIONS #####
#########################################################################

### Extract the values from each predictor for each location in the response data, and put results in a new dataframe
p.data = raster::extract(predictors, response_sp)
sdata = data.frame(response, p.data) #adds all the extracted values to the existing dataframe
head(sdata)
prnames = colnames(p.data)  #get names from new columns
str(sdata)

### Labels to use for environmental variables
# names of predictor columns
# prnames
# # Read csv file with two columns 'variable' with predictor column names and 'label' with labels to use for plotting
# varlabs <- read.csv("Models/varnames.csv")
# # Convert to named vector
# envlab <- varlabs$label
# names(envlab) <- varlabs$variable
# envlab


#########################################################################
#### STEP 3. QUALITY CONTROL ############################################
#########################################################################

# Ensure that all observations (locations) have data values for all variables
sdata = sdata[complete.cases(sdata),]  #'complete.cases' command returns only those rows in the dataframe that have non-NA 
head(sdata)
dim(sdata)
summary(sdata)
# Values for all columns
str(sdata)
nrow(sdata)
missingdata = data.frame(response, p.data)  #creates a new dataframe called 'missingdata'
missingdata = missingdata[!complete.cases(missingdata),]  #!interested in rows that are NOT complete cases
dim(missingdata)

#### THIS IS POTENTIAL TO INCLUDE IF WANT ####
# 
# ## Increase prevalence by sub-sampling absence data --
# 
# # calculate number of presences
# npres <- sdata %>%
#           filter(pa=='1') %>%
#           nrow()
# # Calculate prevalence
# preval <- npres/nrow(sdata)
# preval
#             
# # If prevalence threshold of 5% is not met subsample absences to match it 
# if (preval <0.05) {
# sdata <- sdata %>%
#           filter(pa=='1') %>% 
#           bind_rows(sdata %>%
#                       filter(pa=='0') %>%
#                       sample_n(20*npres))
# }

#################################################

# Save the data frame and labels
save(sdata,#envlab,
     file = paste0('data/SDM2024/processed/',rvar,'_Data.RData')) 

#########################################################################
#### STEP 4. VARIABLE ELIMINATION/SELECTION #############################
#########################################################################

# This is for backwards compatibility for code for now
numvars <- prnames

# Define the number of class levels
numclass <- nlevels(sdata[[1]])

### Data exploration ----

## Summary statistics --
# All
sumtable(sdata, simple.kable = TRUE)
# Environmental variables for presences only
sumtable(sdata[sdata$pa=="1",], simple.kable = TRUE)

## Covariance of environmental variables --
# Correlation
corr <- cor(sdata[-1])
# colnames(corr) <- envlab[colnames(corr)]
# rownames(corr) <- envlab[rownames(corr)]
# Correlation plot
corplot <- ggcorrplot::ggcorrplot(corr, method='circle',type = 'upper',hc.order = TRUE)
corplot


### Preliminary full model to compare variable importance ----

## Build model --
prelRF <- randomForest(pa~.,
                       data=sdata,
                       importance=TRUE)
prelRF

## Extract importance and place in order --
full.importance <- data.table(Predictor=rownames(prelRF$importance),prelRF$importance)
full.importance <- full.importance[order(full.importance[,MeanDecreaseGini],decreasing=T),]
full.importance

## Plot partial dependence --

# Set up object to save plot data to
plotdata <- NULL
predselnf <- numvars # for backwards compatibility for now

# Define class to plot
cl = '1'

# Loop through environmental variables to create data for partial plots
for (j in 1:length(predselnf)) {
    
    pdata <- partial(prelRF,pred.var = predselnf[j],which.class = cl,
                     plot = FALSE,train=sdata,grid.resolution=100,prob = TRUE)
    predname <- predselnf[j]
    temp <- data.frame(predvar=predselnf[j],class=cl,x=pdata[[1]],y=pdata[[2]])
    plotdata <- rbind(plotdata,temp)
  
}

# Round values
plotdata.r <- plotdata %>%
                mutate(y= round(y, 1))

# Create list of partial plots
fullRP.list <- list()

# Loop through each predictor variable to plot partial dependence and add to the list
for (i in predselnf) {
  
  fullRP.list[[i]] <- ggplot(plotdata[plotdata$predvar==i,],aes(x=x,y=y,col=class)) +
    geom_smooth(linewidth=0.8,se=FALSE,span = 0.3,col='#172957') +
    facet_wrap(~ predvar,scales = "free_x", ncol=3) +
    ylim(c(min(c(0,plotdata$y)),max(plotdata$y))) +
    theme(axis.title.y = element_text(vjust=0.5, size=12,colour="black"),
          axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
          axis.text.x  = element_text(vjust=0.5, size=12,colour="black"),
          axis.title.x  = element_blank(),
          plot.title =  element_text(size=12,colour="white", face = "bold",vjust=2),
          strip.background = element_rect(fill="grey90"),
          strip.text.x = element_text(size=12, face="bold"),
          panel.grid.major = element_line(colour="grey80"),
          panel.grid.minor = element_line(colour="grey80"),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key = element_rect(fill = NA),
          legend.text = element_text(size=12,colour="black"))
}

# Write a pdf with all partial plots to check
pdf(paste0("output/SDM2024/combined_plots_",rvar,".pdf"), width = 8, height = 11)
# Loop through the plots and arrange them in a 2x4 grid, 8 plots per page
for (i in seq(1, length(fullRP.list), by = 8)) {
  combined_plot <- wrap_plots(fullRP.list[i:min(i+7, length(fullRP.list))], ncol = 2, nrow = 4)
  print(combined_plot)
}
# Close the PDF device
dev.off()


### Select uncorrelated variables to keep ----

## Correlation matrix with variables in order of full model importance --
vl <- full.importance[[1]] # Variable list
vl <- vl[vl %in% numvars] # Numeric variables only - for back compatibility
cr <- cor(sdata[,vl]) # correlation matrix
# Remove variables correlated to a higher importance variable
for(j in 1:length(cr[1,])){
  if (j == 1){
    pl <- c(names(cr[j,][1]),names( cr[j,][sqrt((cr[j,])^2)<0.65]))
    pl1 <- pl
  } else if (names(cr[j,])[j] %in% pl1){
    rem <- names(cr[j,-c(1:j)][sqrt((cr[j,-c(1:j)])^2)>0.65])
    if (length(rem) != 0L){  
      pl <- pl[!pl %in% rem]
    }
  }
  next
}
# Show list of selected variables
pl

## Calculate Variance Inflation Factors (vif) for selected variables --
# Null model vif function
corvif =  function(dataz) {
  dataz <- as.data.frame(dataz)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1 + rnorm(nrow(dataz)) ,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(data.frame(vif=car::vif(lm_mod)))
}

# Table kept variables and their vif
crval <- as.data.frame(pl)
crval$vif <- corvif(sdata[,pl])
crval

## This process can be redone to decrease allowed correlation if vif values 
## remain too high

### Choose the set of variables to use in model ----

# List selected variables
predsel <- crval[[1]] # keeping this line for backwards compatibility

# List of all variables including response
clms <- c(names(sdata)[1],predsel)

# Data to use in model
mdata <- sdata[,clms]
summary(mdata)


### Save preliminary model and model data ----
save(prelRF,full.importance,crval,clms,mdata,plotdata.r, file=paste0("data/SDM2024/processed/RF_Prelim_",rvar,".RData"))
save(sdata,mdata,clms,predsel,file = paste0("data/SDM2024/processed/",rvar,'_Data.RData'))



#########################################################################
#### STEP 5. BUILD MODEL & VALIDATION ###################################
#########################################################################

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


#### Plot variable importance ----

# Importance plot
imppl <- data.table(Var=rownames(imps[[1]][[1]]))

for (i in 1:10){
  
  imppl <- cbind(imppl,as.data.table(imps[[i]][[1]])[,MeanDecreaseGini])
  
}

setnames(imppl,c('Var','Imp1','Imp2','Imp3','Imp4','Imp5','Imp6','Imp7','Imp8','Imp9','Imp10'))


imppl[,
      c("Mean",'Sd','Se') := 
        .(rowMeans(.SD, na.rm = TRUE), 
          apply(.SD, 1, sd, na.rm = TRUE),
          apply(.SD, 1, plotrix::std.error, na.rm = TRUE)), 
      .SDcols = 2:11]

imppl[,Var:=factor(Var,levels=Var[order(Mean)])]


impplot <-  ggplot(imppl,(aes(x=Var,y=Mean))) +
  geom_bar(stat = 'identity',fill='#66afad', col='#172957',) +
  # scale_x_discrete(labels=envlab[levels(imppl$Var)]) +
  geom_linerange(inherit.aes=FALSE,
                 aes(x=Var, ymin=Mean-Se, ymax=Mean+Se), 
                 colour='#172957', alpha=0.9, linewidth=1.3) +
  ylab(label = 'Mean decrease in Gini coefficient') +
  coord_flip() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_text(vjust=0.5,hjust = 1, size=12,colour="black"),
        axis.text.x  = element_text(vjust=0.5, size=12,colour="black"),
        axis.title.x  = element_text(vjust=-4, size=12,colour="black"),
        plot.title =  element_text(size=12,colour="white", face = "bold",vjust=2),
        strip.background = element_rect(fill="grey90"),
        strip.text.x = element_text(size=12, face="bold"),
        panel.grid.major = element_line(colour="grey80"),
        panel.grid.minor = element_line(colour="grey80"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size=12,colour="black"),
        plot.margin = ggplot2::margin(0.5, 0.5, 1, 0.5, "cm"),)
impplot

ggsave(impplot,
       filename=paste0('output/SDM2024/',rvar,'_VariableImportance.png'),
       device = 'png', width = 15, height=16, units='cm', dpi=300, scale=1)


#### Partial response plots ----

pplotdata <- plotdata # this is here if need to make any changes to plotdata

# Create list for plots
cvRP.list <- list()

# Loop through predictors to create plots
for (i in predsel) {
  
  mxy <- max(pplotdata[pplotdata$class==pcl,'y'])
  
  cvRP.list[[i]] <- ggplot(pplotdata[pplotdata$predvar==i & pplotdata$class==pcl,],aes(x=x,y=y,group=run)) +
    geom_smooth(method='loess',linewidth=0.01,se=FALSE,span = 0.2,col='#66afad') +
    geom_smooth(inherit.aes=FALSE,aes(x=x,y=y),
                method='loess',linewidth=0.9,se=FALSE,span = 0.2,col='#172957') +
    facet_wrap(~ predvar,scales = "free_x",ncol =3) + #labeller = labeller(predvar = envlab)) +
    ylim(c(min(c(0,pplotdata$y)),mxy)) +
    OneB_theme +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
          axis.text.x  = element_text(vjust=0.5, size=12,colour="black"),
          axis.title.x  = element_blank(),
          legend.title = element_blank(),
          legend.position = 'none',
          legend.key = element_rect(fill = NA),
          legend.text = element_text(size=12,colour="black"),
          plot.margin = unit(c(0.5, 0.5, 0,0), "cm"))
}

length(cvRP.list)

# Layout of all plots
cvRP <-  wrap_plots(cvRP.list) + plot_layout(ncol=4)
cvRP

# Save the plot data
save(results,asg.perf,imppl, impplot,pplotdata,cvRP,cvRP.list, 
     file=paste0("data/SDM2024/processed/RF_",rvar,"_Results_Summary.RData"))


#### Raster outputs ----

### Create a raster stack for spatial confidence results
ROutput <- stack()

### Calculate most frequent class and its frequency
# Most frequent class - change to 0 and 1 for absence and presence
MaxClass <- modal(cvpred,freq=FALSE)  # modal() from raster/terra function
ROutput <- addLayer(ROutput,MaxClass)
# Frequency of most frequent class (fraction of runs)
MaxClassF <- modal(cvpred,freq=TRUE)/nruns
ROutput <- addLayer(ROutput,MaxClassF)

### Calculate average probabilities for classes
classsums <- Reduce("+", cvpred.cps)
AvePclass <- classsums / nruns

### Find average probability of maximum frequency class
MaxClassAveProb <- stackSelect(AvePclass, MaxClass+1)
ROutput <- addLayer(ROutput,MaxClassAveProb)

### Calculate new layer for frequency x probability
CombConf <- MaxClassF * MaxClassAveProb
ROutput <- addLayer(ROutput,CombConf)

### number of models predicting presence
cvPA <-  stack(cvpred)
cvSum <-  raster::calc(cvPA,sum)
ROutput <- addLayer(ROutput,cvSum)       

### Rename layers
names(ROutput) <- c("MaxClass","MaxClassF","MaxClassAveProb","CombConf","cvSum")

### Plot layers 
plot(ROutput)

### Export Raster
raster::writeRaster(ROutput$MaxClass, 
                    paste0("output/SDM2024/",rvar,"_raster_output_maxclass.tif"), format="GTiff",overwrite=T)
raster::writeRaster(ROutput$MaxClassF, 
                    paste0("output/SDM2024/",rvar,"_raster_output_maxclassf.tif"), format="GTiff",overwrite=T)
raster::writeRaster(ROutput$MaxClassAveProb, 
                    paste0("output/SDM2024/",rvar,"_raster_output_maxclassaveprob.tif"), format="GTiff",overwrite=T)
raster::writeRaster(ROutput$CombConf, 
                    paste0("output/SDM2024/",rvar,"_raster_output_combconf.tif"), format="GTiff",overwrite=T)
raster::writeRaster(ROutput$cvSum,  
                    paste0("output/SDM2024/",rvar,"VME_raster_output_cvsum.tif"), format="GTiff",overwrite=T)

### Save raster stack to R workspace
save(AvePclass,ROutput,file= paste0("data/SDM2024/processed/",rvar,"_raster_output_all.RData"))


### Extrapolation areas
# library(dsmextra)  # no longer available

covariates.names <- predsel

allpred <- rasterToPoints(predictors)
p.data.all = data.frame(allpred)
aftt_crs <- sp::CRS("+proj=utm +zone=23 +datum=NAD83 +units=m +no_defs")

extrapolation.area <- compute_extrapolation(samples = mdata,
                                            covariate.names = predsel,
                                            prediction.grid = p.data.all,
                                            coordinate.system = aftt_crs)

plot(extrapolation.area$rasters$ExDet$analogue) # analogue areas
plot(extrapolation.area$rasters$ExDet$univariate) # univariate extrapolation
plot(extrapolation.area$rasters$ExDet$combinatorial) # combinatorial extrapolation
plot(extrapolation.area$rasters$mic$analogue) # most important variables causing analogue conditions
plot(extrapolation.area$rasters$mic$univariate) # most important variables causing univariate extrapolation
plot(extrapolation.area$rasters$mic$combinatorial) # most important variables causing combinatorial extrapolation


### Export Raster
writeRaster(extrapolation.area$rasters$ExDet$analogue, 
            paste0("output/SDM2024/",rvar,"_ext.analogue.tif"), format="GTiff",overwrite=T)
writeRaster(extrapolation.area$rasters$ExDet$univariate, 
            paste0("output/SDM2024/",rvar,"_ext.univariate.tif"), format="GTiff",overwrite=T) 
writeRaster(extrapolation.area$rasters$ExDet$combinatorial, 
            paste0("output/SDM2024/",rvar,"_ext.combinatorial.tif"), format="GTiff",overwrite=T) 
writeRaster(extrapolation.area$rasters$mic$analogue, 
            paste0("output/SDM2024/",rvar,"_mic.analogue.tif"), format="GTiff",overwrite=T) 
writeRaster(extrapolation.area$rasters$mic$univariate, 
            paste0("output/SDM2024/",rvar,"_mic.univariate.tif"), format="GTiff",overwrite=T) 
writeRaster(extrapolation.area$rasters$mic$combinatorial, 
            paste0("output/SDM2024/",rvar,"_mic.combinatorial.tif"), format="GTiff",overwrite=T)




