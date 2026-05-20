#Load packages and library files
# library(raster)
# library(rgdal)
# library(maptools)
# library(randomForest)
# library(ranger)
# library(pROC)

######################################################################
###### STEP 0: IMPORT THE RESPONSE AND ENVIRONMENTAL DATA ############
######################################################################
wdir = ("D:/Murillo_WorkSince2015/Publicaciones/Publications/Paragorgia_SDM_paper_2021_Shuangqiang/dsmextra/Telmo_Paper_rasters")
wdir = ("D:/Murillo_WorkSince2015/Publicaciones/Publications/Paragorgia_SDM_paper_2021_Shuangqiang/Data")

setwd(wdir)

# Environmental data
rasterdir = "/Predictors/Predictors_RF_Paragorgia_SDM_Murillo&Kenchington2020Oct19"
rasterdir = "/Predictors/Predictors_RF_Paragorgia_SDM_Murillo&Kenchington2020Oct19_top6"
rasterdir = "/predictors/predictors_used/present"   
rasterdir = "/predictors/predictors_used/RCP45_Resample_names_changed"
rasterdir = "/predictors/predictors_used/RCP85_Resample_names_changed"
rasterdir = "/predictors/predictors_used/Present"  
rasterdir = "/Predictors/Telmo_Paper_rasters_full_extent/Future"  
rasterdir = "/Future"  


predictorfiles = list.files(path = paste(wdir, rasterdir, sep=""), pattern = "\\.tif$", full.names = F)
predictorfiles

# Now read the raster data (create a raster stack)
predictors = c()  #creates an empty vector used as the raster stack
for(x in predictorfiles)  #loop through all entries in 'predictorfiles'
{
  predname = paste(wdir,rasterdir,"/",x,sep="")  #generate a complete file path to 'x'
  predictors = stack(c(predictors,raster(predname)))  #add raster (predname) to stack
}
predictors #confirms raster stack with all raster layers present


# Response
wdir = ("D:/Murillo_WorkSince2015/Publicaciones/Publications/Paragorgia_SDM_paper_2021_Shuangqiang/dsmextra/Telmo_Paper_rasters")
wdir = ("D:/Murillo_WorkSince2015/Publicaciones/Publications/Paragorgia_SDM_paper_2021_Shuangqiang/Data/Response")
setwd(wdir)
responsefile = read.csv ("AA_Lg_Gor_2020_v3_Paragorgia_1datacell_no_Gass_JM.csv", header=TRUE)
responsefile = read.csv ("Paragorgia_arborea_And_PseudoAbs_3km_CorrectedByEffort_JM.csv", header=TRUE)
head(responsefile)
dim(responsefile)

responsedata = data.frame(pa=responsefile$v,x=responsefile$POINT_X, y=responsefile$POINT_Y)
head(responsedata)
dim(responsedata)



response = responsedata[which(!is.na(responsedata[,"pa"])),] 

# Set coordinates (latitude and longitude)
coordinates(response) = ~x+y #converts response from a normal dataframe to a "SpatialPointsDataFrame" by telling R
dim(response)
head(response)

#########################################################################
#### STEP 2. EXTRACT THE VALUES OF PREDICTORS AT RESPONSE LOCATIONS #####
#########################################################################

# Extract the values from each predictor for each location in the response data, and put results in a new dataframe
p.data = extract(predictors, response)
data = data.frame(response, p.data) #adds all the extracted values to the existing dataframe
head(data)
p.names = colnames(p.data)  #get names from new columns
str(data)

#########################################################################
#### STEP 3. QUALITY CONTROL ############################################
#########################################################################

# Ensure that all observations (locatiions) have data values for all variables
data = data[complete.cases(data),]  #'complete.cases' command returns only those rows in the dataframe that have non-NA 
head(data)
dim(data)
#values for all columns
str(data)
nrow(data)
missingdata = data.frame(responsedata, p.data)  #creates a new dataframe called 'missingdata'
missingdata = missingdata[!complete.cases(missingdata),]  #!interested in rows that are NOT complete cases
dim(missingdata)

#########################################################################
#### STEP 4. BUILD THE RANDOM FOREST MODEL ##############################
#########################################################################

# Build the model
formula = as.formula(paste("pa","~",paste(p.names, collapse="+")))

# Warning "Are you sure you want to do regression?"
typeof(data$pa)

# Build the model again, after converting PresenceAbsence to a categorical (factor) variable
data$pa = as.factor(data$pa)

########## USING RANDOM FOREST PACKAGE ##########

model = randomForest(formula, data=data)
model

#########################################################################
#### STEP 5. PREDICT FOR ALL CELLS IN RASTER EXTENT######################
#########################################################################

# Create a map that that shows the probability of presence (from 0 to 1) for each cell in the raster extent
model_FINAL = randomForest(formula, data=data, ntree = 500, type="prob")
map = predict(predictors, model_FINAL, type="prob", index=2)
outfile = "rf_Paragorgia_47cov_updated_115presences.tif"
writeRaster(map, filename=outfile, overwrite=TRUE, format="GTiff", datatype="FLT4S")
plot(map)

wdir = "D:/Murillo_WorkSince2015/Publicaciones/Publications/Paragorgia_SDM_paper_2021_Shuangqiang/Results/RF_predictions"
setwd(wdir)


################################################################
#### STEP 4&5 USING RANGER PACKAGE #############################
################################################################

# ranger cannot predict using the predictors in rasters
allpred<-rasterToPoints(predictors)
dim(allpred)
head(allpred)
new.data = data.frame(allpred)
new.data = new.data[complete.cases(new.data),]

head(new.data)
dim(new.data)
coordinates_area = new.data[,c(1,2)]
head(coordinates_area)

model_ranger = ranger(formula, data=data, num.trees = 500, keep.inbag=TRUE)
response = predict(model_ranger, new.data, type = "response")
Gx_area_map = cbind(coordinates_area, response$predictions)
head(Gx_area_map)
dim(Gx_area_map)

error.jack = predict(model_ranger, new.data, type = "se", se.method = 'jack')
error.infjack = predict(model_ranger, new.data, type = "se", se.method = 'infjack')
Gx_area_map = cbind(coordinates_area, error.jack$se)
mean(error.jack$se)
Gx_area_map = cbind(coordinates_area, error.infjack$se)
mean(error.infjack$se)

#################################################################
#### STEP 6. IMPORTANCE #########################################
#################################################################
#n = nrow(data)
runs = 10

for (run in 1:runs)
{
  
  model = randomForest(formula, data=data, ntree = 500, type="prob", importance=TRUE)
  #save variable importances to dataframe for averaging after loop
  importance = data.frame(model$importance[,4])
  if(run==1){imp.df = data.frame(row.names = rownames(importance))}
  imp.df = cbind(imp.df, importance)
  
}#end of runs


imp = data.frame(row.names=rownames(imp.df), Means = rowMeans(imp.df))
imp


impmatrix<-as.matrix(imp, ncol=2)
sortimpmatrix<-sort(impmatrix[,1], decreasing=TRUE)
sortimpmatrix
top15<-sortimpmatrix[1:15]
top15matrix<-as.matrix(top15)
top15matrix<-top15matrix[order(top15matrix[,1]),]

par(mar=c(5.1, 7.0, 3,3))
mpg=c(5,1,0)
dotchart(top15matrix, color='black', pch=16, xlab="Mean Decrease in Gini Value")   #Save as PNG to avoid label overlap with y axis

# 56 predictors
names(top15)<-c("Ruggedness",
                "Positive Topographic Openness",
                "LS-Factor",
                "Slope",
                "Negative Topographic Openness",
                "Valley Depth",
                "Channel Network Base Level",
                "Surface Temperature Maximum",
                "Bottom Salinity Minimum",
                "Fine-Scale Bathymetric Position Index",
                "Wind Exposition Index",
                "Spring Chlorophyll a Minimum",
                "Bottom Current Maximum",
                "Annual Primary Production Mean",
                "Broad-Scale Bathymetric Position Index")


top15matrix<-as.matrix(top15)
top15matrix<-top15matrix[order(top15matrix[,1]),]

par(mar=c(5.1, 7.0, 3,3))
mpg=c(5,1,0)

library(extrafont)

setwd("C:/Users/DELL/Desktop/Murillo_NewTerm_2020/Publications/Paragorgia_SDM/Results/Figures/Final")

pdf("Importance_FINAL_RF_Model_v2.pdf", family="Times", width=7, height=5)

dotchart(top15matrix, color='black', pch=16, xlab="Mean Decrease in Gini Value") 

dev.off()


#################################################################
#### STEP 7. VALIDATION #########################################
#################################################################

########UNBALANCED DATA#########

### 10-Fold Cross-Validation to Test Model Accuracy ### 

# In cross-validation, all data are used in turn as both calibration and validation data, unlike in manual 70/30 split
# where the model is trained only on 70% of data, and tested with remaining 30%. Considered a more optimal way to test
# model accuracy.
set.seed(2)
n = nrow(data) 
# 10-fold cross-validation
nfold = 10 # Set number of 'folds' for cross-validation 
split = sample(rep(1:nfold, length = n), n) # Split n in 10 sets
resample = lapply(1:nfold, function(x,spl) list(cal=which(spl!=x), val=which(spl==x)), spl=split) # Create values for 
#10 cal/val sets. 
str(resample)  #shows resample is a 'list of 10'

# Create a for loop that takes each combination of calibration and validation data, trains a model using calibration data,
# uses that model to make predictions for the validation data, then calculates the AUC by comparing the predictions to 
# the observations
results = c() # This creates an empty vector to hold the AUC results.
TP=c()
FP=c()
TN=c()
FN=c()
# Here we loop through all cal/val combinations
for(fold in 1:nfold)
{
  # Create the cal/val data sets for the current 'fold'
  cal = data[resample[[fold]]$cal,]
  val = data[resample[[fold]]$val,]
  
  # Fit a model with the current calibration data set
  model = randomForest(formula, data=cal, type="prob") # Note that we are reusing the formula already defined
  
  # Calculate AUC and put it in results
  observations = val$pa
  predictions = predict(model, newdata=val, type="prob")[,2]
  
  if (length(levels(factor(observations))) == length(levels(observations)))
  {
    auc = auc(observations, predictions)   
    results = c(results, auc) # This attaches 'auc' to any previous content of 'results'
    
    # How to creat a confusion matrix from the validated data
    categorical_predictions = predictions
    threshold = 0.051
    categorical_predictions [which(categorical_predictions >= threshold)] = 1  #predictions > 0.5 are presences
    categorical_predictions [which(categorical_predictions < threshold)] = 0  #predictions < 0.5 are absences
    fakeData <- data.frame(obs = factor(observations), pred= factor(categorical_predictions))
    
    # how to create a confusion matrix from a data frame of observations and predictions
    negatives = subset(fakeData, fakeData$obs == 0) # select the negatives
    for (i in 1:ncol(negatives)){negatives[,i]=as.vector(negatives[,i])} # change the columns from
    # factor to vector
    true.negative = ifelse(negatives$obs==negatives$pred,1,0) # get true negatives 
    TN_loop = sum(true.negative)  # sum all the true negatives
    TN = c(TN,TN_loop) # This attaches 'TN_loop' to any previous content of 'TN'
    false.positive = ifelse(negatives$obs!=negatives$pred,1,0) # get false positives
    FP_loop = sum(false.positive) # sum all the false positives
    FP = c(FP,FP_loop) # This attaches 'FP_loop' to any previous content of 'FP'
    
    positives = subset(fakeData, fakeData$obs == 1) # select the postives
    for (i in 1:ncol(positives)){positives[,i]=as.vector(positives[,i])} # change the columns from
    # factor to vector to be able to apply the ifelse function
    true.positive = ifelse(positives$obs==positives$pred,1,0) # get true positives 
    TP_loop = sum(true.positive) # sum all the true positives
    TP = c(TP,TP_loop) # This attaches 'TP_loop' to any previous content of 'TP'
    false.negative = ifelse(positives$obs!=positives$pred,1,0) # get false negatives
    FN_loop = sum(false.negative) # sum all the false negatives
    FN = c(FN,FN_loop) # This attaches 'FN_loop' to any previous content of 'FN'
  }
}

auc = mean(results)
auc

confusion.matrix = matrix(data=c(sum(TN),sum(FN),sum(FP),sum(TP)), ncol=2,nrow=2)

Class.error.negatives = sum(FP)/(sum(TN)+sum(FP))
Class.error.positives = sum(FN)/(sum(TP)+sum(FN))

Accuracy = (sum(TP)+sum(TN))/(sum(TP)+sum(TN)+sum(FN)+sum(FP))
Sensitivity = sum(TP)/(sum(TP)+sum(FN))
Specifity = sum(TN)/(sum(FP)+sum(TN))
TSS = Sensitivity + Specifity - 1


# IF the error 'no case observation' is returned after this code above: more absences than presences, random splitting 
# causes the validation data to only contain absences. Therefore, the AUC can NOT be calculated. See Anders' documents
# for work-around

# Accuracy assessment above only quantifies how good model is at making predictions for areas that are close to those 
# samples that are entered in the calibration dataset

#################################################################
#### STEP 8. PARTIAL PLOTS ######################################
#################################################################

impvar = rownames(imp)[order(imp[,1], decreasing=TRUE)]
out <- list()
for (i in seq_along(impvar)) {
  Test <- NULL
  Test <- partialPlot(model_FINAL, data, impvar[i], xlab=impvar[i],
                      main=paste("Partial Dependence on", impvar[i]), which.class="1")
  out[[impvar[i]]] <- Test
}


# ruggedness_17_geo
ruggedness_17_geoPP<-out[["ruggedness_17_geo"]]
transformy<- exp(ruggedness_17_geoPP$y) / (1 + exp(ruggedness_17_geoPP$y))
ruggedness_17_geo_Prob<-cbind(ruggedness_17_geoPP$x, ruggedness_17_geoPP$y, transformy)
ruggedness_17_geodata<-as.data.frame(ruggedness_17_geo_Prob)
colnames(ruggedness_17_geodata)<-c("x", "y", "transformy")
head(ruggedness_17_geodata)

pdf("PartPlot_Ruggedness.pdf", family="Times", width=5, height=5)

plot(ruggedness_17_geodata$x, ruggedness_17_geodata$transformy, type="n", xlab="Ruggedness", ylab="Presence Probability",cex.axis=1, cex.lab =1.2)
lines(ruggedness_17_geodata$x, ruggedness_17_geodata$transformy, col="black",lwd=3)
par(new=F)

dev.off()


# positive_openness_geo
positive_openness_geoPP<-out[["positive_openness_geo"]]
transformy<- exp(positive_openness_geoPP$y) / (1 + exp(positive_openness_geoPP$y))
positive_openness_geo_Prob<-cbind(positive_openness_geoPP$x, positive_openness_geoPP$y, transformy)
positive_openness_geodata<-as.data.frame(positive_openness_geo_Prob)
colnames(positive_openness_geodata)<-c("x", "y", "transformy")
head(positive_openness_geodata)

pdf("PartPlot_PositiveTopographicOpenness.pdf", family="Times", width=5, height=5)

plot(positive_openness_geodata$x, positive_openness_geodata$transformy, type="n", xlab="Positive Topographic Openness (radians)", ylab="Presence Probability",cex.axis=1, cex.lab =1.2)
lines(positive_openness_geodata$x, positive_openness_geodata$transformy, col="black",lwd=3)
par(new=F)

dev.off()


# slope
slopePP<-out[["slope"]]
transformy<- exp(slopePP$y) / (1 + exp(slopePP$y))
slope_Prob<-cbind(slopePP$x, slopePP$y, transformy)
slopedata<-as.data.frame(slope_Prob)
colnames(slopedata)<-c("x", "y", "transformy")
head(slopedata)

pdf("PartPlot_Slope.pdf", family="Times", width=5, height=5)

plot(slopedata$x, slopedata$transformy, type="n", xlab="Slope (�)", ylab="Presence Probability", cex.axis=1, cex.lab =1.2)
lines(slopedata$x, slopedata$transformy, col="black",lwd=3)
par(new=F)

dev.off()

# LS_factor_5_geo
LS_factor_5_geoPP<-out[["LS_factor_5_geo"]]
transformy<- exp(LS_factor_5_geoPP$y) / (1 + exp(LS_factor_5_geoPP$y))
LS_factor_5_geo_Prob<-cbind(LS_factor_5_geoPP$x, LS_factor_5_geoPP$y, transformy)
LS_factor_5_geodata<-as.data.frame(LS_factor_5_geo_Prob)
colnames(LS_factor_5_geodata)<-c("x", "y", "transformy")
head(LS_factor_5_geodata)

pdf("PartPlot_LS_factor.pdf", family="Times", width=5, height=5)

plot(LS_factor_5_geodata$x, LS_factor_5_geodata$transformy, type="n", xlab="LS-Factor", ylab="Presence Probability", cex.axis=1, cex.lab =1.2)
lines(LS_factor_5_geodata$x, LS_factor_5_geodata$transformy, col="black",lwd=3)
par(new=F)

dev.off()

# negative_openness_geo
negative_openness_geoPP<-out[["negative_openness_geo"]]
transformy<- exp(negative_openness_geoPP$y) / (1 + exp(negative_openness_geoPP$y))
negative_openness_geo_Prob<-cbind(negative_openness_geoPP$x, negative_openness_geoPP$y, transformy)
negative_openness_geodata<-as.data.frame(negative_openness_geo_Prob)
colnames(negative_openness_geodata)<-c("x", "y", "transformy")
head(negative_openness_geodata)

pdf("PartPlot_NegativeTopographicOpenness.pdf", family="Times", width=5, height=5)

plot(negative_openness_geodata$x, negative_openness_geodata$transformy, type="n", xlab="Negative Topographic Openness (radians)", ylab="Presence Probability", cex.axis=1, cex.lab =1.2)
lines(negative_openness_geodata$x, negative_openness_geodata$transformy, col="black",lwd=3)
par(new=F)

dev.off()


# valley_depth_5_geo
valley_depth_5_geoPP<-out[["valley_depth_5_geo"]]
transformy<- exp(valley_depth_5_geoPP$y) / (1 + exp(valley_depth_5_geoPP$y))
valley_depth_5_geo_Prob<-cbind(valley_depth_5_geoPP$x, valley_depth_5_geoPP$y, transformy)
valley_depth_5_geodata<-as.data.frame(valley_depth_5_geo_Prob)
colnames(valley_depth_5_geodata)<-c("x", "y", "transformy")
head(valley_depth_5_geodata)

pdf("PartPlot_ValleyDepth.pdf", family="Times", width=5, height=5)

plot(valley_depth_5_geodata$x, valley_depth_5_geodata$transformy, type="n", xlab="Valley Depth (m)", ylab="Presence Probability", cex.axis=1, cex.lab =1.2)
lines(valley_depth_5_geodata$x, valley_depth_5_geodata$transformy, col="black",lwd=3)
par(new=F)

dev.off()

#################################################################
#### STEP 9. AREAS OF EXTRAPOLATION #############################
#################################################################

# First create an empty raster
extrapolated = raster(predictors, 1) # Make a copy of the first raster in predictors
extrapolated[] = rep(NA, ncell(extrapolated)) # Set all its values to NA

# Then set cells to 1 where extrapolation has occurred
for (i in 1:nlayers(predictors)) # for each predictor...
{
  print (i)
  variable = names(predictors)[i] # Find the name of the i'th environmental variable in our data
  min_observed = min(data[[variable]]) # Find the minimum observed value
  max_observed = max(data[[variable]]) # Find the maximum observed value
  
  r = raster(predictors, i) # Extract the i'th raster from 'predictors'
  extrapolated[r < min_observed] = 1 # For cells where r is less than min_observed, set extrapolated to 1
  extrapolated[r > max_observed] = 1 # For cells where r is greater than max_observed, set extrapolated to 1
}

# Save result to file
outfile = "extrapolated_area_Paragorgia_present_v2.tif"
writeRaster(extrapolated, filename=outfile, overwrite=TRUE, format="GTiff", datatype="INT2S")


