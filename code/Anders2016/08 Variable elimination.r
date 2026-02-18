# 08 Variable elimination
# This R script contains two different ways of eliminating environmental predictor variables
# The first interactively eliminates highly correlated predictors
# The second eliminates predictors without which predictions do not suffer substantially 

#################
# Step 1, setup #
#################

# Import libraries
library(raster)
library(rgdal)
library(maptools)
library(randomForest)
library(pROC)
library(tcltk2)

# Set some directories
setwd("C:/Users/Anders/Dropbox/SDM course/")
rasterdir = "data/rasters/"
csvdir = "data/"

# This is a shortcut if you want to use ALL the raster files (with .tif extensions) in the rasterdir
predictorfiles = list.files(path = paste(rasterdir, sep=""), pattern = "\\.tif$", full.names = F)
predictorfiles = predictorfiles[!predictorfiles == "XLdepth.tif"]

# Now read the raster data
predictors = c()
for(x in predictorfiles)
{
  predname = paste(rasterdir, x, sep="")
  predictors = stack(c(predictors, raster(predname)))
}

#################################################
# Step 2: Eliminate highly correlated variables #
#################################################
# First, set the threshold for maximum allowed correlation between two variables. This is subjective!!
threshold = 0.5

# Unfortunately R does not have a function to directly calculate correlations between rasters, so we need to extract the values to another format
vals = c() # Create empty list for values from the rasters
for (i in 1:nlayers(predictors)) # for each raster...
{
  v = getValues(raster(predictors, i)) # Get the values
  vals = cbind(vals, v) # Add to 'vals'
}
colnames(vals) = names(predictors) # Add names to each list entry
correlations = cor(vals, use="pairwise.complete.obs", method="spearman") # Calculate correlations. I think it makes sense to use Spearman
correlations # Display correlations
for (i in 1:nlayers(predictors)) {correlations[i,i]=0} # Set same-variable correlations to 0
correlations = abs(correlations) # Use absolute values of correlations

# To make the variable elimination interactive, we need to set up some objects for a GUI
# Set up some filenames
tmpfilename = "temporary.txt"
outlistfilename = "elimination_list.txt"
remainingfilename = "remaining_list.csv"

# List of eliminated variables
eliminated = c()

# Find the highest correlation, and the index positions its two variables
maxcor = max(correlations)
maxrow = which(correlations==maxcor,arr.ind=TRUE)[1,1]
maxcol = which(correlations==maxcor,arr.ind=TRUE)[1,2]

# Now comes some code for the GUI widget. You do not need to understand all of this!
# Write some functions for each button
write1 = function()
{
  tmpfile = file(tmpfilename)
  writeLines(colnames(correlations)[maxrow], tmpfile)
  close(tmpfile)
}
write2 = function()
{
  tmpfile = file(tmpfilename)
  writeLines(colnames(correlations)[maxcol], tmpfile)
  close(tmpfile)
}

# Create the GUI widget. Note that I am NOT experienced with GUI coding, so this is probably really poor programming, but it works...
root = tktoplevel()
btn1 = tk2button(root, text=paste(colnames(correlations)[maxrow]), command=write1)
tkpack(btn1)
btn2 = tk2button(root, text=paste(colnames(correlations)[maxcol]), command=write2)
tkpack(btn2)

# Now that the GUI is active, this loop waits for input
# When input is given, it eliminates the variable clikced on, and restarts by showing two new variables
# This continues until no two variables correlate more than the threshold
while (maxcor > threshold)
{
  # Name buttons after the two variables
  tkconfigure(btn1, text=paste(colnames(correlations)[maxrow]), command=write1)
  tkconfigure(btn2, text=paste(colnames(correlations)[maxcol]), command=write2)
  
  # Wait until the results have been written to the temporary file
  while (!file.exists(tmpfilename)) {a=1} # a is a dummy variable
  
  # Read the variable to be eliminated from the temporary file
  tmpfile = file(tmpfilename, "rw")
  eliminate = readLines(tmpfile, 1) # Read one line
  close(tmpfile) # Close the file
  file.remove(tmpfilename) # Delete the file
  
  # Eliminate variable by adding it to the list and removing it from 'correlations'
  eliminated = c(eliminated, eliminate)
  correlations = correlations[-which(colnames(correlations) == eliminate),-which(colnames(correlations) == eliminate)]
  
  # If there's more than one variable left, continue checking  
  if(dim(correlations)[1] > 1)
  {
    maxcor = max(correlations)
    maxrow = which(correlations==maxcor,arr.ind=TRUE)[1,1]
    maxcol = which(correlations==maxcor,arr.ind=TRUE)[1,2]
  }
  
  # Otherwise end
  if(dim(correlations)[1] == 1)
  {
    maxcor = 0
  }
}

# Remove the GUI widget
tkdestroy(root)

# Finally write the results to a file
cornames=colnames(correlations)
write.table(cornames, file=remainingfilename, row.names=FALSE)
outlist = file(outlistfilename)
writeLines(eliminated, outlist)
close(outlist)

##############################################################################
# Step 3: Eliminate variables that don't contribute to improving predictions #
##############################################################################
# We also need a threshold for the jack-knifing, although it is very subjective and imperfect
# This threshold indicates how great a drop in AUC we are willing to allow while still removing a predictor
threshold = 0.01

# Read the response locations and set coordinates
responsefilename = paste(csvdir, "speciesdata.csv", sep="")
responsefile = read.csv(responsefilename)
responsedata = data.frame(species=responsefile$species,x=responsefile$x,y=responsefile$y)
response = responsedata[which(!is.na(responsedata[,"species"])),]
coordinates(response) = ~x+y

# Extract data
p.data = extract(predictors, response)
data = data.frame(response, p.data)
p.names = colnames(p.data)
data = data[complete.cases(data),]
data$species = as.factor(data$species)

# Now we have the data, it is time to do the jack-knifing.
# Although you should ideally do this using cross-validation, we'll use the 70/30 split because it is faster
calsize = 0.70
n = nrow(data)
shuffle = sample(1:n, n)
cal = shuffle[1:(n*calsize)]
val = shuffle[(n*calsize+1):n]
calibration = data[cal,]
validation = data[val,]

# Train the model, and calculate AUC with validation data
formula = as.formula(paste("species","~",paste(p.names, collapse="+")))
model = randomForest(formula, data=calibration, type="prob")
observations = validation$species
predictions = predict(model, newdata=validation, type="prob")[,2]
test.auc = auc(observations, predictions)


jklist = list() #This is where we store the information from the jack-knifing procedure, for reporting later
startlength = length(p.names) # Again, we know this is 50, but deriving it from p.names makes the code more generic
jklist$p.names = array("", startlength, dimnames = "p.names")
jklist$auc = array(0, startlength, dimnames = "Test AUC")
jklist$p.names[1] = paste(p.names, collapse=", ")
jklist$auc[1] = test.auc
jklist

jklist.txt.file = paste("jklist.", threshold, ".txt", sep = "") # Define file to write to
text_to_write = paste(paste(p.names, collapse=","), test.auc) # Create line of text to write, including predictors and AUC
write(text_to_write, file=jklist.txt.file, append=FALSE) # Write line to file

repeat
{ 
  if (length(p.names) == 1) {break} # Break if there's only one predictor left
  
  # Some print statements so we can follow along (this loop takes a while)
  print(cat("Predictors: ", p.names))
  print(cat("Current AUC: ", test.auc))
  
  # Create some data structures to hold results
  jk.auc = array(0, length(p.names), dimnames = "p.names") # Create array to hold AUC values calculated with each predictor removed
  
  for (i in 1:length(p.names))
  {
    jktest.p.names = p.names[-i]
      
    # Fit and predict        
    formula = as.formula(paste("species~", paste(jktest.p.names, collapse = "+")))
    model = randomForest(formula, data=calibration, na.action = na.omit)
    predictions = predict(model, newdata=validation, type="prob")[,2]
    observations = validation$species
  
    jk.auc[i] = auc(observations, predictions)
  }
  
  # Find layer to take out 
  if (max(jk.auc)<(test.auc-threshold)) # If eliminating a predictor results in too poor an AUC value, then break
  {
    break  
  } else { # If eliminating a predictor can improve AUC or at least reduce it by less than the threshols value, eliminate it and repeat
    index = which(jk.auc == max(jk.auc))[1] # Find the predictor whose removal results in the highest AUC
    p.names = p.names[-index] # Remove predictor
    test.auc = max(jk.auc) # Assign new test.auc
    jklist$p.names[startlength-length(p.names)+1] = paste(p.names, collapse=", ") # Remove predictor from p.names so it is not tested again
    jklist$auc[startlength-length(p.names)+1] = test.auc # Add the new AUC value to jklist
  }
  text_to_write = paste(paste(p.names, collapse=","), test.auc) # Create line of text to write, including predictors and AUC
  write(text_to_write, file=jklist.txt.file, append=TRUE) # Write line to file
}

