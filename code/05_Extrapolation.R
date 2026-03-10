
# Extrapolations with dsmextra

covariates.names <- predsel

allpred <- rasterToPoints(predictors)
p.data.all = data.frame(allpred)
aftt_crs <- sp::CRS("+proj=utm +zone=23 +datum=NAD83 +units=m +no_defs")

extrapolation.area <- dsmextra::compute_extrapolation(samples = mdata,
  covariate.names = predsel,
  prediction.grid = p.data.all,
  coordinate.system = aftt_crs)

plot(extrapolation.area$rasters$ExDet$analogue) # analogue areas
plot(extrapolation.area$rasters$ExDet$univariate) # univariate extrapolation
plot(extrapolation.area$rasters$ExDet$combinatorial) # combinatorial extrapolation
plot(extrapolation.area$rasters$mic$analogue) # most important variables causing analogue conditions
plot(extrapolation.area$rasters$mic$univariate) # most important variables causing univariate extrapolation
plot(extrapolation.area$rasters$mic$combinatorial) # most important variables causing combinatorial extrapolation