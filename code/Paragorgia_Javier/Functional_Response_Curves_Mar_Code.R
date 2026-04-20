library(ggplot2)

#For the full study area -JM code modification based on Mar Sacau code
allpred<-rasterToPoints(predictors)
dim(allpred)
head(allpred)
env = data.frame(allpred)
dim(env)
env = env[complete.cases(env),]
dim(env)
head(env)
xy=cbind(env$x,env$y)
head(xy)

# resolution from predictors
wdir = "E:/2021_April_May_lockdown_JM/Publications/Paragorgia_SDM_paper_2021_Shuangqiang/Results/RF_predictions"
setwd(wdir)

# RF present model
rf_probability = raster("rf_Paragorgia_present_8var.tif")
plot(rf_probability)
pred<-extract(rf_probability, xy)
length(pred)
data.response=as.data.frame(cbind(env,pred))
head(data.response)
dim(data.response)


windows()
qplot(data.response$slope,data.response$pred,geom="smooth", span=0.6, xlab=deparse(substitute(Slope)),ylab=(deparse(substitute(Probability))))
qplot(env1$asp_ns_clip,env1$pred1,geom="smooth",span=0.6,ylim=c(0,1),xlab=deparse(substitute(Aspect_NS)), ylab=(deparse(substitute(Probability))))
qplot(env1$speed_clip,env1$pred1,geom="smooth",span=0.6,ylim=c(0,1),xlab=deparse(substitute(Bottom_current_speed)), ylab=(deparse(substitute(Probability))))


# GAM present model
wdir = "E:/2021_April_May_lockdown_JM/Publications/Paragorgia_SDM_paper_2021_Shuangqiang/Results/GAM_predictions"
setwd(wdir)

# GAM present model
gam_probability = raster("GAM_test_present.tif")
plot(gam_probability)
pred<-extract(gam_probability, xy)
length(pred)
data.response=as.data.frame(cbind(env,pred))
head(data.response)
dim(data.response)


##### FIGURES

setwd("C:/Users/murillo-perezj/Desktop/Functional_Response_Curves/GAM")


pdf("RF_slope.pdf", family="Times", width=10, height=8)


p <- ggplot()

# plot the points
#p <- p + geom_point(data = data.response, 
#                    aes(y = pred, x = slope),
#                    shape = 16, 
#                    size = 3)

# plot the curve
p <- p + xlab("Slope") + ylab("Probability")
p <- p + theme(text = element_text(size=30)) 
p <- p + geom_smooth(data = data.response,  
                     aes(x =  slope, 
                         y = pred), span = 0.6, col="red")

p <- p + xlab("Bottom Salinity") + ylab("Probability")
p <- p + theme(text = element_text(size=30)) 
p <- p + geom_smooth(data = data.response,  
                     aes(x =  b_sal_mean_NEMOKrig, 
                         y = pred), span = 3, col="red")

p <- p + theme_classic()
p <- p + theme(text = element_text(size=40)) 
p <- p + theme(axis.text.x = element_text(colour = "black"),
               axis.text.y = element_text(colour = "black"))

p
# smoothing method (function) to use, eg. "lm", "glm", "gam", "loess", "rlm". For method = "auto"
# the smoothing method is chosen based on the size of the largest group (across all panels).
# loess is used for than 1,000 observations; otherwise gam is used with formula = y ~ s(x, bs = "cs").
# Somewhat anecdotally, loess gives a better appearance, but is O(n^2) in memory,
# so does not work for larger datasets.

dev.off()

#Only for the observations
head(data)
dim(data)
xy1=cbind(data$x,data$y)

# RF present model
pred1<-extract(rf_probability, xy1)
length(pred1)
data.response.obs=as.data.frame(cbind(data,pred1))
head(data.response.obs)

# GAM present model
pred1<-extract(gam_probability, xy1)
length(pred1)
data.response.obs=as.data.frame(cbind(data,pred1))
head(data.response.obs)

pdf("GAM_b_cur_mean_obs.pdf", family="Times", width=10, height=8)
pdf("GAM_b_sal_mean_obs.pdf", family="Times", width=10, height=8)
pdf("GAM_b_tmp_mean_obs.pdf", family="Times", width=10, height=8)
pdf("GAM_MLD_obs.pdf", family="Times", width=10, height=8)
pdf("GAM_s_cur_mean_obs.pdf", family="Times", width=10, height=8)
pdf("GAM_s_sal_mean_obs.pdf", family="Times", width=10, height=8)
pdf("GAM_s_tmp_mean_obs.pdf", family="Times", width=10, height=8)
pdf("GAM_slope_obs.pdf", family="Times", width=10, height=8)

p <- ggplot()

# plot the points
p <- p + geom_point(data = data.response.obs, 
                    aes(y = pred1, x = s_tmp_mean_NEMOKrig),
                    shape = 16, 
                    size = 3)

# plot the curve

# Bottom Current
p <- p + xlab("Bottom Current") + ylab("Probability")
p <- p + theme(text = element_text(size=30)) + scale_y_continuous(limits = c(0, NA))
p <- p + geom_smooth(data = data.response.obs,  
                     aes(x =  b_cur_mean_NEMOKrig, 
                         y = pred1), span = 0.6, col="red", se = FALSE)


# Bottom Salinity
p <- p + xlab("Bottom Salinity") + ylab("Probability")
p <- p + theme(text = element_text(size=30)) + scale_y_continuous(limits = c(0, NA))
p <- p + geom_smooth(data = data.response.obs,  
                     aes(x =  b_sal_mean_NEMOKrig, 
                         y = pred1), span = 0.6, col="red")
# Bottom Temperature
p <- p + xlab("Bottom Temperature") + ylab("Probability")
p <- p + theme(text = element_text(size=30)) + scale_y_continuous(limits = c(0, NA))
p <- p + geom_smooth(data = data.response.obs,  
                     aes(x =  b_tmp_mean_NEMOKrig, 
                         y = pred1), span = 0.6, col="red")

#to read how to plot only line and error separately
#https://aosmith.rbind.io/2018/11/16/plot-fitted-lines/

  # to set limits
#https://stackoverflow.com/questions/11214012/set-only-lower-bound-of-a-limit-for-ggplot
  

# MLD
p <- p + xlab("Mixed Layer Depth") + ylab("Probability")
p <- p + theme(text = element_text(size=30)) + scale_y_continuous(limits = c(0, NA))
p <- p + geom_smooth(data = data.response.obs,  
                     aes(x =  MLD_present_Beazley_et_al_2021, 
                         y = pred1), span = 0.6, col="red")
# Surface Current
p <- p + xlab("Surface Current") + ylab("Probability")
p <- p + theme(text = element_text(size=30)) + scale_y_continuous(limits = c(0, NA))
p <- p + geom_smooth(data = data.response.obs,  
                     aes(x =  s_cur_mean_NEMOKrig, 
                         y = pred1), span = 0.6, col="red")
# Surface Salinity
p <- p + xlab("Surface Salinity") + ylab("Probability")
p <- p + theme(text = element_text(size=30)) + scale_y_continuous(limits = c(0, NA))
p <- p + geom_smooth(data = data.response.obs,  
                     aes(x =  s_sal_mean_NEMOKrig, 
                         y = pred1), span = 0.6, col="red")
# Surface Temperature
p <- p + xlab("Surface Temperature") + ylab("Probability")
p <- p + theme(text = element_text(size=30)) + scale_y_continuous(limits = c(0, NA))
p <- p + geom_smooth(data = data.response.obs,  
                     aes(x =  s_tmp_mean_NEMOKrig, 
                         y = pred1), span = 0.6, col="red")
# slope
p <- p + xlab("Slope") + ylab("Probability")
p <- p + theme(text = element_text(size=30)) + scale_y_continuous(limits = c(0, NA))
p <- p + geom_smooth(data = data.response.obs.lm,  
                     aes(x =  slope, 
                         y = pred1), span = 0.6, col="red") 

p <- p + theme_classic()
p <- p + theme(text = element_text(size=40)) 
p <- p + theme(axis.text.x = element_text(colour = "black"),
               axis.text.y = element_text(colour = "black"))
p

dev.off()

qplot(data.response.obs$slope, data.response.obs$pred1,geom="smooth", span=0.6, xlab=deparse(substitute(Slope)),ylab=(deparse(substitute(Probability))))



### to review in another moment

gam.gg = gam(pred1 ~ s(slope, bs = "cs"), data = data.response.obs)
data.response.obs$predgam.gg = predict(gam.gg)
head(data.response.obs)

ggplot(data.response.obs, aes(x = slope, y = pred1)) +
  geom_line(aes(y = predgam.gg), size = 1)

predslm = predict(gam.gg, interval = "confidence")

preds   <- predict(gam.gg, se.fit = TRUE)
head(preds)
predslm <- data.frame(mu   = exp(preds$fit),
                      low  = exp(preds$fit - 1.96 * preds$se.fit),
                      high = exp(preds$fit + 1.96 * preds$se.fit))

head(predslm)
data.response.obs.lm = cbind(data.response.obs, predslm)
head(data.response.obs.lm)


ggplot(data.response.obs.lm, aes(x =  slope, 
                         y = pred1), span = 0.6, col="red") +
  geom_ribbon(data = data.response.obs.lm, aes(ymin = low, ymax = high, color = NULL), alpha = .15)
  

