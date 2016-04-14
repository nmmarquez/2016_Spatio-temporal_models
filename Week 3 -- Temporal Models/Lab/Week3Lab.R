# Spatio-temporal models for ecologists
# Week 3 Lab
# Temporal models
# 


######
# PART 1 - TIME SERIES ANALYSIS
#####
# CO2 concentrations from Mauna Loa Observatory, Hawaii 3/1952 - 2/2016
# mole fraction of CO2, expressed as parts per million
# (ppm) is the number of molecules of CO2 in every one million molecules of dried
# air (water vapor removed)
# See www.esrl.noaa.gov/gmd/ccgg/trends/ for additional details
#

co2.data<- read.table("CO2.csv", sep = ',', header = T)
all.ts<-ts(co2.data$CO2, start = c(1958, 3), frequency = 12)
plot(all.ts, ylab = expression(paste(CO[2], " ppm") ) )

#subset data to end of 2012
co2.ts<- window(all.ts, start = c(1958, 3), end = c(2012, 12) )

#Decompose time series into trend, seasonal, and remainder 
#very simple model: trend is exponential, seasonal component is consistent across years
co2.stl<- stl(log(co2.ts), s.window = "periodic", robust = TRUE, t.window = 1000) 
plot(co2.stl)
summary(co2.stl)

#obtain remainder 
co2.rem<-co2.stl$time.series[,3]

#Evaluate remainder for AR and MA structure
acf(co2.rem)  #does not suggest MA
acf(co2.rem, type = "partial")
co2.ar<- ar(co2.rem, method = "mle")
plot(0:(length(co2.ar$aic) - 1), co2.ar$aic, xlab = "Order", ylab = "AIC")

#fit the AR model
co2.ar.model<- arima(co2.rem, order = c(9,0,0), method = "ML")

#Evaluate residuals 
library("lmtest")
par(mfrow = c(1,2))
hist(co2.ar.model$residuals)
qqnorm(co2.ar.model$residuals)
#Durbin-Watson test for autocorrelation:
dwtest(co2.ar.model$residuals~1)

#Predict series from the seasonal decomposition and AR
#AR component 3 years into the future - predicted + SE
co2.ar.for<- predict(co2.ar.model, n.ahead = 36)

#Add AR to seasonal and trend components from stl
co2.trend.slope<- co2.stl$time.series[2,2] - co2.stl$time.series[1,2]
co2.trend.for<- 1:36 * co2.trend.slope + co2.stl$time.series[658,2]
co2.season.trend.for<- co2.trend.for + rep(co2.stl$time.series[647:658,1], times = 3)
co2.mean<- exp(co2.season.trend.for + co2.ar.for$pred)
co2.low<- exp(co2.season.trend.for + co2.ar.for$pred - 1.96*co2.ar.for$se)
co2.hi<- exp(co2.season.trend.for + co2.ar.for$pred + 1.96*co2.ar.for$se)

co2.obs<- window(all.ts, start = c(2013, 1), end = c(2015, 12) )

#Plot predictions and observed data
plot(co2.mean, ylim = c(390, 405), ylab = "CO2 ppm" )
lines(co2.low, lty = 3)
lines(co2.hi, lty = 3)
lines(co2.obs, col = 2)
legend("topleft", c("Model", "Observed"), lty = c(1,1), col = c(1,2))



######
# PART 2 - STATE SPACE MODELS
#####

library(TMB)

compile("dlm.cpp")
dyn.load(dynlib("dlm"))

DLM_sim <- function(N = 100, seed = 123, sigma.obs = 1, sigma.proc = 0.5, a = -0.5, y1 = 10){

#set.seed(seed)
N <- 100
ytrue <- numeric(length = N)
ytrue[1] <- y1
log.sigma.proc <- log(sigma.proc)
proc.error <- rnorm(N, mean = 0, sd = sigma.proc)
log.sigma.obs <- log(sigma.obs)
for(i in 2:N) {
    ytrue[i] <- a*ytrue[i-1] + proc.error[i-1]
    }
  x <- seq_len(N)
  y <- rnorm(N, mean = ytrue, sd = sigma.obs)
return(list(y = y, ytrue = ytrue) )
}

sim<- DLM_sim(a= 0.5, N=100, sigma.obs = 1, sigma.proc = 0.5, y1 = 3)

plot(sim$ytrue, type = 'b', col = 2, pch = 's', ylim = c(min(sim$y), max(sim$y)))
lines(sim$y)
points(sim$y, pch = 'o')

#Run TMB model on simulated data
  N = length(sim$y)	
  data <- list(y = sim$y)
  parameters <- list(a = 0.5, log_sigma_proc = -1,
  log_sigma_obs = -1, u = rep(mean(sim$y), N))
  obj <- MakeADFun(data, parameters, random = "u", DLL = "dlm")
  obj$hessian <- FALSE
  opt <- do.call("optim", obj)
  sd <- sdreport(obj)

  # extract fixed effects:
  fixed <- summary(sd, "fixed")
  
  # extract estimated process:
  u <- summary(sd, "random")[, "Estimate"]
  u_se <- summary(sd, "random")[, "Std. Error"]
  
 
  


  
