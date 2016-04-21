#Spatiotemporal models for ecologists
#Week 4 Lab
#Geostatistical models

#generate a spatial data set that has a trend and spatial autocorrelation
#trend has a spatially correlated covariate

require(geoR)
require(RandomFields)
require(ape)
require(RColorBrewer)
require(spdep)

set.seed(123)  #so we all have similar results


###########
## In RandomFields GENERATE SPATIALLY AUTOCORRELATED DATA AND THE RANDOM FIELD FROM WHICH IT CAME
#generate locations of observations

# first we simulate some random sampling locations
n <- 200
x <- runif(n=n, min=-1, max=1)
y <- runif(n=n, min=-1, max=1)


#Generate a covariate field, for example depth relative to mean high tide,  - changes linearly with location

#calculate the grid locations
p.grid<- expand.grid(seq(-1, 1, l = 51), seq(-1,1, l = 51) )
#calculate covariate value across the grid
grid.Cov<- 0.2*p.grid[,1] -0.4*p.grid[,2]
#construct a matrix for using image()
grid.Mat<- outer(seq(-1, 1, l = 51), seq(-1, 1, l = 51), FUN = function(x,y) 0.2*x -0.4*y)

par(mfrow = c(1,1))
image(z=grid.Mat,x = seq(-1, 1, l = 51),y =  seq(-1, 1, l = 51), col = brewer.pal(9, "Greens"), xlab = "X Coor", ylab = "Y Coor")
points(x,y,pch = '+')

#Calculate the covariate value at each of the sampling locations
X.cov<- 0.2*x -0.4*y

#generate an autocorrelated set of random effects 
re.cov.model<- RMexp(var = 2.5, scale = 0.4)
re.data <- RFsimulate(model = re.cov.model, x=x, y=y, grid=FALSE)
#to change the color of the plot
RFpar(col= brewer.pal(9, "Blues") )
plot(re.data)

x.seq.cond <- y.seq.cond <- seq(-1.5, 1.5, length=n)

# simulate a field conditional on the above data
re.cond <- RFsimulate(re.cov.model, x=x.seq.cond, y=y.seq.cond, data=re.data)
plot(re.cond, re.data)


#CONSTRUCT A DATA SET WITH CORRELATED COVARIATE

##Generate a spatially autocorrelated mean based on the covariate
#underlying trend
beta0<- 0.3
beta1<- 0.8
y.trend<- beta0 + beta1*X.cov
#underlying mean with spatially autocor random effects
y.mean<- y.trend + re.data$variable1
#observed data
y.obs<- rnorm(n, mean = y.mean, sd = 0.2)

the.data<- data.frame(X=x, Y=y, Obs=y.obs, Cov=X.cov)

#FIT SOME MODELS - 

#TEST FOR SPATIAL AUTOCOR
#test for spatial dependence with Morans I in ape()

#calculate matrix of distances between pairs
loc <- cbind(x,y)
dist<- as.matrix( dist(loc, diag = TRUE, upper = TRUE) )

#spatial autocorrelation of Obs
Moran.I(the.data$Obs, dist)

#Model without spatial autocorrelation
#Run a lm on the data
lm1<- lm(Obs~Cov, data = the.data)
summary(lm1)

#test for spatial dependence in the residuals
Moran.I(lm1$residuals, dist) 

#Look at variogram of residuals
#construct a geodata object

lm1.geo<- list(coords = loc, data = lm1$residuals)

#Compute classical estimator of the variogram
lm1.variog<- variog(lm1.geo, max.dist = max(dist)/2)

#Compute a robust variogram and compare
lm1.variog.r<- variog(lm1.geo, max.dist = max(dist)/2, estimator.type = "modulus")
plot(lm1.variog.r)
points(lm1.variog$u, lm1.variog$v, pch = 15)
legend("bottomright", c("Classical", "Robust"), pch = c(15, 1))


#fit an exponential variogram model to classical empirical variogram using weighted least squares (default)
lm1.variog.fit1<- variofit(lm1.variog, cov.model = "exp", ini.cov.pars = c(2.0, 0.5), max.dist = max(dist)/2)

#ML and REML estimation methods for classical variogram use likfit()
lm1.likfit.ml<- likfit(lm1.geo, cov.model = "exp", ini.cov.pars = c(2.0, 0.5))
lm1.likfit.reml<- likfit(lm1.geo, cov.model = "exp", ini.cov.pars = c(2.0, 0.5), lik.method = "REML")

#Plot it up:
plot(lm1.variog, pch = 15)
lines(lm1.variog.fit1)
lines(lm1.likfit.ml, lty = 2)
lines(lm1.likfit.reml, lty = 3)
legend("bottomright", c("WLS", "ML", "REML"), lty = c(1,2,3))


#CAN ALSO PERFORM ONE STEP ESTIMATION

#Want to estimate the trend and the spatial covariance in 1 step:
#construct a geodata dataset
the.geo<- list(coords = loc, data = the.data$Obs, Cov = the.data$Cov)

#Estimate the trend and covariance parameters together
#Model using ML
geo.ml1<- likfit(the.geo, trend = trend.spatial(~Cov), ini.cov.pars = c(2, 0.5), cov.model = "exp") 
geo.ml1$parameters.summary

geo.reml1<- likfit(the.geo, trend = trend.spatial(~Cov), ini.cov.pars = c(2, 0.5), cov.model = "exp", lik.method = "REML") 

#compare estimates and SE's
geo.reml1$parameters.summary
geo.reml1$AIC

#compare to no spatial autocorr
geo.reml1$nospatial

#estimates of the trend SE's - 
sqrt(diag(geo.reml1$beta.var) )

#compare with Estimate and SE from the linear model
summary(lm1)$coef

#######
#Note - a much larger dataset is needed to be able to disentangle mean from covariance factors (see bottom of lab)
#######

#KRIGING
#calculate covariate value across the grid
p.grid<- expand.grid(seq(-1, 1, l = 51), seq(-1,1, l = 51) )
grid.Cov<- 0.2*p.grid[,1] -0.4*p.grid[,2]

k.field<- krige.conv(the.geo, loc = p.grid, krige = krige.control(obj.model = geo.ml1, type.krige = "ok", trend.d = trend.spatial(~Cov, the.geo), trend.l =  trend.spatial(~grid.Cov ) ) ) 

k.field2<- krige.conv(the.geo, loc = p.grid, krige = krige.control(obj.model = geo.reml1, type.krige = "ok", trend.d = trend.spatial(~Cov, the.geo), trend.l =  trend.spatial(~grid.Cov ) ) ) 

#compare kriged estimates ML and REML
par(mfrow = c(1,2))
image(k.field, main = "ML Model", ylim = c(-1,1), useRaster=TRUE)
image(k.field2, main = "REML Model", ylim = c(-1,1), useRaster=TRUE)

#To Run when have some time...

if(F){
#modify number of sampling locations to evaluate ability to identify trend versus covariance components - 
# first we simulate some random sampling locations
n <- 1000
x <- runif(n=n, min=-1, max=1)
y <- runif(n=n, min=-1, max=1)
loc <- cbind(x,y)
X.cov<- 0.2*x -0.4*y
re.cov.model<- RMexp(var = 2.5, scale = 0.4)
#re.data <- RFsimulate(model = re.cov.model, x=x, y=y, grid=FALSE) 
#RE<- re.data$variable1
re.data<- grf(grid = loc, cov.pars = c(2.5, 0.4))
RE<- re.data$data
beta0<- 0.3
beta1<- 0.8
y.trend<- beta0 + beta1*X.cov
#underlying mean with spatially autocor random effects
y.mean<- y.trend + RE
#observed data
y.obs<- rnorm(n, mean = y.mean, sd = 0.2)
the.data<- data.frame(X=x, Y=y, Obs=y.obs, Cov=X.cov)
#Want to estimate the trend and the spatial covariance in 1 step:
#construct a geodata dataset
the.geo<- list(coords = loc, data = the.data$Obs, Cov = the.data$Cov)
geo.fit<- likfit(the.geo, trend = trend.spatial(~Cov), ini.cov.pars = c(2, 0.5), cov.model = "exp", nugget = ) 
geo.fit$parameters.summary
geo.fit$AIC

#estimates of the trend SE's - 
sqrt(diag(geo.fit$beta.var) )

#compare to no spatial autocorr
geo.fit$nospatial

}

