
#lecture week 4
#Spatio temporal class
#  install.packages("geoR", contriburl="http://leg.ufpr.br/~paulojus/geoR")
library(geoR)




#covariance functions 
curve(cov.spatial(x, cov.pars=c(1, .2)), from = 0, to = 1,
      xlab = "distance", ylab = "C(h)",
      main = "covariance functions", col = 2, lwd = 3)
curve(cov.spatial(x, cov.pars = c(1, .6), cov.model = "sph"), 0, 1,
      add = TRUE, lty = 1, col = 3, lwd = 3)
curve(cov.spatial(x, cov.pars = c(1, .6/sqrt(3)), cov.model = "gau"),
      0, 1, add = TRUE,  col = 4, lwd = 3)
legend("topright", c("exponential", "spherical", "gaussian"),
       lty=c(1,1,1), lwd=c(3,3,3), col = c(2,3,4))


# variograms - from the examples on cov.spatial() in geoR
v.f <- function(x, ...){1-cov.spatial(x, ...)}
curve(v.f(x, cov.pars=c(1, .2)), from = 0, to = 1,
      xlab = "distance", ylab = expression(gamma(h)),
      main = "variograms with equivalent \"practical range\"", col = 2, lwd = 3)
curve(v.f(x, cov.pars = c(1, .6), cov.model = "sph"), 0, 1,
      add = TRUE, lty = 2, col = 3, lwd = 3)
curve(v.f(x, cov.pars = c(1, .6/sqrt(3)), cov.model = "gau"),
      0, 1, add = TRUE,  col = 4, lwd = 3)
legend("topleft", c("exponential", "spherical", "gaussian"),
        lty=c(1,1,1), lwd=c(3,3,3), col = c(2,3,4))
        
lines.variomodel(cov.model = "exp", cov.pars = c(1, 0.2), nugget = 0.15, max.dist = 1, lwd = 3, ylim  = c(0, 1.2) , ylab = expression(gamma(h)), xlab = "h")
 
 
 
        
#simulate some data
#uses exponential model as default with cov.pars = sigma2 & phi
#the underlying covariance function
curve(cov.spatial(x, cov.pars=c(0.5, .2)), from = 0, to = 1,
      xlab = "distance", ylab = "C(h)",
      main = "covariance function", col = 2, lwd = 3)
      
sim <- grf(grid = expand.grid(x = seq(0.0555, 0.944444, l = 8), y = seq(0.0555, 0.944444, l = 8)), cov.pars = c(0.5, 0.2))

#View Gaussian Random Field with regular grid of data:
image(sim, col = terrain.colors(20) )

    
#simulate some data on an irregular set of points
loc<- matrix(nrow = 100, ncol = 2, data = runif(200) )
#uses exponential model as default with cov.pars = sigma2 & phi
sim <- grf(grid = loc, cov.pars = c(0.5, 0.2)) 

#plot locations and levels with colors
points(sim, pt.divide = "quintile", xlab = "Coord X", ylab = "Coord Y")     

#calculate variogram clouds - both methods

cloud1 <- variog(sim, option = "cloud", max.dist = 1)
cloud2 <- variog(sim, option = "cloud", estimator.type = "modulus", max.dist = 1)
par(mfrow = c(1,2))
plot(cloud1, main = "classical estimator")
plot(cloud2, main = "modulus estimator")

#calculate binned variogram clouds
bin1 <- variog(sim, uvec = seq(0, 0.8, l = 8), bin.cloud = T)
bin2 <- variog(sim, uvec = seq(0, 0.8, l = 8), estimator.type = "modulus", bin.cloud = T)

par(mfrow = c(1,2))
plot(bin1, bin.cloud = T, main = "classical estimator")
plot(bin2, bin.cloud = T, main = "modulus estimator")
        
#calculate empirical variograms
vario1 <- variog(sim, uvec = seq(0, 0.8, l = 8))
vario2<- variog( sim, uvec = seq(0, 0.8, l = 8), estimator.type = "modulus")
par(mfrow = c(1,2))
plot(vario1, main = "classical estimator", pch = 15)
plot(vario2, main = "modulus estimator", pch = 15)


#evaluate variogram fits with s100 data set

#fit by eye by plotting a specific model#lines.variomodel(cov.model = "exp", cov.pars = c(0.8,  0.3), nug = 0, max.dist = 0.8, lty = 2)

#fit with various methods, nugget fixed to 0.15, exponential variogram, classic variogram estimator
#Note: beta is an estimate of the mean or trend if made a function of coordinates or covariates
ml <- likfit(s100, ini = c(1, 0.5), fix.nugget = T, nugget = 0.15)reml <- likfit(s100, ini = c(1, 0.5), fix.nugget = T, lik.method = "RML", nugget = 0.15)

#calculate classical empirical variogram for s100
bin<- variog(s100, uvec = seq(0, 1, l = 11), bin.cloud = T)

ols <- variofit(bin, ini = c(1, 0.5), fix.nugget = T, nugget = 0.15, weights = "equal")wls <- variofit(bin, ini = c(1, 0.5), fix.nugget = T, nugget = 0.15)

plot(bin, ylim = c(0, 1.2))
lines(ml, max.dist = 1, lwd = 2)
lines(reml, lwd = 2, max.dist = 1, col = 2)lines(ols, lwd = 2, max.dist = 1, col = 3)lines(wls,  lwd = 2, max.dist = 1, col = 4)legend("bottomright", legend = c("ML", "REML", "OLS","WLS"), lwd = c(2, 2, 2, 2), col = c(1, 2, 3, 4), cex = 0.7)




#look at model
ml1 





       
        
        
               