
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2016_Spatio-temporal_models/Week 2 -- mixed-effects/Lecture 2" )
Use_REML = FALSE

######################
# Simulate data
######################

# Simulate predictors
Factor = rep( 1:10, each=10)
Z = rnorm( length(unique(Factor)), mean=0, sd=1)
X0 = 0

# Simulate response
Y = Z[Factor] + X0 + rnorm( length(Factor), mean=0, sd=1)

######################
# Run in R
######################

library(lme4)
Lme = lmer( Y ~ 1|factor(Factor), REML=Use_REML)

######################
# Run in TMB
######################

library(TMB)

# Compile model
Version = "linear_mixed_model"
compile( paste0(Version,".cpp") )

# Build inputs
Data = list( "n_data"=length(Y), "n_factors"=length(unique(Factor)), "Factor"=Factor-1, "Y"=Y)
Parameters = list( "X0"=-10, "log_SD0"=2, "log_SDZ"=2, "Z"=rep(0,Data$n_factor) )
Random = c("Z")
if( Use_REML==TRUE ) Random = union( Random, "X0")

# Build object
dyn.load( dynlib("linear_mixed_model") )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)  #

# Prove that function and gradient calls work
Obj$fn( Obj$par )
Obj$gr( Obj$par )

# Optimize
start_time = Sys.time()
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, control=list("trace"=1) )
  Opt[["final_gradient"]] = Obj$gr( Opt$par )
  Opt[["total_time"]] = Sys.time() - start_time

# Get reporting and SEs
Report = Obj$report()
  SD = sdreport( Obj )

######################
# Shrinkage estimator
######################

# jim's attempt at replicating this from first principles
Mu = mean(Y)
  Mu_s = tapply( Y, INDEX=Factor, FUN=mean)
Sigma = sd( Mu_s )
  Sigma_s = sd( Y - Mu_s[Factor] )
Weights_hat = c( 1/Sigma^2, length(Y)/length(unique(Factor))/Sigma_s^2 )
  Weights_hat = Weights_hat / sum(Weights_hat)
# Predictions
Mu_s_hat = ( Mu*Weights_hat[1] + Mu_s*Weights_hat[2] )
cbind( Mu_s_hat-Mu, ranef(LMM)[['factor(Site_i)']]+fixef(LMM)['(Intercept)'] )

######################
# Compare estimates
######################

# Global mean
c( fixef(Lme), Report$X0, Mu )

# Random effects
cbind( "True"=Z, ranef(Lme)[['factor(Factor)']], Report$Z, Mu_s-Mu )

# Variances
summary(Lme)
unlist( Report[c("SDZ","SD0")] )

