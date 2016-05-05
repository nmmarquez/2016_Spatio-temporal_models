
setwd("~/Documents/Classes/2016_Spatio-temporal_models/Week 6 -- 2D spatial models/Lecture/")

# devtools::install_github("james-thorson/TMB-debug")

library(TMB)
library(RandomFields)

###################
# Equal distance 2D autoregressive
###################

Dim = c("n_x"=10, "n_y"=10)
loc_xy = expand.grid("x"=1:Dim['n_x'], "y"=1:Dim['n_y'])
Scale = 2
Sigma2 = (0.5) ^ 2
beta0 = 3
prob_missing = 0.2

# Simulate spatial process
RMmodel = RMgauss(var=Sigma2, scale=Scale)
epsilon_xy = array(RFsimulate(model=RMmodel, x=loc_xy[,'x'], y=loc_xy[,'y'])@data[,1], dim=Dim)
image( z=epsilon_xy )

# SImulate counts
c_xy = array(NA, dim=dim(epsilon_xy))
for(x in 1:nrow(c_xy)){
for(y in 1:ncol(c_xy)){
  c_xy[x,y] = rpois(1, exp(beta0 + epsilon_xy[x,y]) )
  if( rbinom(n=1, size=1, prob=prob_missing)==1) c_xy[x,y] = NA
}}
true_abundance =  sum( exp(beta0 + epsilon_xy) )
image(c_xy)
# Compile
Params = list( "beta0"=0, "ln_sigma2"=0, "logit_rho"=0, "epsilon_xy"=array(rnorm(prod(dim(loc_xy))),dim=dim(epsilon_xy)) )
compile( "autoregressive_grid_V1.cpp" )
dyn.load( dynlib("autoregressive_grid_V1") )

######## Version 0 -- Sweep downstream
# Build object
Data = list("Options_vec"=c(0), "c_xy"=c_xy )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_xy", DLL="autoregressive_grid_V1" )
# Optimize
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
par0 = Opt$par
h0 = Obj$env$spHess(random=TRUE)
report0 = Obj$report()
sd0 = sdreport( Obj, bias.correct=TRUE )

######## Version 1 -- Analytic precision matrix
# Build object
Data = list("Options_vec"=c(1), "c_xy"=c_xy )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_xy", DLL="autoregressive_grid_V1" )
# Optimize
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
par1 = Opt$par
h1 = Obj$env$spHess(random=TRUE)
report1 = Obj$report()
sd1 = sdreport( Obj, bias.correct=TRUE )

######## Version 3 -- Built-in function for AR process
# Build object
Data = list("Options_vec"=c(3), "c_xy"=c_xy )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_xy", DLL="autoregressive_grid_V1" )
# Optimize
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
par3 = Opt$par
h3 = Obj$env$spHess(random=TRUE)
report3 = Obj$report()
sd3 = sdreport( Obj, bias.correct=TRUE )

