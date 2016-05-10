rm(list=ls())
set.seed(123)
library(TMB)
setwd( "~/Documents/Classes/2016_Spatio-temporal_models/Week 5 -- 1D spatial models/Lecture/" )

###################
# Equal distance autoregressive
###################

x = 1:100
Rho = 0.8
Sigma2 = (1) ^ 2
n_rep = 3
beta0 = 3

# Simulate spatial process
epsilon_s = rep(NA, length(x))
epsilon_s[1] = rnorm(1, mean=0, sd=sqrt(Sigma2))
for(s in 2:length(x)) epsilon_s[s] = Rho*epsilon_s[s-1] + rnorm(1, mean=0, sd=sqrt(Sigma2))

# SImulate counts
c_si = matrix( nrow=length(x), ncol=n_rep)
for(s in 1:nrow(c_si)){
for(i in 1:ncol(c_si)){
  c_si[s,i] = rpois(1, exp(beta0 + epsilon_s[s]) )
}}

# Compile
Params = list( "beta0"=0, "ln_sigma2"=0, "logit_rho"=0, "epsilon_s"=rnorm(length(x)) )
model <- "autoregressive_V1_edit"
if (file.exists(paste0(model, ".so"))) file.remove(paste0(model, ".so"))
if (file.exists(paste0(model, ".o"))) file.remove(paste0(model, ".o"))
if (file.exists(paste0(model, ".dll"))) file.remove(paste0(model, ".dll"))
compile(paste0(model, ".cpp"))
dyn.load( dynlib(model) )

######## Version 0 -- Stochastic process with automatic sparseness detection
# Build object
Data = list("Options_vec"=c(0), "c_si"=c_si )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_s", DLL=model )
# Optimize
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
par0 = Opt$par
h0 = Obj$env$spHess( random=TRUE )
h0f = Obj$env$spHess()
Report0 <- Obj$report()

######## Version 1 -- Analytic precision matrix
# Build object
Data = list("Options_vec"=c(1), "c_si"=c_si )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_s", DLL=model )
# Optimize
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
par1 = Opt$par
h1 = Obj$env$spHess(random=TRUE)
h1f = Obj$env$spHess()
Report1 <- Obj$report()

######## Version 2 -- Covariance and built-in function
# Build object
Data = list("Options_vec"=c(2), "c_si"=c_si )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_s", DLL=model)
# Optimize
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
par2 = Opt$par
h2 = Obj$env$spHess(random=TRUE)
Report2 <- Obj$report()

######## Version 3 -- Built-in function for AR process
# Build object
Data = list("Options_vec"=c(3), "c_si"=c_si )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_s", DLL=model )
# Optimize
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
par3 = Opt$par
h3 = Obj$env$spHess(random=TRUE)
Report3 <- Obj$report()

######## Version 4 -- ar not scaled
# Build object
Map <- list(ln_sigma2=factor(NA))
Data = list("Options_vec"=c(4), "c_si"=c_si )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_s", DLL=model , map=Map)
# Optimize
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
par3 = Opt$par
h3 = Obj$env$spHess(random=TRUE)
Report4 <- Obj$report()

######## Version 4 -- ar not scaled
# Build object
Map <- list(ln_sigma2=factor(NA))
Data = list("Options_vec"=c(5), "c_si"=c_si )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_s", DLL=model , map=Map)
# Optimize
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
par3 = Opt$par
h3 = Obj$env$spHess(random=TRUE)
Report5 <- Obj$report()


print("method 0 ")
print(Report0$rho)
print("method 1 ")
print(Report1$rho)
print("method 2 ")
print(Report2$rho)
print("method 3 ")
print(Report3$rho)
print("method 4 ")
print(Report4$rho)
print("method 5 ")
print(Report5$rho)
