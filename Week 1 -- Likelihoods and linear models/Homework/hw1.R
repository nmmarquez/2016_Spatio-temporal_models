rm(list=ls())
setwd("~/Documents/Classes/2016_Spatio-temporal_models/Week 1 -- Likelihoods and linear models/Homework/")
set.seed(123)
library(TMB)
library(ggplot2)
library(SpatialDeltaGLMM)

# Compile the cpp code
compile("hw1.cpp")

# Clean the data for modeling
data(EBS_pollock_data)
EBS_pollock_data$intercept <- 1
y <- EBS_pollock_data$catch
X <- as.matrix(EBS_pollock_data[,c("intercept", "lat", "long")])
ncov <- ncol(X) - 1
K <- 10
Partition_i <- sample(x=1:K, size=length(y), replace=TRUE)

# function for running with different distributions and kfold validation
run_model <- function(use_gamma, use_cov, X, y){
    end_col <- ifelse(use_cov, ncol(X), 1)
    PredNLL_k = rep(NA, K)
    for(k in 1:(K+1)){
        if(k == (K+1)){
            pred_bool <- rep(0, length(y))
        }
        else{
            pred_bool <- ifelse(Partition_i==k,1,0)
        }
        dyn.load(dynlib("hw1"))
        Params <- list(beta=rep(0,ncol(as.matrix(X[,1:end_col]))), theta=0, log_sigma=0)
        Data <- list(y=y, X=as.matrix(X[,1:end_col]), pred_bool=pred_bool, 
                     use_gamma=as.numeric(use_gamma))
        Obj <- MakeADFun( data=Data, parameters=Params, DLL="hw1")
        Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr)
        Report <- Obj$report()
        Report$use_gamma <- use_gamma
        if(k != (K+1)){
            PredNLL_k[k] <- Report$pred_jnll
        }
        dyn.unload(dynlib("hw1"))
    }
    return(list(PredNLL_k=PredNLL_k, Report=Report))
}

# model specification list
model_spec <- list(log_normal=list(use_gamma=FALSE, use_cov=FALSE),
                   gamma=list(use_gamma=TRUE, use_cov=FALSE),
                   log_normal_cov=list(use_gamma=FALSE, use_cov=TRUE))

# run the models
models <- lapply(model_spec, function(x) run_model(x$use_gamma, x$use_cov, X, y))

beta_normalize <- function(b){
    ifelse(rep(length(b)==1, ncol(X)), c(b, rep(0, ncov)),b)
}

# get the info from the models
model_nll <- lapply(models, function(x) sum(x$Report$jnll_vec))
model_params <- lapply(models, function(x) c(x$Report$zero_prob, 
                                             x$Report$sigma,
                                             beta_normalize(x$Report$beta)))

log_pred_score <- lapply(models, function(x)
    mean(x$PredNLL_k / table(Partition_i)))

# simulate the data 
N <- nrow(EBS_pollock_data) # number of data points
M <- 100 # number of simulations

# simulate the covariates to use for the analysis
simulate_covs <- function(N){
    cbind(rep(1, N), sapply(1:ncov, function(i) 
        rnorm(N, mean(X[,i+1]), sd(X[,i+1]))))
}

# based on a set of covariates and model form simulate a catch
simulate_catch <- function(covs, use_gamma, beta, sigma, zero_prob){
    lin_pred <- covs %*% beta
    zero_catches <-rbinom(N, 1, 1-zero_prob)
    if(use_gamma){
        catches <- rgamma(N, exp(lin_pred), sigma) * zero_catches
    }
    else{
        catches <- rlnorm(N, lin_pred, sigma) * zero_catches
    }
    return(catches)
}

# M covariate simulations to use with the 3 models
covariate_simulations <- lapply(1:M, function(x) simulate_covs(N))

# run all the catch simulations
catch_simulations <- lapply(models, function(x) sapply(1:M, function(i)
    simulate_catch(covariate_simulations[[i]], x$Report$use_gamma, 
                   beta_normalize(x$Report$beta), x$Report$sigma, 
                   x$Report$zero_prob)))

# no more k fold cross validation
K <- 0

# run each model on each simulation
models <- lapply(catch_simulations, function(model_response) 
    lapply(1:M, function(i)
        lapply(model_spec, function(x)
            run_model(x$use_gamma, x$use_cov, covariate_simulations[[i]], 
                      model_response[,i]))))
