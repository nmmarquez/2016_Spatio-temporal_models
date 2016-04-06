rm(list=ls())
setwd("~/Documents/Classes/2016_Spatio-temporal_models/Week 1 -- Likelihoods and linear models/Homework/")
set.seed(777)
library(TMB)
library(ggplot2)
library(SpatialDeltaGLMM)
library(knitr)

# Compile the cpp code
if (file.exists("hw1.so")) file.remove("hw1.so")
if (file.exists("hw1.o")) file.remove("hw1.o")
if (file.exists("hw1.dll")) file.remove("hw1.dll")
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
K <- 10

# simulate the covariates to use for the analysis
simulate_covs <- function(N){
    cbind(rep(1, N), sapply(1:ncov, function(i) 
        rnorm(N, mean(X[,i+1]), 1)))
}

# based on a set of covariates and model form simulate a catch
simulate_catch <- function(covs, use_gamma, betas, sigma, zero_prob){
    exp_val <- exp(covs %*% betas)
    zero_catches <-rbinom(N, 1, 1-zero_prob)
    if(use_gamma){
        shape <- exp_val
        catches <- rgamma(N, shape, scale=sigma) * zero_catches
    }
    else{
        meanlog <- log(exp_val)
        catches <- rlnorm(N, meanlog, sigma) * zero_catches
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

col_mean_sans_zero <- function(mat, log){
    if (log){
        mat_means <- apply(mat, 2, function(x) mean(log(x[x!=0])))
    }
    else{
        mat_means <- apply(mat, 2, function(x) mean(x[x!=0]))
    }
    mat_means
}

# check out some simulations
hist_sim <- function(model, i=1, log=TRUE, catch_sim = catch_simulations){
    if(log){
        datur <- catch_sim[[model]][,i]
    }
    else{
        datur <- log(catch_sim[[model]][,i]+1)
    }
    hist(datur, main=paste0("Simulation with ", model, " version ", i), 
         xlab="simulated observation values")
}


hist_sim("log_normal", 10)
hist_sim("gamma", 10)
hist_sim("log_normal_cov", 10)
hist(EBS_pollock_data$catch)

hist(col_mean_sans_zero(catch_simulations$log_normal, log=T))
hist(col_mean_sans_zero(catch_simulations$gamma, log=T))
hist(col_mean_sans_zero(catch_simulations$log_normal_cov, log=T))
hist(col_mean_sans_zero(catch_simulations$log_normal, log=F))
hist(col_mean_sans_zero(catch_simulations$gamma, log=F))
hist(col_mean_sans_zero(catch_simulations$log_normal_cov, log=F))

# no more k fold cross validation
K <- 1

# run each model on each simulation
model_sims <- lapply(catch_simulations, function(model_response) 
    lapply(1:M, function(i)
        lapply(model_spec, function(x)
            run_model(x$use_gamma, x$use_cov, covariate_simulations[[i]], 
                      model_response[,i]))))

# intercept estimate bias 
int_est <- lapply(model_sims, function(x) t(sapply(1:length(x), function(i)
    sapply(x[[i]], function(j) j$Report$beta[1]))))


par(mfrow=c(3,3), mar=c(2, 4, 4, 2) + 0.1)
for(i in 1:3){
    for(j in 1:3){
        main_ <- ifelse(i==1, paste0('Simulated with ', names(models)[j]),'')
        ylab_ <- ifelse(j==1, paste0('Modeled as ', names(models)[i]),'')
        font.lab_ <- ifelse(j==1, 2, 1)
        int_ <- int_est[[j]][,i]
        true_ <- models[[j]]$Report$beta[1]
        sd_ <- sd(int_)
        xlim_ <- c(min(c(min(int_, true_)))-sd_, max(c(max(int_, true_)))+sd_)
        plot(density(int_), xlim=xlim_, ylab=ylab_, main=main_,
             font.lab=font.lab_, xlab='', cex.lab=1.2)
        abline(v=models[[j]]$Report$beta[1], col="red", lwd=2)
    }
}
par(mfrow=c(1,1))

# plot difference in predictive OOS jnll
pred_jnll <- lapply(model_sims, function(x) t(sapply(1:length(x), function(i)
    sapply(x[[i]], function(j) mean(j$PredNLL_k / table(Partition_i))))))

par(mfrow=c(3,3), mar=c(2, 4, 4, 2) + 0.1)
for(i in 1:3){
    for(j in 1:3){
        main_ <- ifelse(i==1, paste0('Simulated with ', names(models)[j]),'')
        ylab_ <- ifelse(j==1, paste0('Modeled as ', names(models)[i]),'')
        font.lab_ <- ifelse(j==1, 2, 1)
        jnll_diff <- pred_jnll[[j]][,i] - pred_jnll[[j]][,j]
        plot(pred_jnll[[j]][,j], pred_jnll[[j]][,i], ylab=ylab_, main=main_,
             font.lab=font.lab_, xlab='', cex.lab=1.2)
        abline(a=0, b=1, col="red")
    }
}
par(mfrow=c(1,1))

jnll_df <- as.data.frame(lapply(models, function(x) sum(x$Report$jnll_vec)))
rownames(jnll_df) <- "jnll"
kable(jnll_df, format="markdown")

param_mat <- t(sapply(model_params, function(x) x))
param_mat[param_mat == 0] <- NA
param_df <- as.data.frame(param_mat)
names(param_df) <- c("theta", "sigma", "beta_int", "beta_lat", "beta_long")
kable(param_df, format="markdown")

log_pred_score_df <- as.data.frame(log_pred_score)
rownames(log_pred_score_df) <- "pred_score"
kable(log_pred_score_df, format="markdown")
