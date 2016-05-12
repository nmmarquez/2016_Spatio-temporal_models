rm(list=ls())
set.seed(123)
library(TMB)
library(knitr)
setwd("~/Documents/Classes/2016_Spatio-temporal_models/Week 3 -- Temporal Models/Homework/")

model_name <- "hw3"
if (file.exists(paste0(model_name, ".so"))) file.remove(paste0(model_name, ".so"))
if (file.exists(paste0(model_name, ".o"))) file.remove(paste0(model_name, ".o"))
if (file.exists(paste0(model_name, ".dll"))) file.remove(paste0(model_name, ".dll"))
compile(paste0(model_name, ".cpp"))

M <- 100

DLM_sim <- function(N=100, sigma.obs=.4, sigma.proc=0.4, a=0, y1=10, b=-.5){
    ytrue <- numeric(length = N)
    ytrue[1] <- y1
    proc.error <- rnorm(N, mean = 0, sd = sigma.proc)
    for(i in 2:N) {
        ytrue[i] <- a + b*ytrue[i-1] + proc.error[i-1]
    }
    x <- seq_len(N)
    y <- rnorm(N, mean = ytrue, sd = sigma.obs)
    data.frame(y = y, ytrue = ytrue)
}

run_model <- function(data, use_alpha=FALSE){
    N <- nrow(data)
    random <- c("u")
    dyn.load(dynlib(model_name))
    Params <- list(a=rep(0, use_alpha), b=0, log_sigma_proc=0,
                   log_sigma_obs=0, u=rep(mean(data$y), N))
    Data <- list(y=data$y)
    Obj <- MakeADFun(data=Data, parameters=Params, DLL=model_name, random=random)
    Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr)
    Report <- Obj$report()
    SD <- sdreport(Obj)
    dyn.unload(dynlib(model_name))
    params <- c("b", "sigma_proc", "sigma_obs")
    if (use_alpha){
        params <- c("a", "b", "sigma_proc", "sigma_obs")
    }
    summary(SD)[params, 1]
}