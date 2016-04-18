rm(list=ls())
set.seed(123)
library(TMB)
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

model_spec <- expand.grid(b=c(-.5, 0, .5), sigma.obs=c(.2, .4, .8))

dlm_simulations <- lapply(1:nrow(model_spec), function(i) 
    lapply(1:M, function(j) 
        DLM_sim(b=model_spec$b[i], sigma.obs=model_spec$sigma.obs[i])))

dlm_analysis <- lapply(dlm_simulations, function(x) 
    t(sapply(x, function(y) run_model(y))))

plot_param_hist <- function(values, ylab, xlab, main, true_val=NULL){
    lower_ <- quantile(values, probs=.025)
    upper_ <- quantile(values, probs=.975)
    mean_ <- mean(values)
    hist(values, main=main, ylab=ylab, xlab=xlab)
    if(!is.null(true_val)){
        abline(v=true_val, col="red", lwd=2)
    }
    abline(v=mean_, col="blue", lwd=2)
    abline(v=lower_, col="blue", lty=2, lwd=2)
    abline(v=upper_, col="blue", lty=2, lwd=2)
    legend("topright", legend=c("mean", "quant", "true"), lty = c(1,2,1),
           col=c("blue", "blue", "red"), cex=.9)
}

plot_param_distribution <- function(param="b", lines=F){
    par(mfrow=c(3,3))
    if(param=="b"){
        true_vals <- rep(c(-.5, 0, .5), 3)
    }
    else if(param=="sigma_obs"){
        true_vals <- rep(c(.2, .4, .8), each=3)
    }
    else{
        true_vals <- rep(.4, 9)
    }
    for(i in 1:nrow(model_spec)){
        values_ <- dlm_analysis[[i]][,param]
        ylab_ <- paste0("sigma.obs sim = ", rep(c(.2, .4, .8), each=3))[i]
        main_ <- paste0("b sim = ", rep(c(-.5, 0, .5), 3))[i]
        if (!(i %in% c(1,4,7))){
            ylab_ <- ""
        }
        if (!(i %in% c(1,2,3))){
            main_ <- ""
        }
        plot_param_hist(values_, main=main_, ylab=ylab_, xlab="", 
                        true_val=true_vals[i])
    }
    par(mfrow=c(1,1))
}

plot_param_distribution("b")
plot_param_distribution("sigma_proc")
plot_param_distribution("sigma_obs")

Gompertz_simulations <- lapply(1:M, function(x) 
    DLM_sim(a=1.2, b=.7, sigma.obs=.5, sigma.proc=.2, y1=4))

gomp_sims <- t(sapply(Gompertz_simulations, function(x) 
    run_model(x, use_alpha=TRUE)))

par(mfrow=c(2,2))
plot_param_hist(gomp_sims[,"a"], ylab="", xlab="", main="alpha", true_val=1.2)
plot_param_hist(gomp_sims[,"b"], ylab="", xlab="", main="beta", true_val=.7)
plot_param_hist(gomp_sims[,"sigma_proc"], ylab="", xlab="", main="sigma_proc", true_val=.2)
plot_param_hist(gomp_sims[,"sigma_obs"], ylab="", xlab="", main="sigma_obs", true_val=.5)
par(mfrow=c(1,1))
