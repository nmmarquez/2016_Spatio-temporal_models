rm(list=ls())
library(TMB)
library(INLA)
setwd("~/Documents/Classes/2016_Spatio-temporal_models/mort_project/")
source('./mort_viz/db_access.R')
source("./mort_viz/utilities.R")

reload_model(model)

if (!("usa_data.rda" %in% list.files())){
    df <- download_mort_data(102)
    save(df, file="./usa_data.rda")
}

load("./usa_data.rda")

model <- "gpz_log"

df$age_group <- df$age_group_id - min(df$age_group_id)
df$time_group <- df$year_id - min(df$year_id)
df <- subset(df, sex_id == 1 & age_mean <=15)
option <- 0
print <- TRUE
model_name <- model
time_plot(df, 1)

run_model <- function(df, option=1, model_name=model, print=F){
    N_ <- nrow(df)
    T_ <- length(unique(df$time_group))
    A_ <- length(unique(df$age_group))
    Map <- list()
    Random <- c("epsilon_age", "epsilon_time", "epsilon_age_time")
    if (option <= 0){
        Map[["epsilon_age"]] <- factor(rep(NA, length(unique(df$age_group))))
        Random <- NULL
    }
    if (option <= 1){
        Map[["logit_rho_age"]] <- factor(NA)
        Random <- NULL
    }
    if (option <= 2){
        Map[["epsilon_time"]] <- factor(rep(NA, length(unique(df$time_group))))
        Map[["logit_rho_time"]] <- factor(NA)
        Random <- c("epsilon_age")
    }
    if (option <= 3){
        Map[["epsilon_age_time"]] <- factor(matrix(NA, nrow=A_, ncol=T_))
        dim(Map[["epsilon_age_time"]]) <- c(A_, T_)
        Map[["logit_rho_age2"]] <- factor(NA)
        Map[["logit_rho_time2"]] <- factor(NA)
        Random <- c("epsilon_age", "epsilon_time")
    }
    Random <- NULL
    dyn.load(dynlib(model_name))
    par_num <- ifelse(model == "gpz", 5, 3)
    Params <- list(gpz=rep(0, par_num), log_sigma_obs=0, logit_rho_time=0,
                   epsilon_age=rep(0, length(unique(df$age_group))),
                   epsilon_time=rep(0, length(unique(df$time_group))),
                   logit_rho_age=0, logit_rho_age2=0, logit_rho_time2=0,
                   epsilon_age_time=matrix(0, nrow=A_, T_))
    Data <- list(log_rate_mort=df$log_rate, age=df$age_mean, option=option,
                 age_group=df$age_group, time_group=df$time_group)
    Obj <- MakeADFun(data=Data, parameters=Params, DLL=model_name,
                     silent=!print, map=Map, random = Random)
    Obj$env$tracemgc <- print
    Obj$env$inner.control$trace <- print
    Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr,
                  control=list(iter.max=1000, eval.max=500))
    Report <- Obj$report()
    dyn.unload(dynlib(model_name))
    if(Opt$convergence != 0){
        print("Model had partial or non_convergence!")
    }
    Report
}

silder <- run_model(df, option=0, print=T)


df$log_rate_hat <- inf_term(df$age_mean, N0=silder$N0, 
                            lambda=silder$lambda, c=silder$c)
time_plot(df, 1, preds=T)

df$log_rate_hat <- df$log_rate_hat + ya_term(df$age_mean, eta=silder$eta, 
                                             scale=silder$scale,
                                             ceiling=silder$ceiling)
time_plot(df, 1, preds=T)

df$log_rate_hat <- silder$log_rate_mort_hat
time_plot(df, 1, preds=T)
