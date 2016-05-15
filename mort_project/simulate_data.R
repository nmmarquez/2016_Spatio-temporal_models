rm(list=ls())
library(TMB)
setwd("~/Documents/Classes/2016_Spatio-temporal_models/mort_project/")
source('./mort_viz/db_access.R')
source("./mort_viz/utilities.R")

sim_data <- function(x, a=.1, b=.1, c=.1, d=.1, f=.1){
    a * exp(-1 * b * x) + c + d * exp(x * f)
}

if (!("usa_data.rda" %in% list.files())){
    df <- download_mort_data(102)
    save(df, file="./usa_data.rda")
}

load("./usa_data.rda")

model <- "gpz"

time_plot(df, 1)

reload_model(model)

df$age_group <- df$age_group_id - min(df$age_group_id)
df <- subset(df, sex_id == 1)

run_model <- function(df, option=1, model_name=model, print=F){
    N <- nrow(df)
    param_num <- list("0"=3, "1"=5, "2"=8)[[as.character(option)]]
    Map <- list()
    Random <- c("epsilon_age")
    if (option != 1){
        Map[["epsilon_age"]] <- factor(rep(NA, length(unique(df$age_group))))
        Map[["log_sigma_age"]] <- factor(NA)
        Random <- NULL
    }
    dyn.load(dynlib(model_name))
    Params <- list(log_gpz=rep(0, param_num), log_sigma_obs=0, log_sigma_age=0,
                   epsilon_age=rep(0, length(unique(df$age_group))),
                   logit_rho=0)
    Data <- list(log_rate_mort=df$log_rate, age=df$age_mean, option=option,
                 age_group=df$age_group, unique_age=unique(df$age_mean))
    Obj <- MakeADFun(data=Data, parameters=Params, DLL=model_name,
                     silent=!print, map=Map, random = Random)
    Obj$env$tracemgc <- print
    Obj$env$inner.control$trace <- print
    Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr,)
                  #control=list(iter.max=1000, eval.max=500))
    Report <- Obj$report()
    dyn.unload(dynlib(model_name))
    if(Opt$convergence != 0){
        print("Model had partial or non_convergence!")
    }
    Report
}

#gpzm_sub <- run_model(subset(df, age_mean > 25), option=0)
#gpzm <- run_model(df, option=0)
#silder <- run_model(subset(df, age_mean < 80), option=1)
silder <- run_model(df, option=1, print=T)

test_ages <- sort(unique(df$age_mean))
fitted_vals <- sim_data(test_ages, a=silder$alpha, b=silder$beta,
                        c=silder$gamma, d=silder$delta, f=silder$zeta)
fitted_df <- data.frame(test_ages, log(fitted_vals), year_id=2020)
time_plot(df, 1) + geom_line(aes(test_ages, log.fitted_vals.), fitted_df)
fitted_df$log.fitted_vals. <- fitted_df$log.fitted_vals. + silder$age_fixed
time_plot(df, 1) + geom_line(aes(test_ages, log.fitted_vals.), fitted_df)