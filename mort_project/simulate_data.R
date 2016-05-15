rm(list=ls())
library(TMB)
setwd("~/Documents/Classes/2016_Spatio-temporal_models/mort_project/")
source("./mort_viz/utilities.R")

sim_data <- function(x, a=.1, b=.1, c=.1, d=.1, f=.1){
    a * exp(-1 * b * x) + c + d * exp(x * f)
}

plot(log(sim_data(1:100)))

model <- "gpz"
df <- download_mort_data(102)
time_plot(df, 1)

reload_model(model)

df$age_group <- df$age_group_id - min(df$age_group_id)

run_model <- function(df, option=1, model_name=model){
    N <- nrow(df)
    param_num <- list("0"=3, "1"=5, "2"=8)[[as.character(option)]]
    Map <- list()
    dyn.load(dynlib(model_name))
    Params <- list(log_gpz=rep(0, param_num), log_sigma_obs=0, 
                   age_fixed=rep(0, length(unique(df$age_group))))
    Data <- list(log_rate_mort=df$log_rate, age=df$age_mean, option=option,
                 age_group=df$age_group)
    Obj <- MakeADFun(data=Data, parameters=Params, DLL=model_name, silent=T)
    Obj$env$tracemgc <- FALSE
    Obj$env$inner.control$trace <- FALSE
    Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr,
                  control=list(iter.max=1000, eval.max=500))
    Report <- Obj$report()
    dyn.unload(dynlib(model_name))
    if(Opt$convergence != 0){
        print("Model had partial or non_convergence!")
    }
    Report
}

gpzm_sub <- run_model(subset(df, age_mean > 25), option=0)
gpzm <- run_model(df, option=0)
silder <- run_model(subset(df, age_mean < 80), option=1)
silder <- run_model(df, option=1)

test_ages <- sort(unique(df$age_mean))
fitted_vals <- sim_data(test_ages, a=silder$alpha, b=silder$beta,
                        c=silder$gamma, d=silder$delta, f=silder$zeta)
fitted_df <- data.frame(test_ages, log(fitted_vals), year_id=2020)
time_plot(df, 1) + geom_line(aes(test_ages, log.fitted_vals.), fitted_df)
fitted_df$log.fitted_vals. <- fitted_df$log.fitted_vals. + silder$age_fixed 
time_plot(df, 1) + geom_line(aes(test_ages, log.fitted_vals.), fitted_df)
