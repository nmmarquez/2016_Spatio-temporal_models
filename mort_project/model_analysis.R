rm(list=ls())
library(TMB)
library(INLA)
library(MASS)

setwd("~/Documents/Classes/2016_Spatio-temporal_models/mort_project/")
source('./mort_viz/db_access.R')
source("./mort_viz/utilities.R")

model <- "gpz_log"
reload_model(model)

if (!("data.rda" %in% list.files())){
    df <- download_mort_data()
    save(df, file="./data.rda")
}

load("./data.rda")

# load in data and only look at 1 sex ages under 100
df <- subset(df, sex_id == 1 & age_mean<=100 & year_id>=1990)
df$age_group <- df$age_group_id - min(df$age_group_id)
df$time_group <- df$year_id - min(df$year_id)
df$loc_group <- as.numeric(as.factor(df$location_id)) - 1
option <- 0
print <- TRUE
model_name <- model
time_plot(df, 1)

# run model
run_model <- function(df, option=1, model_name=model, print=F){
    N_ <- nrow(df)
    T_ <- length(unique(df$time_group))
    A_ <- length(unique(df$age_group))
    L_ <- length(unique(df$loc_group))
    Map <- list()
    Map[["phi"]] <- factor(array(NA, dim=c(L_, A_, T_)))
    dim(Map[["phi"]]) <- c(L_, A_, T_)
    Map[["epsilon"]] <- factor(array(NA, dim=c(L_)))
    dim(Map[["epsilon"]]) <- c(L_)
    Map[["logit_rho_loc"]] <- factor(NA)
    Map[["logit_rho_age"]] <- factor(NA)
    Map[["logit_rho_time"]] <- factor(NA)
    Map[["log_sigma_loc"]] <- factor(NA)
    Map[["log_sigma_age"]] <- factor(NA)
    Map[["log_sigma_time"]] <- factor(NA)
    Random <- NULL
    if(option == 1){
        Random <- "phi"
        Map[["logit_rho_age"]] <- NULL
        Map[["log_sigma_age"]] <- NULL
        Map[["logit_rho_time"]] <- NULL
        Map[["log_sigma_time"]] <- NULL
        Map[["phi"]] <- NULL
    }
    if (option == 2){
        Random <- c("epsilon")
        #Map[["logit_rho_loc"]] <- NULL
        #Map[["log_sigma_loc"]] <- NULL
        Map[["epsilon"]] <- NULL
    }
    dyn.load(dynlib(model_name))
    par_num <- ifelse(model == "gpz", 5, 6)
    Params <- list(gpz=rep(0, par_num), log_sigma_obs=0, logit_rho_time=0,
                   logit_rho_age=0, logit_rho_loc=log(.9), log_sigma_time=0,
                   log_sigma_age=0, log_sigma_loc=0, 
                   phi=array(0, dim=c(L_, A_, T_)), epsilon=array(0, dim=c(L_)),
                   secular=0.)
    Data <- list(log_rate_mort=df$log_rate, age=df$age_mean, option=option,
                 age_group=df$age_group, time_group=df$time_group,
                 graph=admin_queens(prec=F), loc_group=df$loc_group)
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

silder <- run_model(df, option=0, print=T)

# get the global values to use for simulation
df$log_rate_hat <- silder$log_rate_mort_hat
time_plot(df, 1, preds=T)


add_df <- subset(df, year_id > 2000)
add_df$year_id <- add_df$year_id + 15
add_df$log_rate <- NA
add_df$rate <- NA
add_df$envelope <- NA
add_df$population <- NA

df2 <- rbind(df, add_df)
df2$time_group <- df2$year_id - min(df2$year_id)

silder3 <- run_model(df2, option=1, print=T)
df2$log_rate_hat <- silder3$log_rate_mort_hat
time_plot(subset(df2, loc_group==8), 1, preds=T)
age_plot(subset(df2, loc_group==8), 1, preds=T)
silder$Q_loc
silder3$rho_age
silder3$rho_time
silder3$sigma_time
silder3$sigma_age
silder3$secular
