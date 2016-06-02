rm(list=ls())
library(TMB)
library(INLA)
library(MASS)
library(rje)
library(sparseMVN)

setwd("~/Documents/Classes/2016_Spatio-temporal_models/mort_project/")
source('./mort_viz/db_access.R')
source("./mort_viz/utilities.R")

model <- "gpz_log"
reload_model(model)
logit_rho_loc <- logit(1 / max(rowSums(as.matrix(admin_queens(prec=F))))-.0001)

if (!("data.rda" %in% list.files())){
    df <- download_mort_data()
    save(df, file="./data.rda")
}

load("./data.rda")

# run model
run_model <- function(df, option=1, model_name=model, print=F, Q=F){
    N_ <- nrow(df)
    year_dat <- (df$time_group - mean(df$time_group)) / sd(df$time_group)
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
    if(option >= 1){
        Random <- "phi"
        Map[["logit_rho_age"]] <- NULL
        #Map[["log_sigma_age"]] <- NULL
        Map[["logit_rho_time"]] <- NULL
        Map[["log_sigma_time"]] <- NULL
        Map[["phi"]] <- NULL
        #Map[["log_sigma_loc"]] <- NULL
    }
    if (option >= 2){
        Random <- c("epsilon")
        Map[["logit_rho_loc"]] <- NULL
        Map[["log_sigma_loc"]] <- NULL
        Map[["epsilon"]] <- NULL
    }
    dyn.load(dynlib(model_name))
    par_num <- ifelse(model == "gpz", 5, 6)
    Params <- list(gpz=rep(0, par_num), log_sigma_obs=0, logit_rho_time=0,
                   logit_rho_age=0, logit_rho_loc=logit_rho_loc, 
                   log_sigma_time=0, log_sigma_age=0, log_sigma_loc=0, 
                   phi=array(0, dim=c(L_, A_, T_)), epsilon=array(0, dim=c(L_)),
                   secular=0.)
    Data <- list(log_rate_mort=df$log_rate, age=df$age_mean, option=option,
                 age_group=df$age_group, time_group=df$time_group,
                 year_dat=year_dat, graph=admin_queens(prec=F), 
                 loc_group=df$loc_group)
    Obj <- MakeADFun(data=Data, parameters=Params, DLL=model_name,
                     silent=!print, map=Map, random = Random)
    Obj$env$tracemgc <- print
    Obj$env$inner.control$trace <- print
    Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr,)
                  #control=list(iter.max=1000, eval.max=500))
    Report <- Obj$report()
    if(Opt$convergence != 0){
        print("Model had partial or non_convergence!")
    }
    if(Q){
        Report$Q <- sdreport(Obj, getJointPrecision=TRUE)
    }
    dyn.unload(dynlib(model_name))
    Report
}

run_full <- function(df, sex, holdout=2006){
    df <- subset(df, sex_id == sex & age_mean<=100 & year_id>=1990)
    df <- df[order(df$location_id, df$age_group_id, df$year_id),]
    df$age_group <- df$age_group_id - min(df$age_group_id)
    df$time_group <- df$year_id - min(df$year_id)
    df$loc_group <- as.numeric(as.factor(df$location_id)) - 1
    df_edit <- df
    df_edit$log_rate[df_edit$year_id >= holdout] <- NA
    silder <- run_model(df_edit, option=1, print=T, Q=F)
    std_error <- var_AR1(max(df_edit$year_id) - holdout + 1, 
                         silder$rho_time, silder$sigma_time)**.5 * 1.96
    df_edit$log_rate_hat <- silder$log_rate_mort_hat
    df_edit$lwr_bound <- silder$log_rate_mort_hat
    df_edit$lwr_bound[df_edit$year_id >= holdout] <- 
        df_edit$lwr_bound[df_edit$year_id >= holdout] - std_error
    df_edit$upr_bound <- silder$log_rate_mort_hat
    df_edit$upr_bound[df_edit$year_id >= holdout] <- 
        df_edit$upr_bound[df_edit$year_id >= holdout] + std_error
    df_edit$log_rate <- df$log_rate
    
    add_df <- subset(df, year_id > 2000)
    add_df$year_id <- add_df$year_id + 15
    add_df$log_rate <- NA
    add_df$rate <- NA
    add_df$envelope <- NA
    add_df$population <- NA
    
    df2 <- rbind(df, add_df)
    df2$time_group <- df2$year_id - min(df2$year_id)
    df2 <- df2[order(df2$location_id, df2$age_group_id, df2$year_id),]
    
    silder2 <- run_model(df2, option=1, print=T, Q=F)
    df2$log_rate_hat <- silder2$log_rate_mort_hat
    
    loc_df <- get_locations()
    
    df_edit <- left_join(df_edit, loc_df)
    df2 <- left_join(df2, loc_df)
    list(df_edit=df_edit, df2=df2)
}

all_df <- list(male=0, female=0)

all_df$female <- run_full(df, 2, 2006)
names(all_df$female) <- c()
all_df$male <- run_full(df, 1, 2006)

save(all_df, file="~/Documents/fish_proj/all_data.rda")
