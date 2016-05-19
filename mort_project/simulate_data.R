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

df <- subset(df, sex_id == 1 & age_mean<=100 & year_id>=1990)
df$age_group <- df$age_group_id - min(df$age_group_id)
df$time_group <- df$year_id - min(df$year_id)
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
    par_num <- ifelse(model == "gpz", 5, 6)
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


df$log_rate_hat <- silder$log_rate_mort_hat
time_plot(df, 1, preds=T)
sum(-1 * dnorm(df$log_rate, df$log_rate_hat, silder$sigma_obs, log=T))# 188
print(silder$rho)
print(silder$m)
print(silder$b)
print(silder$rho)
print(silder$m)
print(silder$b)

df_sub <- subset(df, location_id == 101)

silder <- run_model(df_sub, option=0, print=T)


df_sub$log_rate_hat <- silder$log_rate_mort_hat
time_plot(df_sub, 1, preds=T)
sum(-1 * dnorm(df$log_rate, df$log_rate_hat, silder$sigma_obs, log=T))# 188

print(silder$rho)
print(silder$m)
print(silder$b)


# fixed effects
rho <- silder$rho # age at which young effects trade off for old effects
c <- silder$c # base level mortality (log_rate) when most protected
b <- silder$b # adjustment term sns death rates

N0 <- silder$b # this value plus c is the instantaneous death rate at birth
lambda <- silder$lambda # the slope of decline for early life death
m <- silder$m # the slope of increase for senesence moratlity

Q_geo <- admin_queens() # precision matrix for location
Q_age <- Q_ar1(length(unique(df$age_group_id))) # AR1 on age groups
Q_time <- Q_ar1(length(unique(df$year_id))) # AR1 on time
Q_geo_age <- kronecker(Q_geo, Q_age) # kronecker product geo and age
Q_geo_time <- kronecker(Q_geo, Q_time) # kronecker product geo and time

Sigma_geo <- ginv(Q_geo) # covariance matrix by taking inverse
Sigma_geo_time <- ginv(Q_geo_time) # covariance matrix by taking inverse
Sigma_geo_age <- ginv(Q_geo_age) # covariance matrix by taking inverse

# simulate random effects with the vcov matrices
epsilon_geo_age <- mvrnorm(n=1, mu=rep(0, nrow(Sigma_geo_age)), Sigma=Sigma_geo_age)
epsilon_geo_time <- mvrnorm(n=1, mu=rep(0, nrow(Sigma_geo_time)), Sigma=Sigma_geo_time)
epsilon_geo <- mvrnorm(n=1, mu=rep(0, nrow(Sigma_geo)), Sigma=Sigma_geo)

# add the random effects on to the data frame
df3 <- df[order(df$location_id, df$age_group_id, df$year_id),]
geo_ran <- data.frame(location_id=unique(df3$location_id), re_geo=epsilon_geo)
geo_time_ran <- expand.grid(location_id=unique(df3$location_id), year_id=unique(df3$year_id))
geo_time_ran$re_geo_time <- epsilon_geo_time
geo_age_ran <- expand.grid(location_id=unique(df3$location_id), age_group_id=unique(df3$age_group_id))
geo_age_ran$re_geo_age <- epsilon_geo_age
df3 <- left_join(left_join(left_join(df3, geo_ran), geo_time_ran), geo_age_ran)

# mean effects from fixed values
df3$log_rate_hat <- (N0 * exp(lambda * df3$age_mean) + c) *
    (1 - logit_scaled(df3$age_mean, rho, 3)) + 
    logit_scaled(df3$age_mean, rho, 3) * (m * df3$age_mean + b)
# random terms added (geo, time, age)
df3$log_rate <- df3$log_rate_hat + df3$re_geo + df3$re_geo_time + df3$re_geo_age
time_plot(subset(df3, location_id == 101), 1, preds=T)
multiplot(time_plot(subset(df3, location_id == 101), 1), time_plot(subset(df3, location_id == 102), 1), 
          time_plot(subset(df3, location_id == 6), 1), time_plot(subset(df3, location_id == 7), 1), cols=2)

load("./gbd_shapefile/gbd15.rdata")
source("~/Documents/r_shared/woodson_mapping_suite/plot_data.R")
source("~/Documents/r_shared/woodson_mapping_suite/woodson_pallettes.R")

series_map(chloropleth_map=gbd15[gbd15@data$level==3,],
           data=data.table(df3),
           geog_id="location_id",
           variable="log_rate",
           map_title="Log Mortality Rate",
           series_dimension="year_id",
           series_sequence=c(1990,2000,2010),
           color_ramp=rev(woodson_pallettes("cool_toned")),
           histogram=TRUE)