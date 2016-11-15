rm(list=ls())
library(dplyr)

setwd("~/Documents/Classes/2016_Spatio-temporal_models/mort_project/")
source("./mort_viz/utilities.R")

if (!("data.rda" %in% list.files())){
    df <- download_mort_data()
    save(df, file="./data.rda")
}

load("./data.rda")
df <- subset(df, age_mean<=100 & year_id>=1990)


rw_nll <- function(params, data){
    nll <- 0.0
    d <- params[1]
    sigma <- exp(params[2])
    for(i in 2:length(data)){
        est <- data[i-1] + d
        nll <- nll - dnorm(est, data[i], sigma, TRUE)
    }
    return(nll)
}

random_walk <- function(time_series){
    opt <- optim(c(0, 0), rw_nll, data=time_series)
    return(list(d=opt$par[1], sigma=exp(opt$par[2])))
}

lee_carter_model <- function(sex, year_id){
    y <- year_id
    df_ <- subset(df, sex_id == sex & year_id < y)
    df_ <- df_[order(df_$location_id, df_$age_group_id, df_$year_id),]
    locs <- unique(df_$location_id)
    ages <- unique(df_$age_group_id)
    T_ <- length(unique(df_$year_id))
    A_ <- length(unique(ages))
    L_ <- length(unique(locs))
    
    M <- aperm(array(df_$log_rate, dim=c(T_, A_, L_)), c(3,2,1))
    A_cx <- apply(M, c(1,2), mean)
    svds <- lapply(1:L_, function(c) svd(M[c,,] - A_cx[c,]))
    B_cx <- t(sapply(1:L_, function(c) svds[[c]]$u[,1]))
    K_ct <- t(sapply(1:L_, function(c) svds[[c]]$v[,1]))
    rw_params <- apply(K_ct, 1, random_walk)
    
    d_c <- sapply(rw_params, function(x) x$d)
    df2 <- expand.grid(year_id=1990:2040, age_group_id=ages, location_id=locs)
    keep_cols <- c("location_id", "age_group_id", "year_id", "log_rate")
    df2 <- left_join(df2, subset(df, sex_id == sex)[,keep_cols])
    df_last <- subset(df2, year_id == y)
    df_last$log_rate_last <- df_last$log_rate
    df_last$log_rate <- NULL
    df_last$year_id <- NULL
    df2$d_multiply <- df2$year_id - year_id
    df2 <- left_join(df2, data.frame(d=d_c, location_id=locs))
    b_values <- expand.grid(location_id=locs, age_group_id=ages)
    b_values$b <- c(B_cx)
    df2 <- left_join(df2, b_values)
    df2 <- left_join(df2, df_last)
    df2$log_rate_hat <- df2$d_multiply * df2$d * df2$b + df2$log_rate_last
    df2
}

run_full <- function(sex, holdout){
    loc_df <- get_locations()
    df1 <- left_join(lee_carter_model(sex, holdout), loc_df)
    df2 <- left_join(lee_carter_model(sex, 2015), loc_df)
    list(df_edit=df1, df2=df2)
}

all_df <- list(male=0, female=0)

all_df$female <- run_full(2, 2006)
all_df$male <- run_full(1, 2006)

save(all_df, file="./lee_carter_all_data.rda")
