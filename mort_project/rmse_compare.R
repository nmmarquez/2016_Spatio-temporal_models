rm(list=ls())
library(dplyr)

setwd("~/Documents/Classes/2016_Spatio-temporal_models/mort_project/")

results <- list(siler=NA, lc=NA)

load("./all_data.rda")

results$siler <- all_df

load("./lee_carter_all_data.rda")
results$lc <- all_df

rmse_sex <- function(sex, start=2009, end=2011){
    sexn <- ifelse(sex == 1, "male", "female")
    rmse <- c(siler=NA, lc=NA)
    for(model in names(rmse)){
        df <- results[[model]][[sexn]]$df_edit
        df <- subset(df, year_id >= start & year_id <= end)
        N <- nrow(df)
        error <- sqrt(sum((df$log_rate_hat - df$log_rate)**2 / N))
        rmse[model] <- error
    }
    return(rmse)
}

rmse_all <-function(){
    years <- list(near=c(2009, 2011), distant=c(2013, 2015))
    sex_ids <- list(male=1, female=2)
    lapply(sex_ids, function(x) t(sapply(years, function(y)
        rmse_sex(x, y[1], y[2]))))
}

head(results$lc$male$df_edit)
head(results$siler$male$df_edit)
