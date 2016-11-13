library(ggplot2)
library(data.table)
setwd("~/Documents/Classes/2016_Spatio-temporal_models/")

rm(list=ls())
source("./mort_project/mort_viz/utilities.R")
load("./mort_project/all_data.rda")

df <- as.data.table(subset(all_df$male$df_edit, year_id > 2011))

df$diff <- df$log_rate - df$log_rate_hat

ggplot(df, aes(factor(age_group), diff)) + 
    geom_boxplot() +  
    labs(title="Forecasted Error Across Age Groups", x="Age Group", y="Error")
ggsave("./mort_project/forecasted_error_by_age.png")

load("./mort_project/gbd_shapefile/gbd15.rdata")
source("~/Documents/r_shared/woodson_mapping_suite/plot_data.R")
source("~/Documents/r_shared/woodson_mapping_suite/woodson_pallettes.R")

usa <- gbd15[gbd15$parent_id == 102 & 
                 !(gbd15$loc_name %in% c("Alaska", "Hawaii")),]

series_map(chloropleth_map=usa,
           data=subset(df, age_group_id == 12),
           geog_id="location_id",
           variable="diff",
           map_title="Log Mortality Rate Error",
           series_dimension="year_id",
           series_sequence=c(2015),
           color_ramp=woodson_pallettes("thanksgiving"),
           histogram=FALSE)
ggsave("./mort_project/log_rate_error_map.png")

time_plot(subset(all_df$female$df_edit, loc_group ==0 & age_group > 6), 2,
          title_="Alabama Log Death Rates: Females")
ggsave("./mort_project/al_female_sns_mort.png")

time_plot(subset(all_df$male$df_edit, location_name =="Texas"), 1,
          title_="Texas Log Death Rates: Males")
ggsave("./mort_project/tx_male_all.png")

age_plot(subset(all_df$male$df_edit, location_name =="California"), 1,
         title_="California Log Death Rates: Males", preds=T, draws=T)
ggsave("./mort_project/ca_forecast_age.png")

time_plot(subset(all_df$male$df2, location_name =="California"), 1,
          title_="California Log Death Rates: Males", preds=T)
ggsave("./mort_project/ca_forecast_time.png")
