library(ggplot2)
library(plyr)
library(dplyr)
library(grid)
library(RMySQL)
library(data.table)

source('./db_access.R')


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    plots <- c(list(...), plotlist)
    numPlots <- length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

get_locations <- function(){
    call <- 'SELECT location_id, location_name 
             FROM shared.location_hierarchy_history
             WHERE location_set_version_id= 71 AND 
             location_type="admin0"'
    query_cod_db(call)
}

age_df_download <- function(){
    call <- '
    SELECT age_group_id,
    (age_group_years_start + age_group_years_end) / 2 as age_mean
    FROM shared.age_group
    WHERE age_group_id BETWEEN 2 AND 21;
    '
    query_cod_db(call)
}


download_mort_data <- function(location_id=NULL){
    if (is.null(location_id)){
        location_id <- 'SELECT location_id 
                        FROM shared.location_hierarchy_history
                        WHERE location_set_version_id= 71 AND 
                        location_type="admin0"'
    }
    call <- '
    SELECT sex_id, location_id, age_group_id, year_id, mean_pop as population,
    upper_env_hivdeleted as envelope
    FROM mortality.output WHERE location_id IN (%location_id%)
    AND output_version_id = (SELECT output_version_id 
    	                     FROM mortality.output_version WHERE is_best = 1)
    AND sex_id IN (1, 2)
    AND age_group_id BETWEEN 2 AND 21;
    '
    call <- gsub("%location_id%", location_id, call)
    df <- query_cod_db(call)
    df$rate <- df$envelope / df$population * 10**5
    df$log_rate <- log(df$rate)
    df <- left_join(df, age_df_download())
    df
}


age_plot <- function(df, sex_id){
    sub_df <- df[df$sex_id == sex_id,]
    sub_df <- sub_df[order(sub_df$age_mean, sub_df$year_id),]
    plot.age <- ggplot(sub_df, aes(x=year_id, y=log_rate, 
                                   color=age_mean, group=age_mean)) + 
        geom_line()  + scale_x_continuous("Time") + 
        scale_y_continuous("Log Rate(per 100,000) Mortality")
    plot.age <- plot.age + scale_color_gradientn("Age",colours=rainbow(7)) + 
        theme(legend.margin=unit(-0.02,"npc"), legend.text=element_text(size=8),
              text = element_text(size=20))
    plot.age
}

time_plot <- function(df, sex_id){
    sub_df <- df[df$sex_id == sex_id,]
    sub_df <- sub_df[order(sub_df$age_mean, sub_df$year_id),]
    plot.time <- ggplot(sub_df, aes(x=age_mean, y=log_rate, color=year_id, 
                                    group=year_id)) +
        geom_line() + scale_x_continuous("Age") +
        scale_y_continuous("Log Rate(per 100,000) Mortality")
    plot.time <- plot.time + scale_color_gradientn("Time", colours=rainbow(7)) + 
        theme(legend.justification=c(0,1), legend.position=c(0.05,1),
              legend.direction="horizontal", legend.text=element_text(angle=45),
              legend.title.align=1,
              legend.background = element_rect(fill="transparent"),
              text = element_text(size=20))
    plot.time
}

reload_model <- function(model){
    if (file.exists(paste0(model, ".so"))) file.remove(paste0(model, ".so"))
    if (file.exists(paste0(model, ".o"))) file.remove(paste0(model, ".o"))
    if (file.exists(paste0(model, ".dll"))) file.remove(paste0(model, ".dll"))
    compile(paste0(model, ".cpp"))
}