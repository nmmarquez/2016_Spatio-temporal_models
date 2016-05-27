library(ggplot2)
library(plyr)
library(dplyr)
library(grid)
library(RMySQL)
library(data.table)


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
             parent_id = 102
             AND location_name NOT IN ("Alaska", "Hawaii")'
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

admin_queens <- function(rho=.7){
    load("./gbd_shapefile/gbd15.rdata")
    usa <- gbd15[gbd15$parent_id == 102 & 
                     !(gbd15$loc_name %in% c("Alaska", "Hawaii")),]
    mat <- poly2adjmat(usa)
    diag(mat) <- 0
    n_delta_i <- rowSums(mat)
    mat <- -1 * mat
    diag(mat) <- n_delta_i
    mat
}


Q_ar1 <- function(N, sigma=1, rho=.7){
    Q <- matrix(0, nrow=N, ncol=N)
    Q[1,1] <- 1 + rho**2
    for(i in 2:N){
        Q[i,i] <- 1 + rho**2
        Q[i-1,i] <- -1 * rho
        Q[i,i-1] <- -1 * rho
    }
    (1 / sigma**2) * Q
}


download_mort_data <- function(location_id=NULL){
    if (is.null(location_id)){
        location_id <- 'SELECT location_id
                        FROM shared.location_hierarchy_history
                        WHERE location_set_version_id= 71 AND 
                        parent_id = 102
                        AND location_name NOT IN ("Alaska", "Hawaii")'
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

inf_term <- function(x, N0=1, lambda=-1, c=0){
    N0 * exp(lambda * x) + c
}

logit_scaled <- function(x, eta, scale){
    1 / (1 + exp(scale * (eta - x)))
}

ya_term <- function(x, eta=mean(x), scale=1, ceiling=1){
    ceiling * logit_scaled(x, eta, scale)
}

sns_term <- function(x, rho, m, b){
    (m * x + b) * logit_scaled(x, eta=rho, scale=3)
}

gpz_log <- function(x, N0, lambda, a, k, h, m, b){
    infant_term(x, N0, lambda) + 
        ya_term(x, a, h, k) * (1 - logit_rho(x, rho)) + 
        sns_term(x, m, b)
}


age_plot <- function(df, sex_id, preds=FALSE){
    sub_df <- df[df$sex_id == sex_id,]
    sub_df <- sub_df[order(sub_df$age_mean, sub_df$year_id),]
    plot.age <- ggplot(sub_df, aes(x=year_id, y=log_rate, 
                                   color=age_mean, group=age_mean)) + 
        geom_line()  + scale_x_continuous("Time") + 
        scale_y_continuous("Log Rate(per 100,000) Mortality")
    plot.age <- plot.age + scale_color_gradientn("Age",colours=rainbow(7)) + 
        theme(legend.margin=unit(-0.02,"npc"), legend.text=element_text(size=8),
              text = element_text(size=20))
    if(preds){
        plot.age <- plot.age + geom_path(data=sub_df, 
                                         aes(x=year_id, y=log_rate_hat, 
                                             colour=age_mean, group=age_mean),
                                         linetype="dashed")
    }
    plot.age
}

time_plot <- function(df, sex_id, preds=FALSE){
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
    if(preds){
        plot.time <- plot.time + geom_path(data=sub_df, 
                                         aes(x=age_mean, y=log_rate_hat, 
                                             colour=year_id, group=year_id),
                                         linetype="dashed")
    }
    plot.time
}

reload_model <- function(model){
    if (file.exists(paste0(model, ".so"))) file.remove(paste0(model, ".so"))
    if (file.exists(paste0(model, ".o"))) file.remove(paste0(model, ".o"))
    if (file.exists(paste0(model, ".dll"))) file.remove(paste0(model, ".dll"))
    compile(paste0(model, ".cpp"))
}