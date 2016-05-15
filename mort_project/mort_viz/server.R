rm(list=ls())
library(shiny)
library(ggplot2)
library(xlsx)

source("./utilities.R")

loc_df <- get_locations()

shinyServer(function(input,output){
    loc_id <- reactive({subset(loc_df, location_name == input$loc)$location_id})
    dfur <- reactive({download_mort_data(loc_id())})
    
    output$d1 <- renderPlot({
        dfur2 <- dfur()
        multiplot(age_plot(dfur2, 2), time_plot(dfur2, 2), cols=2)},
        width = 1100, height = 900)
    output$d2 <- renderPlot({
        dfur2 <- dfur()
        multiplot(age_plot(dfur2, 1), time_plot(dfur2, 1), cols=2)},
        width = 1100, height = 900)
})