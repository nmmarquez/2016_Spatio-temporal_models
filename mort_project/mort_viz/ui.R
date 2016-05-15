rm(list=ls())
library(shiny)
library(shinydashboard)
library(RMySQL)

source("./utilities.R")

locations <- sort(get_locations()$location_name)

header <- dashboardHeader(title='Log Mortality Rates')

body <- dashboardBody(
    fluidRow(
        column(width=12,
               tabBox(
                   id='tabvals',
                   width=NULL,
                   tabPanel('Female',plotOutput('d1'), value=1),
                    tabPanel('Male', plotOutput('d2'), value=2)
               )
        ) 
    )
)

sidebar <- dashboardSidebar(
    selectInput('loc', 'Location', locations, selected = "United States")
)

dashboardPage(
    header,
    sidebar,
    body
)