rm(list=ls())
set.seed(123)
library(TMB)
library(knitr)
setwd("~/Documents/Classes/2016_Spatio-temporal_models/Week 5 -- 1D spatial models/Homework/")

model_name <- "hw3"
if (file.exists(paste0(model_name, ".so"))) file.remove(paste0(model_name, ".so"))
if (file.exists(paste0(model_name, ".o"))) file.remove(paste0(model_name, ".o"))
if (file.exists(paste0(model_name, ".dll"))) file.remove(paste0(model_name, ".dll"))
compile(paste0(model_name, ".cpp"))