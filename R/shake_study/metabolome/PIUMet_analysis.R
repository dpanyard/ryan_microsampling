##
no_function()


sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

###load data
load("data/shake_study/metabolome_data_analysis/data_preparation/expression_data")
load("data/shake_study/metabolome_data_analysis/data_preparation/sample_info")
load("data/shake_study/metabolome_data_analysis/data_preparation/variable_info")

sxtTools::setwd_project()
setwd("data/shake_study/metabolome_data_analysis/PIUMet_analysis/")

dim(expression_data)
dim(sample_info)
dim(variable_info)

