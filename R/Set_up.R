# Location of input data 
data_in <- "./Data"

# Location of R-scripts
data_scr <- "./R"

# Location for output
data_out <- "./Output"

# libraries
library(rsample)  # data splitting
library(ranger)   # a fast c++ implementation of the random forest algorithm
library(vip)      # visualize feature importance
library(pdp)      # visualize feature effects
library(ggplot2)
library(dplyr)
library(data.table)
library(cluster)
library(vegan)
library(dendextend) # for handling & plotting dendrograms
library(FD) # calculation gower distance

# source function script
source(file = file.path(data_scr, "functions_used.R"))
