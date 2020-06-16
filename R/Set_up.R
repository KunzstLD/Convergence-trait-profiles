# Location of input data 
data_in <- "./Data"

# Location of R-scripts
data_scr <- "./R"

# Location for output
data_out <- "./Output"

# libraries

# data cleaning and data manipulation
library(dplyr)
library(data.table)

# plotting
library(ggplot2)
library(dendextend) # for handling & plotting dendrograms

# Random Forest
library(rsample)  # data splitting
library(ranger)   # a fast c++ implementation of the random forest algorithm
library(vip)      # visualize feature importance
library(pdp)      # visualize feature effects

# Cluster Analysis
library(cluster)
library(FD) # calculation gower distance

# FuzzyCA and related
library(ade4)
library(vegan)

# source function script
source(file = file.path(data_scr, "functions_used.R"))