# Location of input data
data_in <- "./Data"

# Location of R-scripts
data_scr <- "./R"

# Location of intermediate data
data_cache <- "./Cache"

# Location for output
data_out <- "./Output"

# ---- Libraries ---------------------------------------------------------

# data cleaning and data manipulation
library(dplyr)
library(data.table)
library(Hmisc)
library(zeallot)

# plotting
library(ggplot2)
library(dendextend) # for handling & plotting dendrograms

# Random Forest
library(rsample)  # data splitting
library(ranger)   # a fast c++ implementation of the random forest algorithm
library(vip)      # visualize feature importance
library(pdp)      # visualize feature effects
library(caret)

# Cluster Analysis
library(cluster)
library(FD) # calculation gower distance
library(vegan)

# FuzzyCA and related
library(ade4)
library(vegan)

# source function script
source(file = file.path(data_scr, "functions_used.R"))
