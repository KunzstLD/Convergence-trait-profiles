# Data processing
library(foreign)
library(data.table)
library(dplyr)
library(Hmisc)

# Cluster analysis
library(cluster)
library(vegan)
library(dendextend)

# Random forest
library(ranger)
library(vip)
library(caret)


# Path to .shp/.dbf/.csv etc. files
data_in <- file.path(
    getwd(),
    "Comparison_between_climatic_zones",
    "Data"
)

# Path to cache
data_cache <- file.path(
    getwd(),
    "Comparison_between_climatic_zones",
    "Cache"
)

# Path to other R scripts
data_scr <- file.path(
    getwd(),
    "Comparison_between_climatic_zones",
    "R"
)

# Path to graphics 
data_grp <- file.path(
    getwd(),
    "Comparison_between_climatic_zones",
    "Graphs"
)

# functions script
source("./Comparison_between_climatic_zones/R/functions.R")