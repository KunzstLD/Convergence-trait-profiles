# Location of input data
data_in <- "./Data"

# Location of R-scripts
data_scr <- "./R"

# Location of intermediate data
data_cache <- "./Cache"

# Location for output
data_out <- "./Output"

# Location for paper related output 
data_paper <- "./Paper"


# ---- Libraries -----------------------------------------------------------------------------------

# Data cleaning and data manipulation
library(foreign)
library(dplyr)
library(data.table)
library(Hmisc)
library(zeallot)

# Plotting
library(ggplot2)
library(dendextend) # for handling & plotting dendrograms
library(patchwork)
library(adegraphics)

# Table outputs
library(knitr)
library(reactable)

# Random Forest
library(mlr3)
library(mlr3viz)
library(mlr3learners)
library(mlr3measures)
library(mlr3tuning)
library(caret)

# Variable importance
library(Boruta)

# Cluster Analysis, treatment fuzzy coded variables and related
library(cluster)
library(FD)
library(vegan)
library(ade4)
library(vegan)

# Spatial operations
library(sf)
library(maptools)
library(raster)
library(sp)
library(tiff)
library(ijtiff)
library(rgdal)

# source function script
source(file = file.path(data_scr, "functions_used.R"))
