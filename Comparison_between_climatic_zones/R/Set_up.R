# Used libraries
library(foreign)
library(data.table)
library(dplyr)

# Path to .shp/.dbf/.csv etc. files
data_in <- file.path(getwd(),
                     "Comparison_between_climatic_zones", 
                     "Data")

# Path to cache
data_cache <- file.path(getwd(), 
                        "Comparison_between_climatic_zones",
                        "Cache")