# libraries
library(foreign)
library(data.table)
library(dplyr)

# path to .shp/.dbf files
# needs to be modified by user
data_in <- file.path(getwd(),
                     "Comparison_between_climatic_zones", 
                     "Data")

#### Load data ####

# ecoregions classified to Köppen classifikation 
ecoreg_kg <- read.dbf(file.path(data_in, "Ecoregions.dbf"))

# convert to data.table
setDT(ecoreg_kg)

# trait data
# loaded as .rds (already in data.table format)
# traits_EU <- readRDS(file.path(data_in, "Trait_freshecol_2020_pp_harmonized_ecoregions.rds"))
# loaded as .csv
traits_EU <- read.csv(file.path(data_in, "Trait_freshecol_2020_pp_harmonized_ecoregions.csv"))

# lookup table: 
# lookup <- data.table(ecoregion = c("Tundra", "Taiga"), 
#                    key_col = c("ER21", "ER23"))
lookup <- data.frame(ecoregion = c("Tundra", "Taiga"), 
                   key_col = c("ER21", "ER23"))

# merge with KG classification
# data.table way:
# lookup[ecoreg_kg, `:=`(KoppenClas = i.KoppenClas,
#                        SecKoppenClas = i.SecKC),
#        on = c(ecoregion = "NAME")]
# base way
lookup <- base::merge(
  x = lookup,
  y = ecoreg_kg[, c("KoppenClas", "SecKC", "NAME")],
  by.x = "ecoregion",
  by.y = "NAME",
  all.x = TRUE
)

# get ecoregions columns
er_cols <- grep("ER[0-9]{1,}", names(traits_EU), value = TRUE)

# select from trait data taxa and ecoregions
ecoreg_data <- traits_EU[, c(er_cols, "ID_AQEM")]

# melt is a data.table function but also works with data.frames 
# in case traits_EU is a data.frame (ignore warning for now)
ecoreg_data_lf <- melt(ecoreg_data, 
                       measure.vars = er_cols,
                       variable = "key_col")

# subset only to species that have a classification in ecoregions
# data.table way:
# traits_EU_lf_sb <- traits_EU_lf[!is.na(value), ]
# base way:
ecoreg_data_lf <- ecoreg_data_lf[!is.na(ecoreg_data_lf$value), ]

# merge with lookup
# data.table way:
# traits_EU_lf_sb[test, `:=`(ecoregion = i.ecoregion,
#                            KoppenClas = i.KoppenClas,
#                            SecKoppenClas = i.SecKoppenClas),
#                 on = "key_col"]
# base:
ecoreg_data_lf <- base::merge(x = ecoreg_data_lf,
                              y = lookup,
                              by = "key_col",
                              all.x = TRUE)

# TODO: transform back to wide format
# base:
# something like this could be the solution, but we need to include all variables
# could also create a subset with just the ecoregions and ID_AQEM and then merge back
# to the trait data
test <- ecoreg_data_lf[!is.na(ecoreg_data_lf$ecoregion), ]

# test <- reshape2::dcast(test, ID_AQEM + KoppenClas + SecKC ~ key_col)

# check those who have a different first and second Köppen classification
subset_dfc <- test[test$KoppenClas == "Dfc", c("key_col", "ID_AQEM", "KoppenClas")] 

setDT(traits_EU)
traits_EU[ID_AQEM %in% unique(subset_dfc$ID_AQEM), ]

# How many regions do we have? 
# Species can occur in multiple regions, how to deal with that?
# lookup table ecoregion & key 





