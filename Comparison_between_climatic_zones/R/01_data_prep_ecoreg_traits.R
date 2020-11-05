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

# melt is a data.table function but also works with data.frames 
# in case traits_EU is a data.frame (ignore warning for now)
traits_EU_lf <- melt(traits_EU, measure.vars = er_cols, variable = "key_col")

# subset only to species that have a classification in ecoregions
# data.table way:
# traits_EU_lf_sb <- traits_EU_lf[!is.na(value), ]
# base way:
traits_EU_lf_sb <- traits_EU_lf[!is.na(traits_EU_lf$value), ]

# merge with lookup
# data.table way:
# traits_EU_lf_sb[test, `:=`(ecoregion = i.ecoregion,
#                            KoppenClas = i.KoppenClas,
#                            SecKoppenClas = i.SecKoppenClas),
#                 on = "key_col"]
# base:
traits_EU_lf_sb <- base::merge(x = traits_EU_lf_sb,
                               y = lookup,
                               by = "key_col",
                               all.x = TRUE)

# TODO: transform back to wide format
# base:
# something like this could be the solution, but we need to include all variables
# could also create a subset with just the ecoregions and ID_AQEM and then merge back
# to the trait data
test <- traits_EU_lf_sb[!is.na(traits_EU_lf_sb$ecoregion), ]
reshape2::dcast(test, ID_AQEM + KoppenClas + SecKC ~ key_col)








