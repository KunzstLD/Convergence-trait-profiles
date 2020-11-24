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

# lookup table for ecoregions: 
lookup_ER <- readRDS(file.path(
  getwd(),
  "Comparison_between_climatic_zones",
  "Data",
  "ER_lookup.rds"
))

# Meaning KG classifcation
meaning_KG <- fread(
  file.path(
    getwd(),
    "Comparison_between_climatic_zones",
    "Data",
    "Klimaklassifikationen_Koppen_Geiger.txt"
  )
)

# merge with KG classification
# data.table way:
# lookup[ecoreg_kg, `:=`(KoppenClas = i.KoppenClas,
#                        SecKoppenClas = i.SecKC),
#        on = c(ecoregion = "NAME")]
# base way
lookup <- base::merge(
  x = lookup_ER,
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

# check those who have a different first and second Köppen classification
subset_dfc <- test[test$KoppenClas == "Dfc", c("key_col", 
                                               "ID_AQEM",
                                               "KoppenClas",
                                               "SecKC")] 
subset_bsk <- test[test$KoppenClas == "Bsk", c("key_col", 
                                               "ID_AQEM",
                                               "KoppenClas",
                                               "SecKC")]

# subset trait data
setDT(traits_EU)
traits_dfc <- traits_EU[ID_AQEM %in% unique(subset_dfc$ID_AQEM), ]
traits_bsk <- traits_EU[ID_AQEM %in% unique(subset_bsk$ID_AQEM), ]

# save in list
saveRDS(
  object = list("dfc" = traits_dfc, "bsk" = traits_bsk),
  file = file.path(data_cache, "trait_subsets.rds")
)

#### NEXT STEPS
# How many climate regions do we have? 
# Subset of Trait databases to the respective climate regions
# Species can occur in multiple regions -> First analysis without these species
# Separate Cluster Analysis + RF for these datasets -> adopt R scripts
