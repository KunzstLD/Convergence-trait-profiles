#__________________________________________________________________________________________________
#### Load data ####
# TODO: Add confidence of KGC
#__________________________________________________________________________________________________

# - ecoregions classified to Köppen classification 
ecoreg_kg <- read.dbf(file.path(data_in, "Ecoregions.dbf"))

# convert to data.table
setDT(ecoreg_kg)

# correct KoppenClas
ecoreg_kg[KoppenClas == "Bsk", KoppenClas := "BSk"]

# - trait data
# loaded as .rds (already in data.table format)
# traits_EU <- readRDS(file.path(data_in, "Trait_freshecol_2020_pp_harmonized_ecoregions.rds"))
# loaded as .csv
traits_EU <- fread(file.path(data_in, "Trait_freshecol_2020_pp_harmonized_ecoregions.csv"))

# - lookup table for ecoregions: 
lookup_ER <- readRDS(file.path(
  getwd(),
  "Comparison_between_climatic_zones",
  "Data",
  "ER_lookup.rds"
))
setDT(lookup_ER)

# - Meaning KG classification
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
lookup_ER[ecoreg_kg, `:=`(KoppenClas = i.KoppenClas,
                          SecKoppenClas = i.SecKC),
          on = c(ecoregion = "NAME")]

# get ecoregions columns
er_cols <- grep("ER[0-9]{1,}", names(traits_EU), value = TRUE)

# select from trait data taxa and ecoregions
ecoreg_data <- traits_EU[, c(er_cols, "ID_AQEM"), with = FALSE]

# melt is a data.table function but also works with data.frames 
ecoreg_data_lf <- melt(ecoreg_data, 
                       measure.vars = er_cols,
                       variable = "key_col")

# subset only to species that have a classification in ecoregions
ecoreg_data_lf <- ecoreg_data_lf[!is.na(value), ]

# merge with lookup
# data.table way:
ecoreg_data_lf[lookup_ER, `:=`(ecoregion = i.ecoregion,
                               KoppenClas = i.KoppenClas,
                               SecKoppenClas = i.SecKoppenClas),
               on = "key_col"]

# create subset columns
# (B) Arid
# (C) Temperate  
# (D) Cold
# (ET) Polar 
ecoreg_data_lf[, subset_1 := sub("[a-z]{2}|Sk", "",KoppenClas)]

# subset trait data
traits_arid <-
  traits_EU[ID_AQEM %in% ecoreg_data_lf[subset_1 == "B", ID_AQEM],]
traits_temperate <-
  traits_EU[ID_AQEM %in% ecoreg_data_lf[subset_1 == "C", ID_AQEM],]
traits_cold <-
  traits_EU[ID_AQEM %in% ecoreg_data_lf[subset_1 == "D", ID_AQEM],]
traits_polar <-
  traits_EU[ID_AQEM %in% ecoreg_data_lf[subset_1 == "ET", ID_AQEM],]

# save in list
saveRDS(
  object = list(
    "arid" = traits_arid,
    "temperate" = traits_temperate,
    "cold" = traits_cold,
    "polar" = traits_polar
  ),
  file = file.path(data_cache, "trait_subsets.rds")
)

#### NEXT STEPS
# How many climate regions do we have? 4
# Subset of Trait databases to the respective climate regions
# check those who have a different first and second Köppen classification
# Species can occur in multiple regions -> Two analysis
# Separate Cluster Analysis + RF for these datasets -> adopt R scripts
