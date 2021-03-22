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
traits_EU <- readRDS(file.path(data_in, "Trait_freshecol_2020_pp_harmonized_ecoregions.rds"))

# - lookup table for ecoregions: 
lookup_ER <- readRDS(file.path(data_in, "ER_lookup.rds"))
setDT(lookup_ER)

# - Meaning KG classification
meaning_KG <- fread(file.path(data_in, "Klimaklassifikationen_Koppen_Geiger.txt"))

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
# with major climate regions
# (B) Arid
# (C) Temperate  
# (D) Cold
# (E) Polar 
ecoreg_data_lf[,  c("main_cr", "main_cr_secclas") := .(sub("[a-z]{2}|Sk|T", "", KoppenClas),
                                                       sub("[a-z]{2}|Sk|T|F|Wk", "", SecKoppenClas))]

# Check where main and second Köppen classification deviate
# Use this as a look up table for later
saveRDS(ecoreg_data_lf[main_cr != main_cr_secclas],
        file = file.path(data_cache, "lookup_dev_main_sec_köppen_clas.rds"))

# Get ID's of taxa that only occur in one major climate region
# Helper function for dcast call's
fun_binary_length <- function(y) {
  as.numeric(ifelse(length(y) == 0, 0, 1))
}
kg_cr <- dcast(ecoreg_data_lf, ID_AQEM ~ main_cr, 
      fun.aggregate = fun_binary_length)
kg_cr[, sum_occ_mcr := apply(.SD, 1, sum), .SDcols = c("B", "C", "D", "E")]

# subset trait data
kg_cr[sum_occ_mcr == 1, climateregion := fcase(B == 1, "arid",
                                               C == 1, "temperate",
                                               D == 1, "cold",
                                               E == 1, "polar")]
traits_EU[kg_cr[!is.na(climateregion), ],
          `:=`(climateregion = i.climateregion),
          on = "ID_AQEM"]

traits_eu_cr <- traits_EU[!is.na(climateregion), ]
traits_eu_cr[, spec_per_cr := .N , by = climateregion]

# save in list
saveRDS(
  object = traits_eu_cr,
  file = file.path(
    data_cache,
    "Cache_comparison_climatic_reg",
    "trait_eu_cr.rds"
  )
)

# Add total nr. of taxa with information for the relevant traits & 
# nr. of taxa with ecoregion assignment
traits_eu_cr[, c("total_nr_taxa",
                 "nr_taxa_er_assignment") := .(nrow(traits_EU),
                                               length(unique(ecoreg_data_lf$ID_AQEM)))]
# For European taxa
unique(traits_eu_cr[, .(climateregion, spec_per_cr, nr_taxa_er_assignment, total_nr_taxa)])

# TODO: Create for NOA Taxa and save as html table 

