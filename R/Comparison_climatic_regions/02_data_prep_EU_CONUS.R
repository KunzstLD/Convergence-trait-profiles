#___________________________________________________________________________________________________
#### Occurrence data European taxa ####
# Merge ecoregions with Köppen Geiger classification
# Show how much information is available
# and how many taxa are restricted to one majore climate zone
# TODO: Add confidence of KGC
#___________________________________________________________________________________________________

# - ecoregions classified to Köppen classification 
ecoreg_kg <- read.dbf(file.path(data_in, "Ecoregions.dbf"))

# convert to data.table
setDT(ecoreg_kg)

# correct KoppenClas
ecoreg_kg[KoppenClas == "Bsk", KoppenClas := "BSk"]

# - trait data EU
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


traits_eu_cr <- traits_eu_cr[, .SD,
                             .SDcols = names(traits_eu_cr) %like% "locom|feed|resp|volt|size|bf|ovip|dev|species|genus|family|order"] %>%
  .[order %in% c(
    "Ephemeroptera",
    "Hemiptera",
    "Odonata",
    "Trichoptera",
    "Coleoptera",
    "Plecoptera",
    "Diptera",
    "Megaloptera",
    "Neuroptera"
  ),]

# trait aggregation to family-lvl

# just return rows where for each trait there is an observation
test <- normalize_by_rowSum(x = traits_eu_cr[, .SD, .SDcols = patterns("locom|feed|resp")],
                            non_trait_cols = c("order",
                                               "family",
                                               "genus",
                                               "species")) %>%
  na.omit(.,
          cols = names(.[, -c("family",
                              "order", 
                              "genus",
                              "family")]))



completeness_trait_data(x = traits_eu_cr[, .SD, .SDcols = patterns("volt|feed|resp|locom")],
                        non_trait_cols =
                          c("order",
                            "family",
                            "genus",
                            "species"))
  
traits_eu_cr[order %in% c(
  "Ephemeroptera",
  "Hemiptera",
  "Odonata",
  "Trichoptera",
  "Coleoptera",
  "Plecoptera",
  "Diptera",
  "Megaloptera",
  "Neuroptera"
),]



# save in list
saveRDS(
  object = traits_eu_cr,
  file = file.path(
    data_cache,
    "Cache_comparison_climatic_reg",
    "trait_eu_cr.rds"
  )
)

#__________________________________________________________________________________________________
#### Occurrence data North American taxa ####
# preprocessing see addingKGclass_to_Conus_data.R
#__________________________________________________________________________________________________
occ_noa <- fread(file.path(data_in, "Conus_data_KoppenGeiger.csv"))
occ_noa <- occ_noa[!is.na(Genus), ]

# Extract major cr
# A - Tropical
# B - Arid
# C - Temperate
# D - Cold
# E - Polar
occ_noa[, main_cr := sub("[a-z]{2}|Sh|Wh|Sk|Wk|T|m|f|w", "", Letter_code)]
occ_noa_mcr <- dcast(occ_noa,
                     Genus ~ main_cr,
                     fun.aggregate = fun_binary_length,
                     value.var = "main_cr")

# Taxa restricted to one climate region
occ_noa_mcr[, sum_occ_mcr := apply(.SD, 1, sum), .SDcols = c("A", "B", "C", "D", "E")]

noa_one_mcr <- occ_noa_mcr[sum_occ_mcr == 1, ]
noa_one_mcr[sum_occ_mcr == 1,
            climateregion := fcase(A == 1,
                                   "tropical",
                                   B == 1,
                                   "arid",
                                   C == 1,
                                   "temperate",
                                   D == 1,
                                   "cold",
                                   E == 1,
                                   "polar")]

#__________________________________________________________________________________________________
#### Summary availability taxa restricted to one climateregion ####
#__________________________________________________________________________________________________

# For European taxa
summary_cr_taxa <- rbind(data.table(
  climateregion = c("Total Nr. taxa",
                    "Nr. taxa with ER assignment"),
  spec_per_cr = c(nrow(traits_EU),
                  length(unique(
                    ecoreg_data_lf$ID_AQEM
                  )))
),
unique(traits_eu_cr[, .(climateregion, spec_per_cr)]))
setnames(
  summary_cr_taxa,
  old = c("climateregion", "spec_per_cr"),
  new = c("Variable", "Nr. Taxa EU")
)

# For NOA Taxa 
# Information for 625 Genera
# Select Genera that are restricted to one climateregion:
# cold 122
# temperate 45
# arid 9
# polar 1
summary_cr_taxa[noa_one_mcr[, .N, by = climateregion],
                `Nr. Taxa NOA` := i.N,
                on = c(Variable = "climateregion")]

# Trait data NOA
trait_NOA <- readRDS(file.path(data_in, "Traits_US_LauraT_pp_harmonized.rds"))
summary_cr_taxa[Variable == "Total Nr. taxa", `Nr. Taxa NOA` := length(unique(trait_NOA[!is.na(genus), genus]))]
summary_cr_taxa <- rbind(
  summary_cr_taxa[1:2,],
  list("Genera with occurrence data", traits_eu_cr, nrow(occ_noa_mcr)),
  summary_cr_taxa[3:6,]
)
saveRDS(
  summary_cr_taxa,
  file.path(
  data_cache,
  "Cache_comparison_climatic_reg",
  "summary_cr_taxa.rds"
))
