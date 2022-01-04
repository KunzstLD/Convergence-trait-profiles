# _____________________________________________________________________
# Comparison between climatic regions NOA ----
# _____________________________________________________________________
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
    value.var = "main_cr"
)

# Taxa restricted to one climate region
occ_noa_mcr[, sum_occ_mcr := apply(.SD, 1, sum),
    .SDcols = c("A", "B", "C", "D", "E")
]

noa_one_mcr <- occ_noa_mcr[sum_occ_mcr == 1, ]
noa_one_mcr[
    sum_occ_mcr == 1,
    climateregion := fcase(
        A == 1,
        "tropical",
        B == 1,
        "arid",
        C == 1,
        "temperate",
        D == 1,
        "cold",
        E == 1,
        "polar"
    )
]
noa_one_mcr[, .N, by = climateregion]

# Taxa with possible information
noa_genera_temp <- noa_one_mcr[climateregion == "temperate", ]
noa_genera_cold <- noa_one_mcr[climateregion == "cold", Genus]
saveRDS(noa_genera_temp, 
        file.path(data_cache, "noa_genera_temp.rds"))

# _____________________________________________________________________
## Match trait data ----
# _____________________________________________________________________

# Load trait data
trait_NOA <- readRDS(file.path(data_in, "Traits_US_LauraT_pp_harmonized.rds"))

# Use genus level information from North American trait database
trait_NOA_genus <- trait_NOA[is.na(species) & !is.na(genus), ]
trait_NOA_spec <- trait_NOA[!is.na(species) & !genus %in% trait_NOA_genus$genus, ]

# 32 genera not coverd in trait DB
# noa_one_mcr[!Genus %in% trait_NOA$genus, ]

# data on genus-level
traits_noa_clim <- merge(
    noa_one_mcr,
    trait_NOA_genus,
    by.x = "Genus",
    by.y = "genus"
)
# traits_noa_clim[, .N, by = climateregion]

# data on species-level which needs to be aggregated to genus-level
traits_noa_clim_spec <- merge(noa_one_mcr, 
      trait_NOA_spec, 
      by.x = "Genus", 
      by.y = "genus")
traits_noa_clim_agg <- direct_agg(
    trait_data = traits_noa_clim_spec,
    non_trait_cols = c(
        "Genus",
        "A",
        "B",
        "C",
        "D",
        "E",
        "sum_occ_mcr",
        "climateregion",
        "unique_id",
        "family",
        "order",
        "species"
    ),
    method = median,
    taxon_lvl = "Genus"
)


# _____________________________________________________________________
## Final data prep ----
# select rel columns
# just select aquatic insects
# _____________________________________________________________________
traits_noa_clim <- traits_noa_clim[, .SD,
                                   .SDcols = patterns(
                                       "Genus|family|order|climateregion|feed.*|locom.*|size.*|resp.*|volt.*|bf.*|ovip|dev"
                                   )] %>% 
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
setnames(traits_noa_clim, "Genus", "genus")
newcolorder <- c(
    "genus",
    "family",
    "order",
    "climateregion",
    grep("feed", names(traits_noa_clim), value = TRUE),
    grep("resp", names(traits_noa_clim), value = TRUE),
    grep("volt", names(traits_noa_clim), value = TRUE),
    grep("locom", names(traits_noa_clim), value = TRUE),
    grep("ovip", names(traits_noa_clim), value = TRUE),
    grep("bf", names(traits_noa_clim), value = TRUE),
    "size_small",
    "size_medium",
    "size_large"
)
setcolorder(x = traits_noa_clim, 
            newcolorder)

# prep for aggregated data
# merge information on family and order back
traits_noa_clim_agg[traits_noa_clim_spec,
                    `:=`(order = i.order, 
                         family = i.family,
                         climateregion = i.climateregion),
                    on = "Genus"
]
setnames(traits_noa_clim_agg, "Genus", "genus")

# finally, add data on genus-level and aggregated data together
traits_noa_clim <- rbind(traits_noa_clim,
                         traits_noa_clim_agg)
traits_noa_clim[, .N, by = climateregion]

# save
saveRDS(
    object = traits_noa_clim,
    file = file.path(
        data_cache,
        "Cache_comparison_climatic_reg",
        "noa_glvl.rds"
    )
)
