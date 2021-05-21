# _____________________________________________________________________
#### Comparison between climatic regions NOA ####
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


# Load trait data
trait_NOA <- readRDS(file.path(data_in, "Traits_US_LauraT_pp_harmonized.rds"))

# Use genus level information from North American trait database
trait_NOA_genus <- trait_NOA[is.na(species) & !is.na(genus), ]
# 32 genera not coverd in trait DB
traits_noa_clim <- merge(
    noa_one_mcr,
    trait_NOA_genus,
    by.x = "Genus",
    by.y = "genus"
)

# select rel columns
# just select aquatic insects
traits_noa_clim <- traits_noa_clim[, .SD,
                                   .SDcols = patterns(
                                       "Genus|family|order|climateregion|feed.*|locom.*|size.*|resp.*|volt.*|ovip*|bf.*"
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

# save
saveRDS(
    object = traits_noa_clim,
    file = file.path(
        data_cache,
        "Cache_comparison_climatic_reg",
        "noa_glvl.rds"
    )
)
