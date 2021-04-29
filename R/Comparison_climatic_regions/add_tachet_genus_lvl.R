# TODO: Add information on resistance
# Read in trait data EU with climateregion (cr) assignment. 
# Subset of those that occur only in one cr
# Hypothesis for the comparison between climatic regions: 
# larger in warmer areas
# interstitial, ovip_ter may help against dryness/resistence against droughts
traits_eu_cr <- readRDS(file.path(
  data_cache,
  "Cache_comparison_climatic_reg",
  "trait_eu_cr.rds"
))

# Extract genera that are just restricted to one cr
genus_eu_cr <- traits_eu_cr[, .(genus), by = climateregion] %>%
  dcast(., genus ~ climateregion,
        fun.aggregate = fun_binary_length) 
genus_eu_cr[, sum_occ_mcr := apply(.SD, 1, sum), 
            .SDcols = c("arid", 
                        "cold",
                        "polar",
                        "temperate")]
genera_one_cr <- genus_eu_cr[sum_occ_mcr == 1, genus]

# Aggregate to genus-level
eu_glvl <- direct_agg(
  trait_data = traits_eu_cr,
  non_trait_cols = c("species",
                     "genus",
                     "family",
                     "order",
                     "climateregion"),
  method = median,
  taxon_lvl = "genus"
)
eu_glvl <- eu_glvl[genus %in% genera_one_cr, ]
completeness_trait_data(x = eu_glvl[,.SD,
                                    .SDcols = patterns("size|feed|resp|locom|volt|temp")],
                        non_trait_cols = "genus")

# Delete the Plesioperla sp. (Plesioperla is there with the same information)
eu_glvl <- eu_glvl[genus != "Plesioperla sp.", ]

# Add cr info and taxonomic info
eu_glvl[traits_eu_cr,
        `:=`(climateregion = i.climateregion,
             family = i.family,
             order = i.order),
        on = "genus"]

# __________________________________________________________________________________________________
#### Add traits from tachet for genus lvl aggregated traits ####
# __________________________________________________________________________________________________

# Check for data on genus-lvl in tachet
# and add where trait information is not available
# size
genera_size <- eu_glvl[is.na(temp_eury), genus]

# voltinism
genera_volt <- eu_glvl[is.na(volt_semi), genus]

# feeding mode
genera_feeding <- eu_glvl[is.na(feed_shredder), genus]

# respiration
genera_resp <- eu_glvl[is.na(resp_teg), genus]

# temp preference
genera_temp <- eu_glvl[is.na(temp_eury), genus]

# locomotion
genera_locom <- eu_glvl[is.na(locom_burrow), genus]

# ovip
genera_ovip <- eu_glvl[is.na(ovip_aqu), genus]

# Add size from Tachet (genus entries in trait_EU originate from Tachet)
traits_EU <- readRDS(file.path(data_in, "Trait_freshecol_2020_pp_harmonized_ecoregions.rds"))
eu_glvl[traits_EU[is.na(species) & taxa_adjusted %in% genera_size],
        `:=`(size_small = i.size_small,
             size_medium = i.size_medium,
             size_large = i.size_large),
        on = c(genus = "taxa_adjusted")]

# Add voltinism from Tachet
eu_glvl[traits_EU[is.na(species) & taxa_adjusted %in% genera_volt],
        `:=`(volt_bi_multi = i.volt_bi_multi,
             volt_uni = i.volt_uni,
             volt_semi = i.volt_semi),
        on = c(genus = "taxa_adjusted")]

# Add feeding mode from Tachet
eu_glvl[traits_EU[is.na(species) & taxa_adjusted %in% genera_feeding],
        `:=`(feed_shredder = i.feed_shredder,
             feed_gatherer = i.feed_gatherer,
             feed_predator = i.feed_predator,
             feed_parasite = i.feed_parasite, 
             feed_filter = i.feed_filter, 
             feed_herbivore = i.feed_herbivore),
        on = c(genus = "taxa_adjusted")]

# Add respiration from Tachet
eu_glvl[traits_EU[is.na(species) & taxa_adjusted %in% genera_resp],
        `:=`(resp_gil = i.resp_gil,
             resp_teg = i.resp_teg,
             resp_pls_spi = i.resp_pls_spi),
        on = c(genus = "taxa_adjusted")]

# Add temperature pref from Tachet
eu_glvl[traits_EU[is.na(species) & taxa_adjusted %in% genera_temp,],
        `:=`(temp_cold = i.temp_cold, 
             temp_warm = i.temp_warm,
             temp_eury = i.temp_eury),
        on = c(genus = "taxa_adjusted")]

# Add locomotion from Tachet
eu_glvl[traits_EU[is.na(species) & taxa_adjusted %in% genera_locom,],
        `:=`(locom_burrow = i.locom_burrow,
             locom_crawl = i.locom_crawl,
             locom_sessil = i.locom_sessil,
             locom_swim = i.locom_swim),
        on = c(genus = "taxa_adjusted")]

# Add oviposition from Tachet
eu_glvl[traits_EU[is.na(species) & taxa_adjusted %in% genera_ovip, ],
        `:=`(ovip_aqu = i.ovip_aqu,
             ovip_ovo = i.ovip_ovo,
             ovip_ter = i.ovip_ter),
        on = c(genus = "taxa_adjusted")]

# Nr observations per climateregion
eu_glvl[, .N, by = climateregion]

completeness_trait_data(x = eu_glvl[,.SD, .SDcols = patterns("size|feed|resp|locom|volt|temp|ovip")])

na.omit(eu_glvl[, .SD, .SDcols = patterns("feed|temp|resp|volt|ovip|locom|climateregion|genus|family|order")]) 
na.omit(eu_glvl[, .SD, .SDcols = patterns("feed|resp|climateregion|genus|family|order")]) %>% 
  .[, .N, by = .(climateregion)]
# ? size

# Load data from murria and check availability
murria_traits <- readxl::read_excel(path = file.path(data_in, "murria_trait_table_cp.xlsx"))
setDT(murria_traits)

# At least 64 Genera are covered from Murria
genera_in_eu <- eu_glvl[genus %in% murria_traits$Genus, genus]

# TODO: Check available trait information in murria (e.g. for size)
# for those Genera.
# Harmonise & merge the missing trait information
murria_subset <- murria_traits[Genus %in% genera_in_eu, ]

# Taxa with feeeding mode other discarded (3)
murria_subset <- murria_subset[OTHER == 0, ]

# _________________________________________________________________________
# Harmonize Murria's traits (?)
# TODO: Check if providing SIZE information is actually enough 
# (murria has no information on TEMP preference!)
# _________________________________________________________________________
# Size:
# size_small: size < 9 mm (EU: size < 10 mm)
# size_medium: 9 mm < size > 16 mm (EU: 10 mm < size > 20 mm)
# size_large: size > 16 mm (EU: size > 20 mm)
# _________________________________________________________________________
murria_subset[, size_small := apply(.SD, 1, max, na.rm = TRUE),
         .SDcols = c(
           "SIZE1",
           "SIZE2",
           "SIZE3"
         )
]
murria_subset[, size_large := apply(.SD, 1, max, na.rm = TRUE),
              .SDcols = c("SIZE5",
                          "SIZE6",
                          "SIZE7")]
setnames(murria_subset, "SIZE4", "size_medium")
murria_subset[, c("SIZE1",
                  "SIZE2",
                  "SIZE3",
                  "SIZE5",
                  "SIZE6",
                  "SIZE7") := NULL]

# _________________________________________________________________________

# Feeding mode
murria_subset[, feed_shredder := apply(.SD, 1, max, na.rm = TRUE),
         .SDcols = c(
           "MIN",
           "XYL",
           "SHR"
         )
]
murria_subset[, feed_filter := apply(.SD, 1, max, na.rm = TRUE),
              .SDcols = c("AFF", "PFF")]
setnames(murria_subset,
         old = c("GAT",
                 "PRED",
                 "PAR"),
         new = c("feed_gather",
                 "feed_predator",
                 "feed_parasite"),
         )
murria_subset[, c("MIN",
                  "XYL",
                  "SHR",
                  "AFF",
                  "PFF") := NULL]
# _________________________________________________________________________

# Respiration
murria_subset[, resp_pls_spi := apply(.SD, 1, max, na.rm = TRUE),
              .SDcols = c("SPIR", "PLAS")]
setnames(
  murria_subset,
  old = c("TEG",
          "GILL"),
  new = c("resp_teg",
          "resp_gil")
)
murria_subset[, c("SPIR",
                  "PLAS") := NULL]

# _________________________________________________________________________

# Locom

# _________________________________________________________________________

# Volt 

# _________________________________________________________________________

# Temp

# _________________________________________________________________________

# Ovip


# Normalisation



