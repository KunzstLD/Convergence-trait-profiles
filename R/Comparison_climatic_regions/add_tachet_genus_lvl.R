# TODO: Check Plesioperla & Plesioperla sp.
# TODO: Add information on resistance
completeness_trait_data(x = traits_eu_cr[order == "Diptera", 
                                         .SD, 
                                         .SDcols = patterns("size|feed|resp|locom|volt|temp")],
                        non_trait_cols =
                          c("order",
                            "family",
                            "genus",
                            "species"))

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

# Add cr info and taxonomic info
eu_glvl[traits_eu_cr,
        `:=`(climateregion = i.climateregion,
             family = i.family,
             order = i.order),
        on = "genus"]

# Add size from Tachet (genus entries in trait_EU is tachet DB)
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

# nr observations per climateregion
eu_glvl[, .N, by = climateregion]

completeness_trait_data(x = eu_glvl[,.SD, .SDcols = patterns("size|feed|resp|locom|volt|temp|ovip")])

na.omit(eu_glvl[, .SD, .SDcols = patterns("feed|resp|temp|climateregion|genus|family|order")])

# size
