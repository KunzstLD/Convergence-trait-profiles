# _____________________________________________________________________________
# Data preparation
# TODO: col_two_levels actually selects also for columns with less than
# two levels
# _____________________________________________________________________________

# Read in data
trait_subsets <- readRDS(file.path(data_cache, "trait_subsets.rds"))

# Omit data we (for now) do not consider (development)
trait_subsets[, c("dev_hemimetabol", "dev_holometabol") := NULL]

# ----- Trait Aggregation -----------------------------------------------------

# Columns not to consider for aggregation
non_trait_cols <- grep(
      "ER[0-9]{1,}|species|genus|family|order|tax.*|climateregion|ID.*",
      names(trait_subsets),
      value = TRUE
    )

# Trait aggregation to genus lvl (for comparability reasons with NOA data)
trait_agg <- direct_agg(
    trait_data = trait_subsets[!is.na(species), ],
    non_trait_cols = non_trait_cols,
    method = median,
    taxon_lvl = "genus",
    na.rm = TRUE
  )

# Bind back entries that are already resolved on genus-lvl
trait_agg <- rbind(
  trait_agg,
  trait_subsets[genus %in% c("Boreoheptagyia sp.", "Plesioperla sp."), .SD,
    .SDcols = patterns("genus|feed.*|resp.*|volt.*|locom.*|ovip.*|size.*|bf.*")
  ]
)

# ---- Climateregion & taxon selection -----------------------------------------

# Genera per climateregion
genus_cr <- trait_subsets[, .N, by = .(genus, climateregion)]

# Calculate rel. distribution according to major climatic regions
genus_cr[, rel_distrib := N/sum(N), by = genus]

# Genera only occurring in one climateregion: 37
genus_cr[rel_distrib == 1, ]

# Genera that occurred to more than 50 % in one climate region
# (according to the trait databases)
# Many taxa seem to have equal occurrences in temperate and cold regions
Hmisc::describe(genus_cr[rel_distrib > 0.5, ])

# Get distinct (to more or less one climate region) genera and subset
# trait data, subsequently merge classification into climate zones
genera_distinct <- genus_cr[rel_distrib > 0.5, ] %>%
  .[climateregion %in% c("temperate", "cold"), genus]

trait_agg <- trait_agg[genus %in% genera_distinct, ]
trait_agg[genus_cr[rel_distrib > 0.5, ],
  `:=`(climateregion = i.climateregion),
  on = "genus"
]

# ---- Preprocessing HC --------------------------------------------------------

preproc_data <- list()
for (region in c("temperate", "cold")) {
  data <- trait_agg[climateregion == region, ]
  
  # rm size and body form traits for now (only small coverage)
  data[, c(
    "bf_streamlined",
    "bf_cylindrical",
    "bf_spherical",
    "bf_flattened",
    "size_large",
    "size_medium",
    "size_small",
    "climateregion"
  ) := NULL]

  # ---- Data preparation for HC -----------------------------------------------
  # convert to data.frame -> data table does not support row.names
  setDF(data)

  # add row.names
  row.names(data) <- data$genus
  data$genus <- NULL

  # data with only two levels/integer/factor variables have to be numeric
  col_int_fac <-
    names(Filter(function(y) {
      is.integer(y) | is.factor(y)
    }, data))
  data[, col_int_fac] <- apply(data[, col_int_fac], 2, as.double)

  # save to list
  preproc_data[[region]] <- data
}

# save
saveRDS(
  object = preproc_data,
  file = file.path(
    data_cache,
    "preproc_data_genus.rds"
  )
)
