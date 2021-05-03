#___________________________________________________________________________________________________
#### Distribution taxa EU & NOA
# Find those that are restricted to one climate zone
# check those who have a different first and second Köppen classification
#___________________________________________________________________________________________________

# Load European trait data with climate region assignment
trait_subsets_eu <-
  readRDS(file.path(
    data_cache,
    "Cache_comparison_climatic_reg",
    "trait_subsets.rds"
  ))

# Omit data we (for now) do not consider (development)
trait_subsets_eu[, c("dev_hemimetabol", "dev_holometabol") := NULL]

# ---- Climateregion & taxon selection -----------------------------------------
trait_region <- trait_subsets_eu[, .(species, genus, family, order, climateregion)]
trait_region[, spec_per_cr := .N, by = climateregion]

ggplot(trait_region[!is.na(species), ])+
  geom_pointrange(aes(x = factor(climateregion), y = spec_per_cr, ymax = spec_per_cr, ymin = 0))+
  coord_flip()

# Genera per climateregion
genus_cr <- trait_subsets_eu[, .N, by = .(genus, climateregion)]

# Calculate rel. distribution according to major climatic regions
# e.g. species of Accentrella occurred 2 times in arid, 3 times in temperate, 4 times in cold, and 2 in polar
# Hence, arid = 2/11, ...
genus_cr[, rel_distrib := N/sum(N), by = genus]

# Genera per region
genus_cr[, .N , by = climateregion]

# Genera only occurring in one climateregion: 37
genus_cr[rel_distrib == 1, ]

# Genera that occurred to more than 50 % in one climate region
# (according to the trait databases)
# Many taxa seem to have equal occurrences in temperate and cold regions
Hmisc::describe(genus_cr[rel_distrib >= 0.5, ])



genus_cr[rel_distrib < 0.5, .N, by = climateregion]

ggplot(genus_cr, 
       aes(x = rel_distrib, 
           y = N, 
           color = climateregion))+
  geom_point()

# Get distinct (to more or less one climate region) genera and subset
# trait data, subsequently merge classification into climate zones
genera_distinct <- genus_cr[rel_distrib > 0.5, ] %>%
  .[climateregion %in% c("temperate", "cold"), genus]


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

trait_agg


trait_agg <- trait_agg[genus %in% genera_distinct, ]
trait_agg[genus_cr[rel_distrib > 0.5, ],
          `:=`(climateregion = i.climateregion),
          on = "genus"
]

trait_agg[, .N, by = "climateregion"]