# read in trait subsets
trait_subsets <- readRDS(file.path(data_cache, "trait_subsets.rds"))

# bind together
trait_subsets <- rbindlist(trait_subsets, idcol = "climateregion")

# fetch duplicates
dupl <- trait_subsets[duplicated(ID_AQEM), ID_AQEM]

# rm duplicates
trait_subsets_uq <- trait_subsets[!ID_AQEM %in% dupl, ]

# How often occur taxa in several climate regions?
# prev <- list()
# for(i in seq_along(dupl)) {
#   cr <- test[ID_AQEM == dupl[i], climateregion]
#   prev[[i]] <-
#     list(
#       "Prevalence" = cr %in% c("arid", "temperate", "cold", "polar"),
#       "ID_AQEM" = dupl[i],
#       "climateregion" = cr
#     )
# }
# ind <- lapply(prev, function(y) sum(y$Prevalence) <= 2)
# prev[unlist(ind)]


# check if all taxa are resolved on species-level
# only two are on genus-level
trait_subsets[is.na(species), ]

# Nr of taxa
nr_taxa <- trait_subsets[, .N, by = climateregion]

# How complete is the trait data?


# Which families/orders occur per region?

# save
saveRDS(
  trait_subsets,
  file = file.path(data_cache, "trait_subsets.rds")
)
saveRDS(
  trait_subsets_uq,
  file = file.path(data_cache, "trait_subsets_uq.rds")
)