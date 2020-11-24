# read in trait subsets
trait_subsets <- readRDS(file.path(data_cache, "trait_subsets.rds"))

# check if all taxa are resolved on species-level
less_splvl <- lapply(trait_subsets, function(y) y[is.na(species), .N])

for(i in names(less_splvl)) {
  if (less_splvl[[i]] == 0) {
    cat("all taxa of", i , "are resolved on species-level \n")
  }
  else{
    cat("Taxa of", i, "are not all resolved on species-level")
  }
}

# Nr of taxa
nr_of_taxa <- lapply(trait_subsets, function(y) y[, .N])