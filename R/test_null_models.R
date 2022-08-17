trait_data_bind <- readRDS(file = file.path(data_cache,
                                            "trait_data_ww_bind.rds"))
# test data
trait_data_bind[, id := paste0(continent, "_", family)]
test <-
  trait_data_bind[continent %in% c("AUS", "EU"), .SD[1:5], by = "continent"]


# Continents/regions to loop over
families_p_continent <- test[, .N, by = "continent"]
test[, continent := NULL]

# Vectors & lists for storing and overwriting in the for loop
family_pool <- test[, id]
families_sampled <- NULL
results <- list()

# Specify number of simulations and steps taken in the inner for loop
sampling_times <- 4 # how many simulations
n_continents <-
  length(families_p_continent$continent) # nr of continents
steps <-
  base::split(1:(sampling_times * n_continents),
              cut(1:(sampling_times * n_continents), 4, labels = FALSE))# steps of the inner for loop

# Get a "fair" order of continents and regions
# (i.e. similar times at first, at second, at third, ...)
perm_continents <- unlist(combinat::permn(c("AUS", "EU")))
loop_continents <-
  base::rep(perm_continents,
            sampling_times * n_continents / length(perm_continents))

for (j in 1:sampling_times) {
  sim_data <- list()
  step <- steps[[j]]
  for (i in loop_continents[step]) {
    # Update the family pool after each round so that families do not get assigned twice
    # If family pool is empty, assign original pool
    family_pool <- family_pool[!family_pool %in% families_sampled]
    if (length(family_pool) == 0) {
      family_pool <- test[, id]
    }
    
    # n is the number of families per continent
    n <- families_p_continent[continent == i, N]
    families_sampled <- sample(family_pool, size = n)
    
    # subset trait data
    sim_data[[i]] <- test[id %in% families_sampled,]
  }
  results[[j]] <- sim_data
}
results[[4]]

# Check if some families have been assigned twice
# (TRUE if families are all unique!)
#lapply(results, function(y) all(!y$AUS[, id] %in% y$EU[, id]))
