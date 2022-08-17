test <- trait_data_bind[continent %in% c("AUS", "EU"), .SD[1:5], by = "continent"]

n_families <- test[, .N, by = "continent"]
test[, continent := NULL]


# Vectors for loop
family_pool <- test[, id]
families_sampled <- NULL
sim_data <- list()
results <- list()

# Prepartions for "fair" assignment
# Ensure that famlies for each continent are sampled equally 
# (i.e. similar times at first, at second, at third, ...)
sampling_times <- 100
perm_continents <- unlist(combinat::permn(c("AUS", "EU")))
loop_continents <- base::rep(perm_continents, sampling_times/length(perm_continents))

for(j in 1:sampling_times) {
  for (i in loop_continents) {
    # Update the family pool after each round so that families do not get assigned twice
    # If family pool is empty, assign original pool
    family_pool <- family_pool[!family_pool %in% families_sampled]
    if (length(family_pool) == 0) {
      family_pool <- test[, id]
    }
    
    # n is the number of families per continent
    n <- n_families[continent == i, N]
    families_sampled <- sample(family_pool, size = n)
    
    # subset trait data
    sim_data[[i]] <- trait_data_bind_sim[id %in% families_sampled,]
  }
  results[[j]] <- sim_data
}

