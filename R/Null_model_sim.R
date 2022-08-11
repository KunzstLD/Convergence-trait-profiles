###################################################################################################
# Script for calculating null models
# # - shuffling families across continents/regions maintaining the number 
# of families per contients/regions
# Traits are tightly linked (see the work of Poff et al. 2006) and only a small
# fraction of the theoretical amount of trait combinations is realised. Hence, it 
# doesn't make sense to shuffle traits.
###################################################################################################

# Functional niche overlap:
# Pianka's overlap(?) or Ellipses? 
# Aim would be to reject the hypothesis 
# that the observed differences in niche space (perhaps accounting for central position and distribution) 
# are pure chance 
# - i.e. would arise from any distribution of taxa across continents. 
# The null model should therefore shuffle families across continents/regions, taking the numbers into account.

# 1) Global pool with families 
#   - Several families occur on multiple continents/regions. Hence, use together with their 
# orginal continent for mergining back trait information
# 2) Take random sample from global pool, respect initial numbers 
trait_data_bind <- readRDS(file = file.path(data_cache,
                                            "trait_data_ww_bind.rds"))
trait_data_bind[, id := paste0(continent, "_", family)]

test <- trait_data_bind[continent %in% c("AUS", "EU"), .SD[1:5], by = "continent"]


n_families <- test[, .N, by = "continent"]
test[, continent := NULL]
family_pool <- test[, id]

families_sampled <- NULL
sim_data <- list()

# Sample the global family pool
for(i in n_families$continent) {
  # Update the family pool after each round so that families do not get assigned twice
  family_pool <- family_pool[!family_pool %in% families_sampled]
  n <- n_families[continent == i, N]
  families_sampled <- sample(family_pool, size = n)
  
  # subset trait data
  sim_data[[i]] <- test[id %in% families_sampled,]
}

# Apply clustering
# Calculate Overlap (no PCOA?)
# Delineate TPGs


# TPGs:
# Same shuffling as for the niche overlap
# RDA; Cluster approach from PUP?
# Mean trait profile
# trait_CONT <- readRDS(file.path(data_cache, "trait_dat_grp_assig.rds"))
tpgs_mean <-
  trait_CONT[, .(mean_profile = mean(affinity)), by = c("continent", "group", "trait")] 


# Analysis: climate zones vs. number of trait profiles per family