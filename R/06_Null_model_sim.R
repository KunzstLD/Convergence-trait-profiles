# __________________________________________________________________________________________________
# Null models ----
# # - shuffling families across continents/regions maintaining the number 
# of families per contients/regions
# Traits are tightly linked (see the work of Poff et al. 2006) and only a small
# fraction of the theoretical amount of trait combinations is realised. Hence, it 
# doesn't make sense to shuffle traits.
# __________________________________________________________________________________________________

# Functional niche overlap:
# (Pianka's overlap(?) or Ellipses?) 
# Aim would be to reject the hypothesis that the observed differences in niche space 
# (perhaps accounting for central position and distribution) are pure chance 
# - i.e. would arise from any distribution of taxa across continents. 
# The null model should therefore shuffle families across continents/regions, taking the numbers into account.

# 1) Global pool with families 
#   - Several families occur on multiple continents/regions. Hence, use together with their 
#     original continent for merging back trait information
# 2) Take random sample from global pool, respect initial numbers 


## Simulate null models ----
trait_data_bind <- readRDS(file = file.path(data_cache,
                                            "trait_data_ww_bind.rds"))
trait_data_bind[, id := paste0(continent, "_", family)]
# Create copy of data for the use in for-loop (rm continent column to avoid confusion) 
trait_data_bind_sim <- copy(trait_data_bind)
trait_data_bind_sim[, "continent" := NULL]

n_continent <- trait_data_bind[, .N, by = "continent"][, N]
names_continent <- trait_data_bind[, .N, by = "continent"][, continent]

# Establish sampling times and family pool
sampling_times <- 1000
family_pool <- trait_data_bind_sim[, id]

sim_data <- list()
for(i in 1:sampling_times){
  families_sampled <- sample(family_pool) 
  # Could even apply randomness to the split
  families_sampled <- split(families_sampled, rep(1:5, n_continent))
  names(families_sampled) <- names_continent  
  sim_data[[i]] <- lapply(families_sampled, function(y) trait_data_bind_sim[id %in% y,])
}
# lapply(sim_data[[1]], head)
sim_data <- lapply(sim_data, function(y)
  rbindlist(y, idcol = "continent")) %>%
  rbindlist(., idcol = "iteration")
# saveRDS(sim_data,
#         file.path(data_cache, "simulated_data_nmodels.rds"))

## Convex hulls ----
# Calc PCoA and overlap of convex hulls
sim_data <- readRDS(file.path(data_cache, "simulated_data_nmodels.rds"))
ovl_results <- list()
# start.time <- Sys.time()
for(i in unique(sim_data$iteration)) {
  ovl_results[[i]] <-
    calc_pcoa(x = sim_data[iteration == i,]) %>%
    transf_pcoa_dt(.) %>%
    cbind(., "continent" = sim_data[iteration == i, continent]) %>% 
    calc_cnx_hull(.)
}
# takes ~ 1 hour 14 minutes
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# saveRDS(ovl_results, file = file.path(data_cache, "ovl_results.rds"))

## Ellipses ----
# Calc PCoA and overlap of ellipses
sim_data <- readRDS(file.path(data_cache, "simulated_data_nmodels.rds"))

ellipse_results <- list()
# start.time <- Sys.time()
for(i in unique(sim_data$iteration)) {
  ellipse_results[[i]] <-
    calc_pcoa(x = sim_data[iteration == i,]) %>%
    transf_pcoa_dt(.) %>%
    cbind(., "continent" = sim_data[iteration == i, continent]) %>% 
    calc_ellipses(.) 
}
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# takes ~ 1 hour 48 minutes
# saveRDS(ellipse_results, file = file.path(data_cache, "ellipses_results.rds"))

# TPGs ----
## Clustering ----
# Delineate TPGs
# TPGs: Same shuffling as for the niche overlap
# How much variance explain the different continents/regions (of the mean trait profiles)?
# If with permuted data more variance is explained than with the real data that could be interpreted
# as an indicator for convergence (however, that's not the case). Random mean trait profiles are
# only explained by a small margin by the different continents/regions. At least this deviates from 
# the real dataset (~ 17 % explained by continents/regions)
# - TODO: Check out cluster approach from PUP?
# trait_CONT <- readRDS(file.path(data_cache, "trait_dat_grp_assig.rds"))
sim_data <- readRDS(file.path(data_cache, "simulated_data_nmodels.rds"))

# Clustering & optimal number of groups
tpgs_sim <- list()
for (i in unique(sim_data$iteration)) {
  x <- sim_data[iteration == i,]
  # get tpgs
  tpgs <- split(x, x[, continent]) %>%
    lapply(., calc_clustering) %>%
    add_tpgs_td(cl_obj = .)
  # merge tpgs back
  tpgs_sim[[i]] <- x[tpgs,
                     group := i.group,
                     on = c(id = "family")]
}
tpgs_sim <- rbindlist(tpgs_sim, idcol = "iteration")
# saveRDS(tpgs_sim, file = file.path(data_cache, "tpgs_sim.rds"))

# Calculate mean trait profiles
tpgs_sim <- readRDS(file.path(data_cache, "tpgs_sim.rds"))
tpgs_sim[, iteration := NULL]
tpgs_sim_lf <- melt(
  tpgs_sim,
  id.vars = c("iteration",
              "continent",
              "family",
              "order",
              "id",
              "group"),
  variable.name = "trait",
  value.name = "affinity"
)
tpgs_sim_lf_mean <- tpgs_sim_lf[, .(mean_profile = mean(affinity)),
                                by = c("iteration",
                                       "continent",
                                       "group",
                                       "trait")]
# saveRDS(tpgs_sim_lf_mean, file.path(data_cache, "mean_tps_simulated.rds"))

## dbRDA with mean trait profiles "real" data ----
trait_CONT <- readRDS(file.path(data_cache, "trait_dat_grp_assig.rds"))
mean_tps_real <-
  trait_CONT[, .(mean_profile = mean(affinity)), by = c("continent", "group", "trait")] 
mean_tps_real_wf <- dcast(mean_tps_real, continent + group ~ trait,
                          value.var = "mean_profile")
# Normalize (not sure if this is necessary, but doesn't hurt if normalise once more)
normalize_by_rowSum(mean_tps_real_wf,
                    non_trait_cols = c("continent", "group"))

# create id and convert to data.frame for distance matrix calculations
mean_tps_real_wf[, id := paste0(continent, "_", group)]

# prepare real value for comparison with simulated values 
expl_var_real <- calc_dbrda(t_data = mean_tps_real_wf)
expl_var_real <- as.data.table(expl_var_real)
setnames(expl_var_real, "expl_var_real", "var_constr")
expl_var_real[, type := "real"]

## dbRDA with mean trait profiles simulated data ----
mean_tps_sim <- readRDS(file.path(data_cache, "mean_tps_simulated.rds"))
mean_tps_sim_wf <- dcast(mean_tps_sim, continent+iteration+group~trait, 
      value.var = "mean_profile")
normalize_by_rowSum(mean_tps_sim_wf,
                    non_trait_cols = c("continent", 
                                       "group", 
                                       "iteration"))
mean_tps_sim_wf[, id := paste0(continent, "_", group)]

# Calc. expl. variance of continents/regions 
expl_var_sim <- list()
for (i in unique(mean_tps_sim_wf$iteration)) {
  expl_var_sim[[i]] <- as.data.table(x = calc_dbrda(mean_tps_sim_wf[iteration == i,]))
}
expl_var_sim <- rbindlist(expl_var_sim, idcol = "iteration")
setnames(expl_var_sim, "V1", "var_constr")
expl_var_sim[, type := "simulated"]
expl_var <- rbind(expl_var_sim, expl_var_real, fill = TRUE)
saveRDS(expl_var, file.path(data_cache, "expl_var.rds"))


## RF on simulated trait datasets ----
tpgs_sim <- readRDS(file.path(data_cache, "tpgs_sim.rds"))
tpgs_sim[, iteration := NULL]
tpgs_sim[, group := as.factor(group)]

mc <- getOption("mc.cores", 8)
# start.time <- Sys.time()
rf_results <- mclapply(unique(tpgs_sim$iteration), function(it)
  calc_rf_tpgs(x = tpgs_sim[iteration == it, ]))
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# saveRDS(rf_results, file.path(data_cache, "sim_rf_results.rds"))
rf_results <- readRDS(file.path(data_cache, "sim_rf_results.rds"))
rf_results_tf <-
  lapply(rf_results, function(y)
    as.data.table(unlist(y), keep.rownames = TRUE))
rf_results_tf <- rbindlist(rf_results_tf, idcol = "iteration")
rf_results_tf[, V1 := sub("most_imp_vars\\.", "", V1)]
rf_results_tf[, continent := sub("(.+)(\\.)(.+)", "\\1", V1)]
rf_results_tf[, V1 := sub("(.+)(\\.)(.+)", "\\3", V1)]
setnames(rf_results_tf,
         c("V1", "V2"),
         c("trait", "imp_score"))
imp_traits_sim <- rf_results_tf[order(-imp_score), .SD[1:5],
                                by = c("iteration", "continent")] %>%
  .[order(iteration, continent),]
saveRDS(imp_traits_sim, file.path(data_cache, "sim_five_most_imp_traits.rds"))
