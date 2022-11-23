# TODO: Calculate Overlap with Convex Hulls and with Ellipses -> Difference to 
# real data via boxplot?
trait_data_bind <- readRDS(file = file.path(data_cache,
                                            "trait_data_ww_bind.rds"))

# Create copy of data for the use in for-loop (rm continent column to avoid confusion) 
trait_data_bind_sim <- copy(trait_data_bind)
trait_data_bind_sim[, "continent" := NULL]

n_continent <- trait_data_bind[, .N, by = "continent"][, N]
names_continent <- trait_data_bind[, .N, by = "continent"][, continent]

# Establish sampling times and family pool
sampling_times <- 100
family_pool <- trait_data_bind_sim[, id]

sim_data <- list()
for(i in 1:10){
  families_sampled <- sample(family_pool) 
  # Could even apply randomness to the split
  families_sampled <- split(families_sampled, rep(1:5, n_continent))
  names(families_sampled) <- names_continent  
  sim_data[[i]] <- lapply(families_sampled, function(y) trait_data_bind_sim[id %in% y,])
}
# lapply(sim_data[[1]], head)

test <- lapply(sim_data[1:10], function(y)
  rbindlist(y, idcol = "continent")) %>%
  rbindlist(., idcol = "iteration")

test_pcoa <- list()
for(i in unique(test$iteration)[1:5]) {
  test_pcoa[[i]] <-
    calc_pcoa(x = test[iteration == i,]) %>%
    transf_pcoa_dt(.) %>%
    cbind(., "continent" = test[iteration == i, continent]) #%>% 
    # calc_cnx_hull(.)
}
test_hull <- list()
for(i in unique(test$iteration)[1:5]) {
  test_hull[[i]] <-
    calc_pcoa(x = test[iteration == i,]) %>%
    transf_pcoa_dt(.) %>%
    cbind(., "continent" = test[iteration == i, continent]) %>% 
    calc_cnx_hull(.) # for this to work, comment out overlap caclulation in this function
}
ggplot(test_pcoa[[5]], aes(x = A1, y = A2)) +
  geom_point(alpha = 0.15) +
  scale_color_d3(name = "Continent",
                 labels = c("AUS", "EUR", "NA", "NZ", "SA")) +
  geom_polygon(
    data = test_hull[[5]],
    alpha = 0.05,
    aes(color = continent)
  ) +
#  facet_wrap(~continent)
  theme_bw()