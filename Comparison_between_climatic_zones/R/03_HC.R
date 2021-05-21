# ________________________________________________________
# Analysis HC with approrpiate distance
# for fuzzy coded traits

# TODO: Analyse per Order?
# TODO: create custom dist function for fuzzy coded and binary variables?

# load data
data <- readRDS(file.path(
  data_cache,
  "preproc_data_genus.rds")
)

# Hierarchical clustering
hc_output <- list()

set.seed(1234)
for (i in c("temperate")) { # , "cold" 

  # Create dist. matrix
  # for now use Orloci's Chord distance
  dist_mat <- decostand(data[[i]], "norm", na.rm = TRUE) %>%
    vegdist(., "euclidean", na.rm = TRUE) %>%
    as.matrix()

  # standardise distances to [0,1]
  dist_mat <- dist_mat/sqrt(2)

  # Some NA values, need to be exclude for some species
  # for which only few trait information is available
  orig_dist_mat <- dist_mat
  dist_mat <- na.omit(dist_mat)
  dist_mat <- dist_mat[, colnames(dist_mat) %in% rownames(na.omit(dist_mat))]

  # Optimal number of groups
  # Using the gap statistic
  gap <- clusGap(
    x = dist_mat,
    FUN = mycluster_hc,
    K.max = 10,
    B = 500
  )

  # Stores optimal nr of groups, gap stat. and dist. matrix
  optimal_nog <- maxSE(gap$Tab[, "gap"],
    gap$Tab[, "SE.sim"],
    method = "Tibs2001SEmax"
  )

  hc_output[[i]] <- list(
    "distance_matrix" = dist_mat,
    "orig_distance_matrix" = orig_dist_mat,
    "gap_statistic" = gap,
    "optimal_nog" = optimal_nog
  )
}

saveRDS(
  object = hc_output,
  file = file.path(
    data_cache,
    "hc_output_temperate.rds"
  )
)