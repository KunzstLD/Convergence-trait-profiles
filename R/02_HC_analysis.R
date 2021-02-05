# ________________________________________________________
# Analysis HC with approrpiate distance
# for fuzzy coded traits
# TODO: Analyse per Order
# TODO: Use Ade4!
# ________________________________________________________

traits_ww <- readRDS(file.path(
  data_cache,
  "preproc_traits.rds")
)

# Hierarchical clustering
hc_output <- list()

set.seed(1234)
for (region in c("AUS", "EU", "NOA", "NZ")) {
  data <- traits_ww[[region]]

  # Create dist. matrxi using Overlap index Manly ---------------------
  # TODO (check what this approach does)
  dist_mat <- dist.ktab(data, type = "F") %>%
    as.matrix(.)

  # Contribution of grouping features to the global distance
  contrib_global_dist <- kdist.cor(data, type = "F")

  # Clustering ---------------------------------------------------------
  hc_taxa <- hclust(as.dist(dist_mat), method = "ward.D")

  # Get labels of dendrogram
  dend_label <- hc_taxa %>%
    as.dendrogram() %>%
    labels()

  # Optimal number of groups using the gap statistic
  gap <- clusGap(
    x = dist_mat,
    FUN = mycluster_hc,
    K.max = 15,
    B = 500
  )

  optimal_nog <- maxSE(gap$Tab[, "gap"],
    gap$Tab[, "SE.sim"],
    method = "Tibs2001SEmax"
  )

  # Store results --------------------------------------------------------
  hc_output[[region]] <- list(
    "distance_matrix" = dist_mat,
    "contribu_global_dist" = contrib_global_dist,
    "hc_element" = hc_taxa,
    "labels_dendrogram" = dend_label,
    "gap_statistic" = gap,
    "optimal_nog" = optimal_nog
  )
}

# TODO: Optimal number of groups? Gap statstic, works also well when data
# fall into "one cluster" (i.e. indication that there is no cluster structure
# if this is the case)
# -> Read on Gap stats a bit more; check Brown paper for what they did
# TODO: Internal cluster validation to test for the presence of a
# hierarchical structure in the data
# Use Approach from Nemec -> SIGTREE
# Or any other measure? If we find an optimal number of groups
# doesn't this confirm our assumption on the structure of the data
# TODO: Comparison Hierarchical clustering with distance matrix
# Cophenetic correlation?
# Different linkage methods, which to choose?
# TODO: test the influence of ward.D, ward.D2 and
# other, see also Everitt and
# https://cran.r-project.org/web/packages/dendextend/vignettes/Cluster_Analysis.html#animals---attributes-of-animals

saveRDS(
  object = hc_output,
  file = file.path(
    data_cache,
    "hc_output_ww.rds"
  )
)
