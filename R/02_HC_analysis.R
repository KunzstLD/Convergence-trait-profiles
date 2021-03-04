# ________________________________________________________
# Analysis HC with approrpiate distance
# for fuzzy coded traits
# TODO: Analyse per Order
# TODO: Use Ade4!

# From Legendre:
# hierarchical: 
# "inferior-ranking clusters become members of larger, higher-ranking clusters"
# vs. non-hierarchical: 
#"They produce a single partition that optimizes within-group homogeneity,
# instead of a hierarchical series of partitions optimizing the hierarchical
# attribution of objects to clusters"
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

plot(hc_output$NZ$gap_statistic)
hc_output$NZ$optimal_nog

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
# From Legendre:
# "At each clustering step, Ward’s method finds the pair of objects or clusters whose fusion
# increases as little as possible the sum, over all groups formed so far, of the squared
# distances between objects and cluster centroids; that sum is the total within-group
# sum-of-squares"

saveRDS(
  object = hc_output,
  file = file.path(
    data_cache,
    "hc_output_ww.rds"
  )
)
