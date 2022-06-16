# ___________________________________________________________________________
# Hierarchical cluster analysis ----

# Re HCA (from Legendre)
# hierarchical:
# "inferior-ranking clusters become members of larger, higher-ranking clusters"
# vs. non-hierarchical:
# "They produce a single partition that optimizes within-group homogeneity,
# instead of a hierarchical series of partitions optimizing the hierarchical
# attribution of objects to clusters"
# ___________________________________________________________________________
traits_ww <- readRDS(file.path(
  data_cache,
  "preproc_traits.rds"
))

# Hierarchical clustering
hc_output <- list()

set.seed(1234)
for (region in c("AUS", "EU", "NOA", "NZ", "SA")) {
  dat <- traits_ww[[region]]
  
  # Create dist. matrix using Overlap index Manly ---------------------
  # TODO (check what this approach does)
  # if scann == FALSE, than method 2 from dist.prop is used
  # (see function dist.ktab)
  # which is the Overlap index Manly
  dist_mat <- dist.ktab(dat, type = "F") # %>%
  #   as.matrix(.)
  
  # Contribution of grouping features to the global distance
  contrib_global_dist <- kdist.cor(dat, type = "F")
  
  # Clustering ---------------------------------------------------------
  # hc_taxa <- hclust(as.dist(dist_mat), method = "ward.D")
  hc_taxa <- hclust(dist_mat, method = "ward.D")
  
  # Get labels of dendrogram
  dend_label <- hc_taxa %>%
    as.dendrogram() %>%
    labels()
  
  # Optimal number of groups using the gap statistic
  # The gap statstic, works also well when data
  # fall into "one cluster" (i.e. indication that there is no cluster structure
  # if this is the case)
  gap <- clusGap(
    x = as.matrix(dist_mat),
    FUN = mycluster_hc,
    K.max = 15,
    B = 500
  )
  
  optimal_nog <- maxSE(gap$Tab[, "gap"],
                       gap$Tab[, "SE.sim"],
                       method = "Tibs2001SEmax")
  
  # Ordination & single linkage clustering --------------------------------
  pcoa <- dudi.pco(dist_mat, scannf = FALSE)
  hc_single <- hclust(dist_mat, method = "single")
  
  # Quality trait space
  q <- coranking(Xi = dist_mat, 
                 X = pcoa$li[, 1:2]) 
  nx <- coRanking::R_NX(q)
  qual_ordin <- coRanking::AUC_ln_K(nx)
  
  
  # Store results ------------------------------------------------------
  hc_output[[region]] <- list(
    "distance_matrix" = dist_mat,
    "contribu_global_dist" = contrib_global_dist,
    "hc_wardD" = hc_taxa,
    "labels_dendrogram" = dend_label,
    "gap_statistic" = gap,
    "optimal_nog" = optimal_nog,
    "pcoa" = pcoa,
    "quality_pcoa" = qual_ordin,
    "hc_single" = hc_single
  )
}
saveRDS(
  object = hc_output,
  file = file.path(
    data_cache,
    "hc_output_ww.rds"
  )
)

# ___________________________________________________________________________
# How would the results differ if another agglomeration method is used? -----
# So far Ward's agglomeration method has been used:
# "At each clustering step, Ward’s method finds the pair of objects
#  or clusters whose fusion increases as little as possible the sum,
# over all groups formed so far, of the squared distances between
# objects and cluster centroids; that sum is the total within-group
# sum-of-squares"(Legendre)
# other, see also Everitt and
# https://cran.r-project.org/web/packages/dendextend/vignettes/Cluster_Analysis.html#animals---attributes-of-animals

# ?Internal cluster validation
# to test for the presence of a
# hierarchical structure in the data
# Use Approach from Nemec -> SIGTREE
# Or any other measure? If we find an optimal number of groups
# doesn't this confirm our assumption on the structure of the data
# ___________________________________________________________________________
hclust_methods <-
  c(
    "ward.D",
    "single",
    "complete",
    "average",
    "mcquitty",
    "median",
    "centroid",
    "ward.D2"
  )

inv_traits_dendlist <- dendlist()
dend_pl <- list()
for (i in hclust_methods) {
  tmp_dend <- dist.ktab(traits_ww[["AUS"]], type = "F") %>%
    as.matrix(.) %>%
    as.dist(.) %>%
    hclust(., method = i) %>%
    as.dendrogram(.)

  dend_pl[[i]] <- tmp_dend  
  inv_traits_dendlist <- dendlist(inv_traits_dendlist, tmp_dend)
}
names(inv_traits_dendlist) <- hclust_methods

# Compare cophenetic correlation
# Cophenetic distance:
# "distance (or similarity) level at which objects x 1 and x 2 become
# members of the same cluster during the course of clustering"
cophenetic_cors <- cor.dendlist(inv_traits_dendlist)
corrplot::corrplot(cophenetic_cors, "pie", "lower")

# Compare agglomeration methods with optimal number of groups
# Fowlkes-Mallows index:
# 1 if two clusters conform
# 0 when they do not conform
cors_clustering_sol <- cor.dendlist(inv_traits_dendlist,
  method = "FM_index",
  k = hc_output$AUS$optimal_nog
)
corrplot::corrplot(cors_clustering_sol, "pie", "lower")