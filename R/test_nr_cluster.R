# Test Optimal number of clusters stability
traits_ww <- readRDS(file.path(
  data_cache,
  "preproc_traits.rds"
))

# Hierarchical clustering
list_nog <- list()

for (iter in 1:1000) {
    for (region in c("NOA")) {
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
            method = "Tibs2001SEmax"
        )
    }
    list_nog[[iter]] <- optimal_nog
}
sum(unlist(list_nog) == 8)
sum(unlist(list_nog) == 11)
