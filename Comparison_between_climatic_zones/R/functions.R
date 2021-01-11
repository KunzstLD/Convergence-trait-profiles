# _________________________________________________________________________
#### Statistical Analysis
# _________________________________________________________________________

# Clustering --------------------------------------------------------------
mycluster_hc <- function(x, k) {
  list(cluster = cutree(hclust(as.dist(x),
                               method = "ward.D2"),
                        k = k))
}

# RF Analysis -------------------------------------------------------------
# custom prediction function
custom_pred <- function(object, newdata) {
  pred <- predict(object, newdata)$predictions
  avg <- purrr::map_df(as.data.frame(pred), mean)
  return(avg)
}



# _________________________________________________________________________

# TODO: Writing a custom distance function for fuzzy coded and binary data
# or, implement Orlocis Chord distance yourself to figure out which variable contributed
# the most to the distances

# bin_cols <- names(Filter(function(y) is.integer(y) | is.factor(y), data))
# fc_cols <- names(data)[!names(data) %in% bin_cols]

# custom_dist <- function(x,
#                         fc.cols,
#                         binary.cols,
#                         asym.bin) {
#   # for fuzzy coded data: Orloci's chord distance
#   # na.rm deletes missing observation pairwise
  
#   dist_fc <-
#     vegdist(decostand(x[, fc.cols], "norm", na.rm = TRUE),
#             "euclidean", na.rm = TRUE)
#   # divide result by sqrt(2) to standardize to range [0, 1]
  
  
#   # binary data: Gower distance with asymm variables
#   dist_bin <- cluster::daisy(x[, binary.cols],
#                              type = list(asymm = asym.bin))
  
#   list("dist_mat_fc" = dist_fc,
#        "dist_mat_bin" = dist_bin)
# }

# test <- custom_dist(
#   x = data,
#   fc.cols = fc_cols,
#   binary.cols = bin_cols,
#   asym.bin = bin_cols
# )