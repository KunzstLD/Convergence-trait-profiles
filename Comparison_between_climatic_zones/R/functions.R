# _________________________________________________________________________
#### Statistical Analysis
# _________________________________________________________________________

# Helper functions data preparation ---------------------------------------

#### create individual pattern of trait name (not category!) ####
# i.e. feed_herbivore, feed_shredder -> feed
create_pattern_ind <- function(x, non_trait_cols) {
  if (missing(non_trait_cols)) {
    trait_names_pattern <- sub("\\_.*|\\..*", "", names(x)) %>%
      unique() %>%
      paste0("^", .)
  } else{
    pat <- paste0(non_trait_cols, collapse = "|")
    # get trait names & create pattern for subset
    trait_names_pattern <-
      grep(pat, names(x), value = TRUE, invert = TRUE) %>%
      sub("\\_.*|\\..*", "", .) %>%
      unique() %>%
      paste0("^", .)
  }
  trait_names_pattern
}

#### check for completeness of trait dataset ####
completeness_trait_data <- function(x, non_trait_cols) {
  trait_names_pattern <- create_pattern_ind(
    x = x,
    non_trait_cols = non_trait_cols
  )

  # test how complete trait sets are
  output <- matrix(ncol = 2, nrow = length(trait_names_pattern))
  for (i in seq_along(trait_names_pattern)) {
    # vector containing either 0 (no NAs) or a number (> 0) meaning that all
    # entries for this trait contained NA
    vec <-
      x[, apply(.SD, 1, function(y) {
        base::sum(is.na(y))
      }),
      .SDcols = names(x) %like% trait_names_pattern[[i]]
      ]

    # How complete is the dataset for each individual trait?
    output[i, ] <-
      c(
        (length(vec[vec == 0]) / nrow(x)) %>%
          `*`(100) %>%
          round(),
        trait_names_pattern[[i]]
      )
  }
  return(as.data.frame(output))
}


# Helper functions Clusering & RF -----------------------------------------

#### Helper F: Clustering ####
mycluster_hc <- function(x, k) {
  list(cluster = cutree(hclust(as.dist(x),
                               method = "ward.D2"),
                        k = k))
}

#### Helper F: RF Analysis ####
# custom prediction function
custom_pred <- function(object, newdata) {
  pred <- predict(object, newdata)$predictions
  avg <- purrr::map_df(as.data.frame(pred), mean)
  return(avg)
}

# Functions for testing RF ----------------------------------------------

#### Meta Clustering ####
# creates distance matrix
# calculates gap statistic and optimal number of groups
# creates dendrogram and provides a data.table with clusters and labels
meta_hclustering <- function(x) {
  dist_mat <- decostand(x,
    "norm",
    na.rm = TRUE
  ) %>%
    vegdist(., "euclidean", na.rm = TRUE) %>%
    as.matrix()

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

  # create dendrogram and get labels
  hc_taxa <- hclust(as.dist(dist_mat), method = "ward.D2")

  # get dendrogram
  dend_taxa <- as.dendrogram(hc_taxa)

  # grouping of taxa
  data_cluster <-
    data.table(
      family = names(cutree(dend_taxa,
        k = optimal_nog
      )),
      groups = cutree(dend_taxa, k = optimal_nog)
    )
  # output
  list(
    "data_cluster" = data_cluster,
    "hc_element" = hc_taxa,
    "optimal_nog" = optimal_nog
  )
}

#### Meta RF ####
meta_rf <- function(train,
                    test) {

  # Number of features
  n_features <- length(setdiff(names(train), "groups"))

  # Grid for different hyperparameters
  hyper_grid <- expand.grid(
    mtry = c(1, 5, 10, 15, 20, n_features - 1),
    node_size = seq(1, 10, by = 3),
    sample_size = c(.632, .80),
    OOB_error = NA,
    rmse = NA
  )

  # RF
  for (j in seq_len(nrow(hyper_grid))) {
    # train model
    model <- ranger(
      formula = groups ~ .,
      data = train,
      seed = 123,
      verbose = FALSE,
      num.trees = n_features * 10,
      mtry = hyper_grid$mtry[j],
      min.node.size = hyper_grid$node_size[j],
      sample.fraction = hyper_grid$sample_size[j]
    )

    # add OOB error to grid
    hyper_grid$OOB_error[j] <- model$prediction.error
  }

  # Feature importance
  # Use best tuning parameters
  best_set <- hyper_grid[order(hyper_grid$OOB_error), ][1, ]

  # Re-run model with impurity-based variable importance
  m3_ranger_impurity <- ranger(
    formula = groups ~ .,
    data = train,
    num.trees = 100,
    mtry = best_set$mtry,
    min.node.size = best_set$node_size,
    sample.fraction = best_set$sample_size,
    importance = "impurity",
    verbose = FALSE,
    seed = 123
  )

  # Re-run model with permutation-based variable importance
  m3_ranger_permutation <- ranger(
    formula = groups ~ .,
    data = train,
    num.trees = 100,
    mtry = best_set$mtry,
    min.node.size = best_set$node_size,
    sample.fraction = best_set$sample_size,
    importance = "permutation",
    verbose = FALSE,
    seed = 123
  )

  # Predictions
  # Training data
  res_train <- predict(m3_ranger_impurity, train)
  pred_train <- confusionMatrix(res_train$predictions, train$groups)

  # Test data
  res_test <- predict(m3_ranger_impurity, test)

  u <- union(res_test$predictions, test$groups)
  tab <- table(
    factor(res_test$predictions, u),
    factor(test$groups, u)
  )
  pred_test <- confusionMatrix(tab)

  list(
    "rf_impurity" = m3_ranger_impurity,
    "rf_permutation" = m3_ranger_permutation,
    "pred_train" = pred_train,
    "pred_test" = pred_test)
}

# _________________________________________________________________________

#### Trait aggregation ####

# Direct aggregation to specified taxonomic level
direct_agg <- function(trait_data,
                       non_trait_cols,
                       method,
                       taxon_lvl,
                       na.rm = TRUE) {
  # get names of trait columns
  pat <- paste0(non_trait_cols, collapse = "|")
  trait_col <- grep(pat, names(trait_data), value = TRUE, invert = TRUE)

  # aggregate to specified taxon lvl
  # Before applying this function, subset that no NA values occur in data
  # (otherwise all NA entries are viewed as a group & aggregated as well)
  agg_data <- trait_data[,
    lapply(.SD, method, na.rm = na.rm),
    .SDcols = trait_col,
    by = taxon_lvl
  ]
  agg_data
}

# _________________________________________________________________________
# TODO: Writing a custom distance function for fuzzy coded and binary data
# or, implement Orlocis Chord distance yourself to figure out which 
# variable contributed the most to the distances

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