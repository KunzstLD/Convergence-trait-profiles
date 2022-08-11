# _________________________________________________________________________
#### Data handling ####
# _________________________________________________________________________

# loads data from input directory, stores them in a list
# and assigns a name (filename)
load_data <- function(path, pattern){
  files <- list.files(path = path, pattern = pattern)
  data <- lapply(files, function(y) readRDS(file = file.path(path, y)))
  data <- setNames(data, files)
  data
}

# Helper function for dcast calls
# Add to fun.aggregate argument for assigning a 1 to casted values when there length
# is greater than 0
fun_binary_length <- function(y) {
  as.numeric(ifelse(length(y) == 0, 0, 1))
}

# In case for subsetting (currently only for one variable can be subsetted)
subset_trait_data <- function(data, trait, trait_value){
  data  <- data[get(trait) == trait_value, ]
  data
}

# check if columnNR of elements in a list is the same
check_columnNr <- function(x){
  ck <- lapply(x, length) %>% unlist()
  if(abs(max(ck) - min(ck)) == 0){
    print("Files have the same number of columns")
  }
}

# check if colnames of data.frames stored in a list are the same 
# TODO: Improve output, but message when colnames differ
check_colNames <- function(x) {
  col_names_first_element <- lapply(x, names)[[1]]
  lapply(x, function(y) {
    all(names(y) %in% col_names_first_element)
  })
}

# create individual pattern of trait name (not category!)
# i.e. feed_herbivore, feed_shredder -> feed
create_pattern_ind <- function(x, non_trait_cols) {
  if (missing(non_trait_cols)) {
    trait_names_pattern <- sub("\\_.*|\\..*", "", names(x)) %>%
      unigooque() %>%
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

# check for completeness of trait datasets
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

# Normalization of trait scores 
# All trait states of one trait are divided by their row sum
# Hence, trait affinities are represented as "%" or ratios
normalize_by_rowSum <- function(x, 
                                non_trait_cols, 
                                na.rm = TRUE) {
  # get trait names & create pattern for subset
  trait_names_pattern <- create_pattern_ind(x = x,
                                            non_trait_cols = non_trait_cols)
  
  # loop for normalization (trait categories for each trait sum up to 1)
  for (cols in trait_names_pattern) {
    # get row sum for a specific trait
    x[, rowSum := apply(.SD, 1, sum, na.rm = na.rm),
      .SDcols = names(x) %like% cols]
    
    # get column names for assignment
    col_name <- names(x)[names(x) %like% cols]
    
    # divide values for each trait state by
    # the sum of trait state values
    x[, (col_name) := lapply(.SD, function(y) {
      round(y / rowSum, digits = 2)
    }),
    .SDcols = names(x) %like% cols]
  }
  # del rowSum column
  x[, rowSum := NULL]
  return(x)
}


# _________________________________________________________________________
#### Data extraction ####
# _________________________________________________________________________

# Helper FUN: most important traits to distinguish TPGs
extract_vi <- function(data) {
  lapply(data, function(y)
    y$rf_permutation$variable.importance) %>%
    lapply(., function(y)
      names(y[order(-y)])[1:5])
}

# Helper FUN: change in prediction error
extract_test_error <- function(data) {
  lapply(data, function(y)
    y$pred_test$overall[["Accuracy"]]) %>%
    unlist()
} 

# Helper FUN: change in training error
extract_train_error <- function(data) {
  lapply(data, function(y)
    y$pred_train$overall[["Accuracy"]]) %>%
    unlist()
}


# __________________________________________________________________________________________________
#### Statistical Analysis ####
# __________________________________________________________________________________________________

# Clustering --------------------------------------------------------------
mycluster_hc <- function(x, k) {
  list(cluster = cutree(hclust(as.dist(x),
                               method = "ward.D"),
                        k = k))
}

# RF Analysis ------------------------------------------------------------

# Custom prediction function
custom_pred <- function(object, newdata) {
  pred <- predict(object, newdata)$predictions
  avg <- purrr::map_df(as.data.frame(pred), mean)
  return(avg)
}

# Meta rf function with hyperparameter tuning 
meta_rf <- function(train,
                    test) {
  
  # Number of features
  n_features <- length(setdiff(names(train), "group"))
  
  # Grid for different hyperparameters
  hyper_grid <- expand.grid(
    mtry = c(1, 5, 10, 15, 20, n_features - 1),
    node_size = seq(1, 10, by = 3),
    sample_size = c(.632, .80),
    num_trees = 100,
    OOB_error = NA,
    rmse = NA
  )
  
  # RF
  for (j in seq_len(nrow(hyper_grid))) {
    # train model
    model <- ranger(
      formula = group ~ .,
      data = train,
      seed = 123,
      verbose = FALSE,
      num.trees = hyper_grid$num_trees[[j]],
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
    formula = group ~ .,
    data = train,
    num.trees = best_set$num_trees,
    mtry = best_set$mtry,
    min.node.size = best_set$node_size,
    sample.fraction = best_set$sample_size,
    importance = "impurity",
    verbose = FALSE,
    seed = 123
  )
  
  # Re-run model with permutation-based variable importance
  m3_ranger_permutation <- ranger(
    formula = group ~ .,
    data = train,
    num.trees = best_set$num_trees,
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
  pred_train <- confusionMatrix(res_train$predictions, train$group)
  
  # Test data
  res_test <- predict(m3_ranger_impurity, test)
  
  u <- union(res_test$predictions, test$group)
  tab <- table(
    factor(res_test$predictions, u),
    factor(test$group, u)
  )
  pred_test <- confusionMatrix(tab)
  
  list(
    "rf_impurity" = m3_ranger_impurity,
    "rf_permutation" = m3_ranger_permutation,
    "pred_train" = pred_train,
    "pred_test" = pred_test)
}

### Multivariate Welch-Test ####
# https://github.com/alekseyenko/WdStar

dist.ss2 = function(dm2, f){ #dm2 is matrix of square distances; f factor
  K = sapply(levels(f), function(lev) f==lev)
  t(K)%*%dm2%*%K/2
}


generic.distance.permutation.test <-
  function(test.statistic,
           dm,
           f,
           nrep = 999,
           strata = NULL) {
    N = length(f)
    generate.permutation = function() {
      f[sample(N)]
    }
    
    if (!is.null(strata)) {
      # map elements of each strata back to their positions in the factor variable
      strata.map = order(unlist(tapply(seq_along(f), strata, identity)))
      generate.permutation = function() {
        p = unlist(tapply(f, strata, sample)) # permute within strata
        p[strata.map]
      }
    }
    
    stats = c(test.statistic(dm, f),
              replicate(nrep,
                        test.statistic(dm, generate.permutation())))
    
    p.value = sum(stats >= stats[1]) / (nrep + 1)
    statistic = stats[1]
    list(p.value = p.value,
         statistic = statistic,
         nrep = nrep)
  }

WdS <- function(dm, f) {
  # This method computes Wd* statistic for distance matrix dm and factor f
  ns = table(f)
  SS2 = dist.ss2(as.matrix(dm) ^ 2, f)
  s2 = diag(SS2) / ns / (ns - 1)
  W = sum(ns / s2)
  
  idxs = apply(combn(levels(f), 2), 2, function(idx)
    levels(f) %in% idx)
  
  Ws = sum(apply(idxs, 2,
                 function(idx)
                   sum(ns[idx]) / prod(s2[idx]) *
                   (sum(SS2[idx, idx]) / sum(ns[idx]) - sum(diag(
                     SS2[idx, idx]
                   ) / ns[idx]))))
  k = nlevels(f)
  h = sum((1 - ns / s2 / W) ^ 2 / (ns - 1))
  Ws / W / (k - 1) / (1 + (2 * (k - 2) / (k ^ 2 - 1)) * h)
}

WdS.test <- function(dm, f, nrep = 999, strata = NULL) {
  generic.distance.permutation.test(
    WdS,
    dm = dm,
    f = f,
    nrep = nrep,
    strata = strata
  )
}


# _________________________________________________________________________
#### Plotting ####
# _________________________________________________________________________

# Dendrogram plot ---------------------------------------------------------
fun_dendrog_pl <- function(hc,
                           optimal_nog,
                           labels,
                           hang_height = 0.001) {
  hc %>% 
    as.dendrogram() %>%
    color_branches(k = optimal_nog) %>%
    hang.dendrogram(hang_height = hang_height) %>%
    set("labels_cex", 0.7) %>%
    dendextend::ladderize() %>%
    set("labels", labels) 
}

# Heatmap plot -----------------------------------------------------------
# Create heatmap with ggplot for TPGs and Grouping features for each continent
# Columns of data are hardcoded, maybe change in the future
fun_heatmap_single_cont <- function(data) {
  ggplot(data, aes(x = family,
                   y = trait,
                   fill = affinity)) +
    geom_tile() +
    facet_grid(
      factor(grouping_feature) ~ group,
      scales = "free",
      space = "free",
      labeller = as_labeller(grouping_feature_names)
    ) +
    scale_fill_gradient(name = "Trait affinity",
                        low = "#FFFFF1",
                        high = "#012345") +
    labs(x = "Family", y = "Trait") +
    theme_bw() +
    theme(
      axis.title = element_text(size = 22),
      axis.text.x = element_text(
        family = "Roboto Mono",
        size = 16,
        angle = 90,
        hjust = 1,
        vjust = 0.2
      ),
      axis.text.y = element_text(family = "Roboto Mono",
                                 size = 16),
      legend.title = element_text(family = "Roboto Mono",
                                  size = 16),
      legend.text = element_text(family = "Roboto Mono",
                                 size = 16),
      strip.text = element_text(family = "Roboto Mono",
                                size = 16),
      plot.title = element_text(family = "Roboto Mono", 
                                size = 20),
      panel.grid = element_blank()
    )
}

# Plotting function for tpg that occur across the tested continents
fun_heatmap_tpg <- function(data,
                            facet_names) {
  ggplot(data, aes(x = family,
                   y = trait,
                   fill = affinity)) +
    geom_tile() +
    facet_grid(
      . ~ continent,
      scales = "free",
      space = "free",
      labeller = as_labeller(facet_names)
    ) +
    scale_fill_gradient(name = "Trait affinity",
                        low = "#FFFFF1",
                        high = "#012345") +
    labs(x = "Family", y = "Trait") +
    theme_bw() +
    theme(
      axis.title = element_text(size = 20),
      axis.text.x = element_text(
        family = "Roboto Mono",
        size = 16,
        angle = 90,
        hjust = 1,
        vjust = 0.2
      ),
      axis.text.y = element_text(family = "Roboto Mono",
                                 size = 16),
      legend.title = element_text(family = "Roboto Mono",
                                  size = 16),
      legend.text = element_text(family = "Roboto Mono",
                                 size = 16),
      strip.text = element_text(family = "Roboto Mono",
                                size = 15),
      plot.title = element_text(family = "Roboto Mono",
                                size = 16),
      panel.grid = element_blank()
    )
}

# Helper function to plot boruta results --------------------------------------

# Does require data.table
# Takes as input a result from attStats() on an object created by boruta()
fun_boruta_results <- function(data) {
  data$traits <- rownames(data)
  setDT(data)
  data[, traits := reorder(traits, meanImp)]
  data[, decision := factor(decision, levels = c("Confirmed", "Tentative", "Rejected"))]
  
  ggplot(data,
         aes(
           x = as.factor(traits),
           y = meanImp,
           color = decision
         )) +
    geom_point(size = 2.5) +
    scale_color_brewer(palette = "Dark2") + 
    labs(x = "Traits", y = "Mean permutation importance score", 
         color = "Decision") +
    theme_bw() +
    theme(
      axis.title = element_text(size = 16),
      axis.text.x = element_text(
        family = "Roboto Mono",
        size = 14,
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(family = "Roboto Mono",
                                 size = 14),
      legend.title = element_text(family = "Roboto Mono",
                                  size = 16),
      legend.text = element_text(family = "Roboto Mono",
                                 size = 14),
      strip.text = element_text(family = "Roboto Mono",
                                size = 14),
      plot.title = element_text(family = "Roboto Mono", 
                                size = 16),
      panel.grid = element_blank()
    )
}

# __________________________________________________________________________________________________
#### Trait aggregation ####
# __________________________________________________________________________________________________

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
