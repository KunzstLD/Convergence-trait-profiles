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
      y / rowSum
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
#### Clustering ####
# __________________________________________________________________________________________________
mycluster_hc <- function(x, k) {
  list(cluster = cutree(hclust(as.dist(x),
                               method = "ward.D2"),
                        k = k))
}

# __________________________________________________________________________________________________
#### RF Analysis ####
# __________________________________________________________________________________________________

# Custom prediction function
custom_pred <- function(object, newdata) {
  pred <- predict(object, newdata)$predictions
  avg <- purrr::map_df(as.data.frame(pred), mean)
  return(avg)
}

# Meta rf function with hyperparameter tuning 
# old, now with mlr3
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

# __________________________________________________________________________________________________
### Multivariate Welch-Test ####
# https://github.com/alekseyenko/WdStar
# __________________________________________________________________________________________________
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

##### Dendrogram plot ####
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

##### Heatmap plot ####
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

##### Helper function to plot boruta results ####

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

# __________________________________________________________________________________________________
#### Null Models ####
# __________________________________________________________________________________________________

# Calculate pcoa for simulated datasets
# id column should be created before
calc_pcoa <-
  function(x,
           traits = c("feed", "resp", "volt", "locom", "size", "bf")) {
    trait_patterns <- paste0(traits, collapse = ".*|")
    x <- x[, .SD, .SDcols = patterns(paste0(trait_patterns, "|id"))]
    setDF(x)
    
    taxa_reg_names <- x$id
    row.names(x) <- taxa_reg_names
    x$id <- NULL
    
    # Konvert to ktab object
    vec <- sub("\\_.*", "\\1", names(x))
    blocks <- rle(vec)$lengths
    x <- prep.fuzzy(x, blocks)
    x <- ktab.list.df(list(x))
    dist <- dist.ktab(x, type = "F")
    
    # PcOA
    pcoa <- dudi.pco(dist, scannf = FALSE, nf = 25)
  }

# Transform pcoa element obtained from "calc_pcoa" into a data.table
transf_pcoa_dt <- function(pcoa_obj){
  pcoa_scores <- pcoa_obj$li[1:2]
  pcoa_scores$id <- rownames(pcoa_scores)
  setDT(pcoa_scores)
}

# Calculate convex hull for first two PCoA axes
# and subsequently overlap between the continents/regions
calc_cnx_hull <- function(scores) {
  hull <- scores %>%
    group_by(continent) %>%
    slice(chull(A1, A2))
  setDT(hull)
  hull_split <- split(hull[, .(A1, A2)], f = hull$continent)
  list(
    "Overlap" = Overlap(hull_split),
    "Overlap_symmetric" = Overlap(hull_split, symmetric = TRUE)
  )
}

# Calculate ellipses for first two PCoA axes
# Maximum likelihood approach
calc_ellipses <- function(scores,
                          ellipses = c("1.AUS",
                                       "1.EU",
                                       "1.NOA",
                                       "1.NZ",
                                       "1.SA")) {
  # we just have one community
  pcoa_siber <- createSiberObject(scores[, .(
    iso1 = A1,
    iso2 = A2,
    group = continent,
    community = 1
  )])
  
  # Calculate all possible permutations of ellipses
  perm_ellipses <- gtools::permutations(n = 5,
                                        r = 2,
                                        v = ellipses)
  perm_ellipses <- as.data.frame(perm_ellipses)
  rownames(perm_ellipses) <- paste0(perm_ellipses$V1,
                                    "_",
                                    perm_ellipses$V2)
  
  # Calc overlap
  ellipse95_overlap <- list()
  for (i in 1:nrow(perm_ellipses)) {
    ellipse95_overlap[[i]] <- maxLikOverlap(
      perm_ellipses[i, "V1"],
      perm_ellipses[i, "V2"],
      pcoa_siber,
      p.interval = 0.95,
      n = 100
    )
  }
  names(ellipse95_overlap) <- rownames(perm_ellipses)
  ellipse95_overlap
}


# Hierarchical clustering & optimal number of groups
# method is ward.D2
# uses funciton mycluster_hc
calc_clustering <-
  function(x,
           traits = c("feed", "resp", "volt", "locom", "size", "bf")) {
    trait_patterns <- paste0(traits, collapse = ".*|")
    x <- x[, .SD, .SDcols = patterns(paste0(trait_patterns, "|id"))]
    setDF(x)
    
    taxa_reg_names <- x$id
    row.names(x) <- taxa_reg_names
    x$id <- NULL
    
    # Konvert to ktab object
    vec <- sub("\\_.*", "\\1", names(x))
    blocks <- rle(vec)$lengths
    x <- prep.fuzzy(x, blocks)
    x <- ktab.list.df(list(x))
    dist <- dist.ktab(x, type = "F")
    
    # HC & optimal number of clusters
    hc <- hclust(dist, method = "ward.D2")
    dend <- as.dendrogram(hc)
    
    gap <- clusGap(
      x = as.matrix(dist),
      FUN = mycluster_hc,
      K.max = 15,
      B = 500
    )
    
    optimal_nog <- maxSE(gap$Tab[, "gap"],
                         gap$Tab[, "SE.sim"],
                         method = "Tibs2001SEmax")
    list("hc" = hc,
         "dend" = dend,
         "optimal_nog" = optimal_nog)
  }

# Add TPGs to simulated datasets
add_tpgs_td <- function(cl_obj) {
  results <- list()
  for (i in names(cl_obj)) {
    results[[i]] <-
      data.table(
        family = names(
          cutree(
            cl_obj[[i]]$dend,
            k = cl_obj[[i]]$optimal_nog,
            order_clusters_as_data = FALSE
          )
        ),
        group = cutree(
          cl_obj[[i]]$dend,
          k = cl_obj[[i]]$optimal_nog,
          order_clusters_as_data = FALSE
        )
      )
  }
  results <- rbindlist(results, idcol = "continent")
}

# Calculate dbrda for mean trait profiles
calc_dbrda <-
  function(t_data,
           traits = c("feed", "resp", "volt", "locom", "size", "bf")) {
    trait_patterns <- paste0(traits, collapse = ".*|")
    x <- t_data[, .SD, .SDcols = patterns(paste0(trait_patterns, "|id"))]
    setDF(x)
    
    taxa_reg_names <- x$id
    row.names(x) <- taxa_reg_names
    x$id <- NULL
    
    # Convert to ktab object
    vec <- sub("\\_.*", "\\1", names(x))
    blocks <- rle(vec)$lengths
    x <- prep.fuzzy(x, blocks)
    x <- ktab.list.df(list(x))
    dist <- dist.ktab(x, type = "F")
    
    # dbrda with continents
    dbrda_res <- dbrda(formula = dist ~ continent, data = t_data)
    dbrda_res <- summary(dbrda_res)
    expl_var <- dbrda_res$constr.chi / dbrda_res$tot.chi
  }

# RF with tpgs as dependent variable
calc_rf_tpgs <- function(x,
                         traits = c("feed", "resp", "volt", "locom", "size", "bf")) {
  most_imp_vars <- list()
  # ls_instance <- list()
  # scores_test <- list()
  # scores_train <- list()
  trait_pat <- paste0(traits, collapse = ".*|")
  for (cont in unique(x$continent)) {
    set.seed(1234)
    dat <- x[continent == cont, ]
    dat <-
      dat[, .SD, .SDcols = patterns(paste0(trait_pat, "|group"))]
    
    # Split in train and test data (stratified)
    ind <- createDataPartition(dat$group, p = 0.7)
    train <- dat[ind[[1]],]
    test <- dat[-ind[[1]],]
    
    # Create tasks
    task_train <- TaskClassif$new(
      id = paste0(cont, "_train"),
      backend = train,
      target = "group"
    )
    
    # Specify stratification for CV
    task_train$col_roles$stratum <- "group"
    
    # Create random forest learner
    # list of possible learners: mlr3learners
    rf_learner <- lrn("classif.ranger",
                      predict_type = "prob",
                      importance = "permutation")
    
    # Set up search space
    # rf_learner$param_set
    search_space <- ps(
      mtry = p_int(lower = 1, upper = length(task_train$data()) - 1),
      min.node.size = p_int(lower = 1, upper = 10)#,
      #  sample.fraction = p_dbl(lower = 0.632, upper = 0.8)
    )
    
    # Resampling
    resampling <- rsmp("cv",
                       folds = 5L)
    
    # Performance measure for resampling (multivariate Brier score)
    mbrier_metric <- mlr3::msr("classif.mbrier")
    
    # When to terminate tuning
    evals <- trm("evals", n_evals = 100)
    instance = TuningInstanceSingleCrit$new(
      task = task_train,
      learner = rf_learner,
      resampling = resampling,
      measure = mbrier_metric,
      search_space = search_space,
      terminator = evals
    )
    
    # Optimization via grid search
    tuner <- tnr("grid_search", resolution = 10)
    tuner$optimize(instance)
    # ls_instance[[cont]] <- instance
    
    # Train rf on train dataset with optimized parameters
    rf_learner$param_set$values <-
      instance$result_learner_param_vals
    rf_learner$train(task_train)
    pred_train <- rf_learner$predict_newdata(newdata = train)
    # scores_train[[cont]] <- pred_train$score(mlr3::msr("classif.mbrier"))
    
    # Check model on test dataset
    pred_test <- rf_learner$predict_newdata(newdata = test)
    # scores_test[[cont]] <- pred_test$score(mlr3::msr("classif.mbrier"))
    
    # Retrieve most important variables
    most_imp_vars[[cont]] <- rf_learner$importance()
  }
  list("most_imp_vars" = most_imp_vars)
} 



