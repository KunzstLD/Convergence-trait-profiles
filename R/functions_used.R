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


# _________________________________________________________________________
#### Statistical Analysis ####
# _________________________________________________________________________

# Clustering --------------------------------------------------------------
mycluster_hc <- function(x, k) {
  list(cluster = cutree(hclust(as.dist(x),
                               method = "ward.D2"),
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
  ggplot(data, aes(x = trait_label,
                   y = family,
                   fill = affinity)) +
    geom_tile() +
    facet_grid(
      group ~ factor(grouping_feature),
      scales = "free",
      space = "free",
      labeller = as_labeller(grouping_feature_names)
    ) +
    scale_fill_gradient(name = "Trait affinity",
                        low = "#FFFFF1",
                        high = "#012345") +
    labs(x = "Trait", y = "Family") +
    theme_bw() +
    theme(
      axis.title = element_text(size = 12),
      axis.text.x = element_text(
        family = "Roboto Mono",
        size = 11,
        angle = 60,
        hjust = 1
      ),
      axis.text.y = element_text(family = "Roboto Mono",
                                 size = 11),
      legend.title = element_text(family = "Roboto Mono",
                                  size = 12),
      legend.text = element_text(family = "Roboto Mono",
                                 size = 11),
      strip.text = element_text(family = "Roboto Mono",
                                size = 11),
      panel.grid = element_blank()
    )
}

# Plotting function for tpg that occur across the tested continents
fun_heatmap_tpg <- function(data) {
  ggplot(data, aes(x = trait,
                   y = family,
                   fill = affinity)) +
    geom_tile() +
    facet_grid(continent ~ .,
               scales = "free",
               space = "free") +
    scale_fill_gradient(name = "Trait affinity",
                        low = "#FFFFF1",
                        high = "#012345") +
    labs(x = "Trait", y = "Family") +
    theme_bw() +
    theme(
      axis.title = element_text(size = 12),
      axis.text.x = element_text(family = "Roboto Mono",
                                 size = 11),
      axis.text.y = element_text(family = "Roboto Mono",
                                 size = 11),
      legend.title = element_text(family = "Roboto Mono",
                                  size = 12),
      legend.text = element_text(family = "Roboto Mono",
                                 size = 11),
      strip.text = element_text(family = "Roboto Mono",
                                size = 11),
      panel.grid = element_blank()
    )
}

