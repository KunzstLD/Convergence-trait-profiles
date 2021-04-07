#___________________________________________________________________________________________________
#### Variable importance for the traits to determine which traits drive TPG selection ####
# 1) RF normal mode with permutation importance
# From the project description: 
# -Convergence of trait profiles: 
#   - Do we find the same two traits occur among the three most import traits?
#   - Repeat the procedure by removing most important variable
#   -> Does the prediction error change?
#   -> Which variables are now most important?
# Consider collinearity!
# Finding interactions!
# Second analyses using iRF!
#___________________________________________________________________________________________________

# Load trait data with grouping assignment from the cluster analysis
trait_AUS <- readRDS(file = file.path(data_cache, "trait_AUS_with_groups.rds"))
trait_EU <- readRDS(file = file.path(data_cache, "trait_EU_with_groups.rds"))
trait_NOA <- readRDS(file = file.path(data_cache, "trait_NOA_with_groups.rds"))
trait_NZ <- readRDS(file = file.path(data_cache, "trait_NZ_with_groups.rds"))

trait_dat <- list(
  "AUS" = trait_AUS,
  "EU" = trait_EU,
  "NOA" = trait_NOA,
  "NZ" = trait_NZ
)

# Check distribution within groups
obs_group <- lapply(trait_dat, function(y) y[, .N, by = group])
obs_group <- rbindlist(obs_group, idcol = "continent")  
obs_group[, group := factor(group, levels = c(1:10))]
ggplot(obs_group, aes(x = group, y = N)) +
  geom_pointrange(aes(ymin = 0, ymax = N)) +
  facet_wrap(as.factor(continent) ~., 
             scales = "free") +
  ylim(0, 20)


# "group" should be a factor variable 
trait_dat <- lapply(trait_dat, function(y) y[, group := factor(make.names(group))])

# --- Calculate RF ---------------------------------------------------------------------------------
results_rf <- list()
for(region in c("AUS", "EU", "NOA", "NZ")) {
  
  dat <- copy(trait_dat[[region]])
  dat <- dat[, -c("family", "order")]

  # Create training and test data set
  ind <- createDataPartition(dat$group, p = 0.75)
  train <- dat[ind[[1]], ]
  test <- dat[-ind[[1]], ]
  
  # Calculate RF results
  set.seed(1234)
  results_rf[[region]] <- meta_rf(train = train,
                             test = test)
}

# prepare results
rf_summary <- data.frame(extract_vi(results_rf))
setDT(rf_summary)
rf_summary <- melt(
  rf_summary,
  measure.vars = c("AUS", "EU", "NOA", "NZ"),
  variable.name = "continent",
  value.name = "most_important_traits"
)

acc_results <- data.frame(extract_test_error(results_rf),
                          extract_train_error(results_rf))
setDT(acc_results, keep.rownames = "continent")

rf_summary[acc_results,
           `:=`(test_accuracy = i.extract_test_error.results_rf.,
                train_accuracy = i.extract_train_error.results_rf.),
           on = "continent"]



# use caret
for(region in c("AUS", "EU", "NOA", "NZ")) {
  
  dat <- trait_dat[["AUS"]]
  dat[, c("family", "order") := NULL]
  n_features <- length(dat)
  
  ind <- createDataPartition(dat$group, p = 2/3)
  train <- dat[ind[[1]],]
  test <- dat[-ind[[1]],]
  
  # createFolds() does stratified k-fold cv
  # From the docs: "The random sampling is done within the levels of y (=outcomes) 
  # when y is a factor in an attempt to balance the class distributions within the splits."
  strat_folds <- createFolds(y = train$group, 
              k = 5)
  
  # Specify type of resampling, 5 -fold CV
  fitControl <- trainControl(method = "cv",
                             number = 5,
                             classProbs = TRUE,
                             index = strat_folds)
  
  # Tuning grid
  grid <- expand.grid(
    mtry = seq(2, n_features-1, 2),
    splitrule = c("gini"),
    min.node.size = c(1:10)
  )
  
  set.seed(1234)
  rf_fit <- train(
    group ~ .,
    data = train,
    method = "ranger",
    metric = "Kappa",
    trControl = fitControl,
    tuneGrid = grid,
    num.trees = (n_features - 1) * 10,
    verbose = FALSE,
    importance = "permutation"
   )
  
}

# --- Brier score ---------------------------------
# Do not use Accuracy but rather Brier score 
# to evaluate prediction strength of RF Model

# TODO Custom metric for CV
# function(data, 
#          lev = NULL,
#          model = NULL){
#   
# }

# On OOB data
rf_fit$finalModel

# On test dataset
pred_test <- predict(rf_fit, newdata = test, type = "prob")
setDT(pred_test)

pred_results_test <- data.table(
  ground_truth = as.numeric(sub("X", "", test$group)),
  pred_group = apply(pred_test, MARGIN = 1, function(y)
    which(y == max(y)))
)
pred_results_test[, o_t := ifelse(ground_truth == pred_group, 1, 0)]                              
pred_results_test[, pred_group_label := paste0("X", pred_group)]
pred_results_test[, f_t := apply(pred_test, MARGIN = 1, function(y)
  y[which(y == max(y))])]

pred_results_test[, mean((f_t - o_t)^2)]


# --- Variable selection with wrapper algorithm Boruta ---------------------------------------------
res <- list()
for (region in c("AUS", "EU", "NOA", "NZ")) {
  # Preprocessing
  data <- copy(trait_dat[[region]])
  data <- data[, -c("family", "order")]
  data[, group := factor(as.character(group))]
  
  res[[region]] <- Boruta(group ~ .,
                          data = data,
                          maxRuns = 100)
}
# TODO: show values of all iterations
boruta_res <- lapply(res, attStats)

for(region in names(boruta_res)) {
  plot <- fun_boruta_results(data = boruta_res[[region]]) +
    ggtitle(paste0(region, ": Variable selection with Boruta"))
  
  ggplot2::ggsave(
    filename = file.path(
      data_paper,
      "Graphs",
      paste0("Variable_selec_boruta_", region, ".png")
    ),
    width = 40,
    height = 20,
    units = "cm"
  )
}
