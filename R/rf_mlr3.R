library(mlr3)
library(mlr3viz)
library(mlr3learners)
library(mlr3measures)
library(mlr3tuning)

# Load trait datasets and combine
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

# Example AUS dataset
dat <- trait_dat[["AUS"]]
dat[, c("family", "order") := NULL]

# TODO Split in train and test data (stratified)
ind <- createDataPartition(dat$group, p = 2/3)
train <- dat[ind[[1]], ]
test <- dat[-ind[[1]], ]

# Create tasks
task_train <- TaskClassif$new(id = "AUS_train",
                              backend = train,
                              target = "group")
task_test <- TaskClassif$new(id = "AUS_test",
                             backend = test,
                             target = "group")

# Specify stratification for cv 
task_train$col_roles$stratum <- "group"

# Create random forest learner
# list of possible learners: mlr3learners
rf_learner <- lrn("classif.ranger",
                  predict_type = "prob",
                  importance = "permutation")

# Set up search space
# rf_learner$param_set
search_space <- ps(
  mtry = p_int(lower = 1, upper = length(task_train$data())-1), 
  min.node.size = p_int(lower = 1, upper = 10)#,
#  sample.fraction = p_dbl(lower = 0.632, upper = 0.8)
)

# Resampling
resampling <- rsmp("cv", 
                   folds = 5L)
# resampling$instantiate(task)

# Check if stratification worked
# dt <- merge(resampling$instance, task$data()[, row_id := .I], by = "row_id")
# # overall class distribution
# dt[, class_distrib := (.N / nrow(dt)) * 100,
#    by = group]
# dt[, class_distrib := round(class_distrib, 2)]
# # class distribution in folds
# dt[, n_fold := .N, by = fold] 
# dt[, class_distrib_fold := (.N/n_fold)*100, by = .(fold, group)]
# dt[, class_distrib_fold := round(class_distrib_fold, 2)]
# unique(dt[, .(group, fold, class_distrib_fold, class_distrib)] %>% 
#          .[order(group, fold), ])

# Performance measure for resampling
mbrier_metric <- msr("classif.mbrier")

# When to terminate tuning
evals20 <- trm("evals", n_evals = 100)
instance = TuningInstanceSingleCrit$new(
  task = task_train,
  learner = rf_learner,
  resampling = resampling,
  measure = mbrier_metric,
  search_space = search_space,
  terminator = evals20
)
instance

# Optimization via grid search
tuner <- tnr("grid_search", resolution = 10)
set.seed(1234)
tuner$optimize(instance)

# Inspect optimal values
instance$result_learner_param_vals
instance$archive$data[order(mtry, min.node.size), ]
instance$archive$data[order(classif.mbrier), ]

# Train on full dataset with optimized parameters
rf_learner$param_set$values = instance$result_learner_param_vals
rf_learner$train(task_train)

# Check model on test dataset
predict_test = rf_learner$predict(task_test)
predict_test$score(msr("classif.mbrier"))

# Retrieve most important variables 
most_imp_vars <- rf_learner$importance()

# Run again with Boruta for quality check with same hyperparam values

# trait interactions



