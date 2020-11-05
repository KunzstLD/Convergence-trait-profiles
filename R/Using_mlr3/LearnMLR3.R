##################################################
# https://mlr3book.mlr-org.com/index.html
##################################################


##### Using trait data ####

# libraries
library(mlr3)
library(mlr3learners)

# data
data <- merge(data,
  data_cluster[[dataset]][, c("taxa", "groups")],
  by.x = "family",
  by.y = "taxa"
)

saveRDS(data, file = file.path(getwd(), "Output", "mlr_ex_data.rds"))
data <- readRDS(file = file.path(getwd(), "Output", "mlr_ex_data.rds"))

#### preprocessing ####
# remove taxonomical information
data <- data[, -grep("order|family", names(data))]

# response variable has to be encoded as factor or character
data$groups <- as.factor(data$groups)

# create task:
task_data <- TaskClassif$new(id = "trait_data", 
							backend = data, 
							target = "groups")

#### Modelling ####

# create learner 
learner <- lrn("classif.ranger", importance = "permutation")
train_set <- sample(task_data$nrow, 0.7 * task_data$nrow)
test_set <- setdiff(seq_len(task_data$nrow), train_set)

# View parameters: 
learner$param_set

# TODO: Works, but where are the hyperparameters
# Why no split into train and test set? 
# Check not automated version!

# Hyperparameter tuning
library(paradox)
library(mlr3tuning)
tune_ps <- ParamSet$new(list(
  ParamInt$new("mtry", lower = 1, upper = task_data$ncol-1),
  ParamInt$new("min.node.size", lower = 1, upper = 10)
))
terminator <- trm("evals", n_evals = 10)
tuner <- tnr("random_search")


# AutoTuner inherits from Learner base class
at <- AutoTuner$new(
  learner = learner,
  resampling = rsmp("holdout"),
  measure = msr("classif.ce"),
  search_space = tune_ps,
  terminator = terminator,
  tuner = tuner
)
at

grid <- benchmark_grid(
  task = task_data,
  learner = list(at, lrn("classif.rpart")),
  resampling = rsmp("cv", folds = 3)
)

# avoid console output from mlr3tuning
logger <- lgr::get_logger("bbotk")
logger$set_threshold("warn")

bmr <- benchmark(grid)
bmr$aggregate(msrs(c("classif.ce", "time_train")))

# training with best hyperparamters
learner$train(task, row_ids = train_set)
print(learner$model)

# extract importance
learner$importance()

# predicting
prediction <- learner$predict(task, row_ids = test_set)
print(prediction)
as.data.table(prediction)

# confusion matrix
prediction$confusion


#######################################################################