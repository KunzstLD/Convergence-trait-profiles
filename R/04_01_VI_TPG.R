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

# Check distribution within groups (familie richness)
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
# lapply(trait_dat, dim)

# --- Calculate RF ---------------------------------------------------------------------------------

most_imp_vars <- list()
ls_instance <- list()
scores_test <- list()
scores_train <- list()

trait_dat_cp <- copy(trait_dat)
trait_dat_cp <- lapply(trait_dat_cp, function(y) y[, c("family", "order") := NULL])

for(region in names(trait_dat_cp)) {
  
  set.seed(1234)
  dat <- trait_dat_cp[[region]]
  
  # TODO Split in train and test data (stratified)
  ind <- createDataPartition(dat$group, p = 0.7)
  train <- dat[ind[[1]],]
  test <- dat[-ind[[1]],]
  
  # Create tasks
  task_train <- TaskClassif$new(id = paste0(region, "_train"),
                                backend = train,
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
    mtry = p_int(lower = 1, upper = length(task_train$data()) - 1),
    min.node.size = p_int(lower = 1, upper = 10)#,
    #  sample.fraction = p_dbl(lower = 0.632, upper = 0.8)
  )
  
  # Resampling
  resampling <- rsmp("cv",
                     folds = 5L)
  
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
  
  # Performance measure for resampling (multivariate Brier score)
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
  
  # Optimization via grid search
  tuner <- tnr("grid_search", resolution = 10)
  tuner$optimize(instance)
  ls_instance[[region]] <- instance

  # Train rf on train dataset with optimized parameters
  rf_learner$param_set$values <- instance$result_learner_param_vals
  rf_learner$train(task_train)
  pred_train <- rf_learner$predict_newdata(newdata = train)
  scores_train[[region]] <- pred_train$score(msr("classif.mbrier"))
  
  # Check model on test dataset
  pred_test <- rf_learner$predict_newdata(newdata = test)
  scores_test[[region]] <- pred_test$score(msr("classif.mbrier"))
  
  # Retrieve most important variables
  most_imp_vars[[region]] <- rf_learner$importance()
}

# Get optimal param. values and associated mbrier score
# ls_instance$AUS$archive$data[order(classif.mbrier),][1,]
# ls_instance$EU$archive$data[order(classif.mbrier),][1,]
# ls_instance$NOA$archive$data[order(classif.mbrier),][1,]
# ls_instance$NZ$archive$data[order(classif.mbrier),][1,]

# Inspect prediction re mbrier score
scores_test
scores_train

# save
scores_test <- do.call(rbind, scores_test) %>%
  as.data.frame()
names(scores_test)[names(scores_test) == "classif.mbrier"] <- "mbrier_score_test"
scores_test$continent <- rownames(scores_test)
scores_train <- do.call(rbind, scores_train) %>%
  as.data.frame()
names(scores_train)[names(scores_train) == "classif.mbrier"] <- "mbrier_score_train"
perf_summary <- cbind(scores_test, scores_train)
setDT(perf_summary)
setcolorder(perf_summary, 
            neworder = "continent")
saveRDS(perf_summary, file.path(data_cache, "perf_summary.rds"))

# Most important variables
most_imp_vars <- lapply(most_imp_vars, function(y) as.data.frame(y)) %>% 
  do.call(rbind, .)
saveRDS(most_imp_vars, file = file.path(data_cache, "most_imp_vars.rds"))


# --- Variable selection with wrapper algorithm Boruta ---------------------------------------------
output <- list()
for (region in c("AUS", "EU", "NOA", "NZ")) {
  
  set.seed(1234)
  dat <- trait_dat_cp[[region]]
  output[[region]] <- Boruta(group ~ .,
                          data = dat,
                          maxRuns = 100)
}

# TODO: show values of all iterations
boruta_res <- lapply(output, attStats)
saveRDS(boruta_res, 
        file = file.path(data_cache, "boruta_res.rds"))

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

# ?trait interactions

# --- Brier score ----------------------------------------------------------------------------------
# Do not use Accuracy but rather Brier score 
# to evaluate prediction strength of RF Model

# o_t <- list()
# idx <- as.numeric(sub("X", "",test$group))
# o_t <- matrix(data = 0,
#               ncol = length(unique(test$group)),
#               nrow = nrow(test))
# # Create o_t matrix based on the ground truth of the test set
# for(i in 1:nrow(o_t)){
#   o_t[i, idx[i]] <- 1
# }
# 
# # calculate Brier score
# res <- vector(mode = "numeric")
# for(i in 1:nrow(pred_test)) {
#   res[[i]] <- sum((pred_test[i,] - o_t[i,]) ^ 2)
# }
# mean(res)

# How to implement in caret? Or use other framework?
# https://stackoverflow.com/questions/66235293/how-to-obtain-brier-score-in-random-forest-in-r
# https://cran.r-project.org/web/packages/randomForestSRC/randomForestSRC.pdf
