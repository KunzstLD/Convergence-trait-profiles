# Multinominal classification random forest  -------------------------------------------------
# TODO: Update RF tuning
# Create a script for graphics
# Function script

# merge cluster form HC to data
data <- merge(data, data_cluster[[i]][, c("taxa", "groups_podani")],
              by.x = "family", by.y = "taxa")
names(data)[names(data) %in% "groups_podani"] <- "Classf"
# response variable has to be encoded as factor or character
data$Classf <- as.factor(data$Classf)

# remove taxonomical information
data <- data[, -grep("order|family", names(data))]

# create training and test data set
set.seed(123)
split <- initial_split(data[, -1], prop = .63, strata = "Classf")
Traits_train <- training(split)
Traits_test  <- testing(split)

# test & train data
m_ranger <- ranger(
  formula = Classf ~ .,
  data = Traits_train,
  num.trees = 500,
  verbose = FALSE,
  seed = 123
)

# confusion matrix (error matrix)
# m_ranger$confusion.matrix

#### Tuning ####
# Grid for different parameters
hyper_grid <- expand.grid(
  num_trees   = c(250, 500),
  # number of variables to possibly split at
  mtry        = seq(1, 24, by = 3),
  node_size   = seq(1, 10, by = 3),
  sample_size = c(.632, .80),
  splitrule   = c("gini", "extratrees"),
  OOB_error   = 0
)


# test RF
# start.time <- Sys.time()
for (j in 1:nrow(hyper_grid)) {
  # train model
  model <- ranger(
    formula         = Classf ~ .,
    data            = Traits_train,
    seed            = 123,
    verbose         = FALSE,
    num.trees       = hyper_grid$num_trees[j],
    mtry            = hyper_grid$mtry[j],
    min.node.size   = hyper_grid$node_size[j],
    sample.fraction = hyper_grid$sample_size[j],
    splitrule       = hyper_grid$splitrule[j]
  )
  
  # add OOB error to grid
  hyper_grid$OOB_error[j] <- model$prediction.error
}
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken

# plot with respect to mtry ()
# plot(hyper_grid$mtry, hyper_grid$OOB_error)

# two configurations lead to lowest OOB error
# hyper_grid[order(hyper_grid$OOB_error),]
# multiple columns:
# hyper_grid[with(hyper_grid, order(OOB_error)),]
# dplyr way
# hyper_grid %>%
# dplyr::arrange(OOB_error) %>%
# head(10)
# Tune further!


#### Feature importance
# use best tuning parameters
best_set <- hyper_grid[order(hyper_grid$OOB_error), ][1, ]

# re-run model with impurity-based variable importance
m3_ranger_impurity <- ranger(
  formula         = Classf ~ .,
  data            = Traits_train,
  num.trees       = best_set$num_trees,
  mtry            = best_set$mtry,
  min.node.size   = best_set$node_size,
  sample.fraction = best_set$sample_size,
  splitrule       = best_set$splitrule,
  importance      = 'impurity',
  verbose         = FALSE,
  seed            = 123
)

# re-run model with permutation-based variable importance
# m3_ranger_permutation <- ranger(
#   formula         = Classf ~ .,
#   data            = Traits_train,
#   num.trees       = best_set$num_trees,
#   mtry            = best_set$mtry,
#   min.node.size   = best_set$node_size,
#   sample.fraction = best_set$sample_size,
#   splitrule       = best_set$splitrule,
#   importance      = 'permutation',
#   verbose         = FALSE,
#   seed            = 123
# )

# list object for saving results
m3_ranger_permutation[[i]] <- ranger(
  formula         = Classf ~ .,
  data            = Traits_train,
  num.trees       = best_set$num_trees,
  mtry            = best_set$mtry,
  min.node.size   = best_set$node_size,
  sample.fraction = best_set$sample_size,
  splitrule       = best_set$splitrule,
  importance      = 'permutation',
  verbose         = FALSE,
  seed            = 123
)

# how many variables not used for a split
# length(which(m3_ranger_impurity$variable.importance == 0))
## [1] 0

# how many variables that whose values were randomized had no influence
# on the performance of the rf model
# length(which(m3_ranger_permutation$variable.importance == 0))
## [1] 0

# ## most important variables according to impurtiy-based VI in distinguishing TPGs
# p1 <- vip(m3_ranger_impurity,
#           num_features = 24,
#           bar = FALSE) +
#   ggtitle(paste("Impurity-based variable importance", name_dataset)) +
#   geom_point(size = 3) +
#   theme_bw() +
#   theme(
#     axis.title = element_text(size = 15),
#     axis.text.x = element_text(family = "Fira", size = 11),
#     axis.text.y = element_text(family = "Fira", size = 11)
#   )
# 
# ## most important variables according to permutation-based VI in distinguishing TPGs
# p2 <- vip(m3_ranger_permutation,
#           num_features = 24,
#           bar = FALSE) +
#   ggtitle(paste("Permutation-based variable importance", name_dataset)) +
#   geom_point(size = 3) +
#   theme_bw() +
#   theme(
#     axis.title = element_text(size = 15),
#     axis.text.x = element_text(family = "Fira", size = 11),
#     axis.text.y = element_text(family = "Fira", size = 11)
#   )
# 
# ## plot
# png(
#   file = file.path(
#     data_out,
#     "Graphs",
#     paste0("Variable_importance_", name_dataset, ".png")
#   ),
#   width = 1500,
#   height = 1300,
#   res = 100
# )
# gridExtra::grid.arrange(p1, p2, nrow = 1)
# dev.off()
# 
# #### Feature effects
# # # Select certain features and inspect on which category this variable
# # # has most (or least) influence
# # # Use of Partial dependence plots (PDP) & ICE curves
# # # -> probabilty model is needed here (gives class probabilities instead of class assignment)
# #
# # # probability model
# m3_ranger_prob <- ranger(
#   formula         = Classf ~ .,
#   data            = Traits_train,
#   num.trees       = best_set$num_trees,
#   mtry            = best_set$mtry,
#   min.node.size   = best_set$node_size,
#   sample.fraction = best_set$sample_size,
#   splitrule       = best_set$splitrule,
#   importance      = 'permutation',
#   probability     = TRUE,
#   verbose         = FALSE,
#   seed            = 123
# )
# #
# # # fetch most important variable (according to permutation based importance)
# most_import_var <-
#    m3_ranger_permutation$variable.importance[order(-m3_ranger_permutation$variable.importance)][1] %>% names
# #
# # # partial dependence
# pd <- pdp::partial(
#   m3_ranger_prob,
#   pred.var = most_import_var,
#   pred.fun = custom_pred,
#   train = Traits_train
# )
# #
# # # For classification where the machine learning model outputs probabilities,
# # # the partial dependence plot displays the probability for a certain class given
# # # different values for feature(s).
# # # An easy way to deal with multiple classes is to draw one line or plot per class
# pdp1 <-
#   ggplot(pd, aes(pd[, most_import_var], yhat, color = factor(yhat.id))) +
#   geom_point(show.legend = FALSE, size = 3) +
#   facet_wrap(~ yhat.id, nrow = 2) +
#   labs(x = most_import_var) +
#   theme_bw() +
#   theme(
#     axis.title = element_text(size = 15),
#     axis.text.x = element_text(family = "Fira", size = 11),
#     axis.text.y = element_text(family = "Fira", size = 11)
#   ) +
#   expand_limits(y = 0)
# 
# png(
#   file = file.path(
#     data_out,
#     "Graphs",
#     paste0("Partial_dependence_", name_dataset, ".png")
#   ),
#   width = 1100,
#   height = 1300,
#   res = 100
# )
# pdp1
# dev.off()
# #
# # #### Predict for the test (unseen) data
# pred_class <-
#   predict(m3_ranger_impurity, Traits_test[, -grep("Classf", names(Traits_test))])
# 
# ## Class predictions for the remaining taxa
# # pred_class$predictions
# 
# ## Asses performance on test data
# u <- union(pred_class$predictions, Traits_test$Classf)
# t <-
#   table(factor(pred_class$predictions, u),
#         factor(Traits_test$Classf, u))
# pred_on_test_data[[i]] <-
#   caret::confusionMatrix(t)$overall[["Accuracy"]]

# caret::confusionMatrix(factor(pred_class$predictions), factor(Traits_test$Classf))$overall[["Accuracy"]]
