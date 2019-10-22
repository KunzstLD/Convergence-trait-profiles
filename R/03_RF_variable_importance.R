# Multinominal classification random forest  -------------------------------------------------
# TODO: Update RF tuning
# Automation 
# Create a script for graphics

# merge cluster form HC to data
# data <- merge(data, data_cluster[[i]][, c("taxa", "groups_podani")],
#               by.x = "family", by.y = "taxa")
data <- merge(data, 
              data_cluster[, c("taxa", "groups_podani")],
              by.x = "family", 
              by.y = "taxa")

# rename response variable
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

# number of features
n_features <- length(setdiff(names(data), "Classf"))

# test & train data
# m_ranger <- ranger(
#   formula = Classf ~ .,
#   data = Traits_train,
#   num.trees = 500,
#   verbose = FALSE,
#   seed = 123
# )

# confusion matrix (error matrix)
# m_ranger$confusion.matrix

# ---------------------------------------------------------
#### Tuning ####
# ---------------------------------------------------------

# Grid for different hyperparameters
hyper_grid <- expand.grid(
  mtry        = c(1, 10, 15, 22, n_features -1),
  node_size   = seq(1, 10, by = 3),
  sample_size = c(.632, .80),
  OOB_error   = NA,
  rmse = NA
)

# test RF
for (j in 1:nrow(hyper_grid)) {
  # train model
  model <- ranger(
    formula         = Classf ~ .,
    data            = Traits_train,
    seed            = 123,
    verbose         = FALSE,
    num.trees       = n_features * 10,
    mtry            = hyper_grid$mtry[j],
    min.node.size   = hyper_grid$node_size[j],
    sample.fraction = hyper_grid$sample_size[j]
  )
  
  # add OOB error to grid
  hyper_grid$OOB_error[j] <- model$prediction.error
}

# top ten models accoridng to OOB error
# two configurations lead to lowest OOB error
hyper_grid[order(hyper_grid$OOB_error),] %>%
  head(., 10)

# default prediction error:
m_ranger$prediction.error
# Could also use the OOB RMSE

# For multiple columns use
# hyper_grid[with(hyper_grid, order(OOB_error)),]
# or the dplyr way
# hyper_grid %>%
# dplyr::arrange(OOB_error) %>%
# head(10)

# plot with respect to mtry ()
# plot(hyper_grid$mtry, hyper_grid$OOB_error)

# ---------------------------------------------------------
#### Feature importance
# ---------------------------------------------------------

# use best tuning parameters
best_set <- hyper_grid[order(hyper_grid$OOB_error), ][1, ]

# re-run model with impurity-based variable importance
m3_ranger_impurity <- ranger(
  formula         = Classf ~ .,
  data            = Traits_train,
  num.trees       = n_features * 10,
  mtry            = best_set$mtry,
  min.node.size   = best_set$node_size,
  sample.fraction = best_set$sample_size,
  importance      = 'impurity',
  verbose         = FALSE,
  seed            = 123
)

# re-run model with permutation-based variable importance
m3_ranger_permutation <- ranger(
  formula         = Classf ~ .,
  data            = Traits_train,
  num.trees       = n_features * 10,
  mtry            = best_set$mtry,
  min.node.size   = best_set$node_size,
  sample.fraction = best_set$sample_size,
  importance      = 'permutation',
  verbose         = FALSE,
  seed            = 123
)

# list object for saving results
# m3_ranger_permutation[[i]] <- ranger(
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

# Zero variables were not used for a split
# length(which(m3_ranger_impurity$variable.importance == 0))

# Zero variables whose values were randomized had no influence on the performance of the rf model
# length(which(m3_ranger_permutation$variable.importance == 0))

## most important variables according to impurtiy-based VI in distinguishing TPGs
# TODO: Automate
p1 <- vip(m3_ranger_impurity,
          num_features = n_features,
          bar = FALSE) +
  ggtitle(paste("Impurity-based variable importance", "EU")) +
  geom_point(size = 3) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(family = "Fira", size = 11),
    axis.text.y = element_text(family = "Fira", size = 11)
  )
 
## most important variables according to permutation-based VI in distinguishing TPGs
p2 <- vip(m3_ranger_permutation,
          num_features = 24,
          bar = FALSE) +
  ggtitle(paste("Permutation-based variable importance", "EU")) +
  geom_point(size = 3) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(family = "Fira", size = 11),
    axis.text.y = element_text(family = "Fira", size = 11)
  )
 
## Plot variable importance
# Most important variables in distinguishing TOG seems to be the pattern of development
# (and respiration tegument)
# Maybe this analysis just using species with the same pattern of development
# would be more meaningful?
png(
  file = file.path(
    data_out,
    "Graphs",
    paste0("Variable_importance_", "EU", ".png")
  ),
  width = 1500,
  height = 1300,
  res = 100
)
gridExtra::grid.arrange(p1, p2, nrow = 1)
dev.off()

# #### Feature effects ####
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

# ------------------------------------------------------------------
#### Predict for the test (unseen) data
# ------------------------------------------------------------------
pred_class <-
  predict(m3_ranger_impurity, Traits_test[, -grep("Classf", names(Traits_test))])

# Asses performance on test data
caret::confusionMatrix(factor(pred_class$predictions), factor(Traits_test$Classf))$overall[["Accuracy"]]

# u <- union(pred_class$predictions, Traits_test$Classf)
# t <-
#   table(factor(pred_class$predictions, u),
#         factor(Traits_test$Classf, u))
# pred_on_test_data[[i]] <-
#   caret::confusionMatrix(t)$overall[["Accuracy"]]