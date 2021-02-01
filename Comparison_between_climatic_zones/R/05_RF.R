# ________________________________________________________
# Multinominal classification random forest
# TODO: Update RF tuning
# Create a script for graphics
# ________________________________________________________

# load data
Hmisc::describe(data_cluster$groups)

# response variable has to be encoded as factor or character
data_cluster[, groups := as.factor(groups)]

# separate taxonomical information
data_cluster <- data_cluster[, -c(
    "ID_AQEM",
    "species",
    "genus",
    "family",
    "order"
)]

# create subset for testing
cl_subset <- data_cluster[, .SD,
    .SDcols = names(data_cluster) %like% "groups|feed.*|locom.*"
]

# create training and test data set
set.seed(123)

ind <- createDataPartition(cl_subset$groups, p = 0.7)
train <- cl_subset[ind[[1]], ]
test <- cl_subset[-ind[[1]], ]

# ---------------------------------------------------------
#### Tuning ####
# ---------------------------------------------------------

# number of features
n_features <- length(setdiff(names(cl_subset), "groups"))

# Grid for different hyperparameters
hyper_grid <- expand.grid(
    mtry = c(1, 5, n_features - 1),
    node_size = seq(1, 10, by = 3),
    sample_size = c(.632, .80),
    OOB_error = NA,
    rmse = NA
)

# test RF
for (j in 1:nrow(hyper_grid)) {
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

# top ten models accoridng to OOB error
# two configurations lead to lowest OOB error
hyper_grid[order(hyper_grid$OOB_error), ] %>%
    head(., 10)

# plot with respect to mtry ()
plot(hyper_grid$mtry, hyper_grid$OOB_error)

# ---------------------------------------------------------
#### Feature importance
# ---------------------------------------------------------

# use best tuning parameters
best_set <- hyper_grid[order(hyper_grid$OOB_error), ][1, ]

# re-run model with impurity-based variable importance
m3_ranger_impurity <- ranger(
    formula = groups ~ .,
    data = train,
    num.trees = n_features * 10,
    mtry = best_set$mtry,
    min.node.size = best_set$node_size,
    sample.fraction = best_set$sample_size,
    importance = "impurity",
    verbose = FALSE,
    seed = 123
)

# re-run model with permutation-based variable importance
m3_ranger_permutation <- ranger(
    formula = groups ~ .,
    data = train,
    num.trees = n_features * 10,
    mtry = best_set$mtry,
    min.node.size = best_set$node_size,
    sample.fraction = best_set$sample_size,
    importance = "permutation",
    verbose = FALSE,
    seed = 123
)

# Variables were not used for a split
length(which(m3_ranger_impurity$variable.importance == 0))
# Variables whose values were randomized 
# had no influence on the performance of the rf model
length(which(m3_ranger_permutation$variable.importance == 0))

# most important variables according to impurtiy-based VI in distinguishing TPGs
# TODO: Automate
impo_imp <- vip(m3_ranger_impurity,
    num_features = 5,
    bar = FALSE
) +
    ggtitle(paste("Impurity-based variable importance temperate EU")) +
    geom_point(size = 3) +
    theme_bw() +
    theme(
        axis.title = element_text(size = 15),
        axis.text.x = element_text(family = "Fira", size = 11),
        axis.text.y = element_text(family = "Fira", size = 11)
    )

# most important variables according to permutation-based VI in distinguishing TPGs
impo_perm <- vip(m3_ranger_permutation,
    num_features = 5,
    bar = FALSE
) +
    ggtitle("Permutation-based variable importance temperate EU") +
    geom_point(size = 3) +
    theme_bw() +
    theme(
        axis.title = element_text(size = 15),
        axis.text.x = element_text(family = "Fira", size = 11),
        axis.text.y = element_text(family = "Fira", size = 11)
    )

# Plot variable importance
# png(
#   file = file.path(
#     data_out,
#     "Graphs",
#     paste0("Variable_importance_", dataset, ".png")
#   ),
#   width = 1500,
#   height = 1300,
#   res = 100
# )
# gridExtra::grid.arrange(impo_imp,
#                         impo_perm,
#                         nrow = 1)
# dev.off()
#
# # get most important variables
# most_important_variables[[dataset]] <-
#   data.frame(
#     importance_impurity = impo_imp$data$Variable[1:10],
#     importance_permutation = impo_perm$data$Variable[1:10],
#     region = dataset
#   )

# ------------------------------------------------------------------
#### Predict for the test (unseen) data
# ------------------------------------------------------------------
pred_class <-
    predict(m3_ranger_impurity, test[, -c("groups")])

# Asses performance on test data
caret::confusionMatrix(factor(pred_class$predictions),
                       factor(test$Classf))$overall[["Accuracy"]]

u <- union(pred_class$predictions, test$groups)
t <- table(
    factor(pred_class$predictions, u),
    factor(test$groups, u)
)
pred_on_test_data <- caret::confusionMatrix(t)$overall[["Accuracy"]]