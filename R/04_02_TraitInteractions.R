# __________________________________________________________________________________________________
# Searching for trait interactions
# using iRF 
# https://github.com/sumbose/iRF
# https://www.stat.berkeley.edu/~kkumbier/vignette.html
# __________________________________________________________________________________________________

library(iRF)
set.seed(1234)

trait_data <- copy(trait_AUS)
trait_data <- trait_data[, -c("family", "order")]
trait_data[, group := factor(as.character(group))]

# Create training and test data set
ind <- createDataPartition(trait_data$group, p = 0.75)
train <- trait_data[ind[[1]],]
test <- trait_data[-ind[[1]],]

p <- length(trait_data)-1
sel.prob <- rep(1/p, p)

# iteratively grow RF, use Gini importance of features as weights
rf <- list()
for (iter in 1:4) {
  rf[[iter]] <- randomForest(
    x = train[, 1:25],
    y = train$group,
    xtest = test[, 1:25],
    ytest = test$group,
    mtry.select.prob = sel.prob,
    keep.forest = TRUE
  )
  
  # update selection probabilities for next iteration
  sel.prob <- rf[[iter]]$importance / sum(rf[[iter]]$importance)
}


# Feature usage on nodes
rforest <- readForest(rand.forest = rf[[3]], x = as.matrix(train[, 1:25]))
head(rforest$tree.info, n=10)

# sparse matrix
rforest$node.feature[1:10, c(1:5, p + 1:5)]

# Finding prevalent sets of features and their high-order combinations used to define these nodes
# by random intersection trees (RIT)
# The following with the following command runs RIT on all class 1 nodes, sampling each leaf node 
# with probability proportional to the number of observations in each leaf node:
class1.nodes <- rforest$tree.info$prediction - 1 == 1
wt <- rforest$tree.info$size.node[class1.nodes]
RIT(rforest$node.feature[class1.nodes,], 
    weights=wt,
    depth=5, 
    branch=2,
    n_trees=100)

# Selecting stable interactions (using bootstrapping)
fit <- iRF(x = as.matrix(train[, 1:25]),
           y = train$group,
           xtest = as.matrix(test[, 1:25]),
           ytest = test$group, 
           n.iter=5, 
       #    n.core=n.cores,
           select.iter = TRUE,
           n.bootstrap=10
)
interaction_res <- fit$interaction
interaction_res[order(-prevalence), ]
# Read for interpretation: https://arxiv.org/abs/1810.07287.

