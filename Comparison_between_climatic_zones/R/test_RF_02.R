# _____________________________________________________________
#### Pipeline ####
# Runs:
#  - Hierarchical clustering
#  - RF
# _____________________________________________________________

# Clustering --------------------------------------------------
output_clustering <- meta_hclustering(x = trait_eu)

# RF ----------------------------------------------------------

# Preprocessing
data <- merge(
    x = trait_eu_cp,
    y = output_clustering$data_cluster,
    by = "family"
)
data <- data[, -c("family", "order")]

# groups as factor
data[, groups := as.factor(as.character(groups))]

# check distribution within groups
Hmisc::describe(data$groups)

# create training and test data set
set.seed(123)

ind <- createDataPartition(data$groups, p = 0.7)
train <- data[ind[[1]], ]
test <- data[-ind[[1]], ]

output_rf <- meta_rf(
    train = train,
    test = test
)