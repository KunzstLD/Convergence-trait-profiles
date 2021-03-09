#__________________________________________________________________________________________________
# Variable importance for the traits 
# 1) RF normal mode with permutation importance
# From the project description: 
# -Convergence of trait profiles: 
#   - Do we find the same two traits occur among the three most import traits?
#   - Repeat the procedure by removing most important variable
#   -> Does the prediction error change?
#   -> Which variables are now most important?
# - ? Repeat with randomly introduced variables?
# Consider collinearity!
# Finding interactions!
# Second analyses using iRF!
#__________________________________________________________________________________________________

# Load data
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
lapply(trait_dat, function(y) Hmisc::describe(y$group))

# RF with removal of variable with the highest importance value
# to see if the correlation structure influences the importance 
# ranking of the variables
var_rm <- NULL
results_rf <- list()
output <- list()

# TODO improve code (two loops is relatively slow!)
# TODO: Either use RFE or Boruta?
for(region in c("AUS", "EU", "NOA", "NZ")) {
  for (i in 1:5) {
    # Preprocessing
    data <- copy(trait_dat[[region]])
    data <- data[, -c("family", "order")]
    data[, group := factor(as.character(group))]
    
    # Rm most important variable
    # First run no variable is selected
    if (is.null(var_rm)) {
      data <- data
    } else{
      data <- data[, -..var_rm]
    }
    
    # Create training and test data set
    ind <- createDataPartition(data$group, p = 2 / 3)
    train <- data[ind[[1]],]
    test <- data[-ind[[1]],]
    
    # Calculate RF results
    set.seed(1234)
    results_rf[[i]] <- meta_rf(train = train,
                               test = test)
    
    # Select most important variable
    vi_perm <- results_rf[[i]]$rf_permutation$variable.importance
    var_rm <- names(vi_perm[order(-vi_perm)][1])
  }
  output[[region]] <- results_rf
}
c(impt_aus, impt_eu, impt_noa, impt_nz) %<-% output

# Combining the results (5 most important traits; prediction error of test and training data)
rf_vimp <- data.table(
  "five_most_important_traits" = c(
    extract_vi(impt_aus),
    extract_vi(impt_eu),
    extract_vi(impt_noa),
    extract_vi(impt_nz)
  ),
  "test_error" = c(
    extract_test_error(impt_aus),
    extract_test_error(impt_eu),
    extract_test_error(impt_noa),
    extract_test_error(impt_nz)
  ),
  "train_error" = c(
    extract_train_error(impt_aus),
    extract_train_error(impt_eu),
    extract_train_error(impt_noa),
    extract_train_error(impt_nz)
  ),
  "run" = rep(c(
    "initial", rep("rm_highest_ranking_var", times = 4)
  ),
  times = 4),
  "database" = rep(c("AUS", "EU", "NOA", "NZ"), each = 5)
  
)

kable(rf_vimp[, .(database, five_most_important_traits, run)], 
      format = "markdown",
      caption = "Five most important traits to distinguish TPGs for each continent. After the 
      initial run the highest ranking variable has been removed and a random forest has been fitted 
      again. This procedure has been repeated 4 times.")