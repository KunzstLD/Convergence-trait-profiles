# _____________________________________________________________
# Testing if RF's VI is susceptible to collinearity
# I would expect that the most important variables
# are changing with each RF run due to collinearity
# _____________________________________________________________

# Load setup script
source("./Comparison_between_climatic_zones/R/Set_up.R")

# Load Trait_EU_agg.rds or any similar dataset
trait_eu <- readRDS("./Data/Trait_EU_agg.rds")

# Rm dev variables and order
trait_eu[, c("dev_hemimetabol", "dev_holometabol") := NULL]

# Save taxonomical information
trait_eu_cp <- copy(trait_eu)

# Rm order column
trait_eu[, order := NULL]

# Create a data.frame and assign families as row.names
setDF(trait_eu)
row.names(trait_eu) <- trait_eu$family
trait_eu$family <- NULL

# Source analysis script
source(file.path(data_scr, "test_RF_02.R"))

# Remove highest ranked variables (VI) and run again
vi_perm_full_run <- output_rf$rf_permutation$variable.importance
most_imp_var <- names(vi_perm_full_run[order(-vi_perm_full_run)][1])

# store results of predictions


# ___________________________________________________________________________
# - Repeat the procedure by removing most important variable
#   -> Does the prediction error change?
#   -> Which variables are now most important?
# - ? Repeat with randomly introduced variables?
trait_eu <- trait_eu[, names(trait_eu)[names(trait_eu) != most_imp_var]]
trait_eu_cp <- trait_eu_cp[, -..most_imp_var]
source(file.path(data_scr, "test_RF_02.R"))

# Some changes between results of full sized dataset compared to
# when removing the most important variable identified in the first
# analysis
vi_perm_lo_size <- output_rf$rf_permutation$variable.importance

pred_train_lo_size <- output_rf$pred_train$overall[["Accuracy"]]
pred_test_lo_size <- output_rf$pred_test$overall[["Accuracy"]]


# ____________________________________________________________________________
# - Repeat by rm whole grouping feature size
trait_eu <- trait_eu[, names(trait_eu)[!names(trait_eu) %like% "size"]]
trait_eu_cp <- trait_eu_cp[, .SD, .SDcols = !names(trait_eu_cp) %like% "size"]
source(file.path(data_scr, "test_RF_02.R"))

vi_perm_lo_gf <- output_rf$rf_permutation$variable.importance

pred_train_lo_gf <- output_rf$pred_train$overall[["Accuracy"]]
pred_test_lo_gf <- output_rf$pred_test$overall[["Accuracy"]]

# Summarizing results
data.table(
    VI_full = names(vi_perm_full_run[order(-vi_perm_full_run)]),
    VI_full_val = vi_perm_full_run[order(-vi_perm_full_run)],
    VI_lo_size = c(names(vi_perm_lo_size[order(-vi_perm_lo_size)]), NA),
    VI_lo_size_val = c(vi_perm_lo_size[order(-vi_perm_lo_size)], NA),
    VI_lo_size_compl = c(
        names(vi_perm_lo_gf[order(-vi_perm_lo_gf)]),
        NA, NA, NA
    ),
    VI_lo_size_compl_val = c(
        vi_perm_lo_gf[order(-vi_perm_lo_gf)],
        NA, NA, NA
    )
)
