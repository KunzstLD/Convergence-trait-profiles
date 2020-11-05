# _________________________________________________________
#### Analysis without subsetting any dataset ####

# TODO:
# Questions to adress:
# - Which method for HC/ complete linkage?
# - Create look-up table with traits
# - Update Aggregated datasets

# _________________________________________________________

# load data
file_list <- load_data(pattern = "*.rds", path = data_in)

# checks if they have the same colnames
file_list %>% check_colNames()

# rm dev
file_list <- lapply(
  file_list,
  function(y) y[, c("dev_holometabol", "dev_hemimetabol") := NULL]
)

# rename file list
names(file_list) <- sub(
  "([A-z]{5})(\\_)([A-z]{2,})(\\_)([A-z]{3})(.+)",
  "\\3", 
  names(file_list)
)
names(file_list)[3] <- "NOA"

# initialize list for clustered output data
data_cluster <- list()

# initialize list for predictions from RF on test data
pred_on_test_data <- list()

# initialize list to store most important variables inside
most_important_variables <- list()

# get perm_importance plots
impo_perm <- list()

for (dataset in names(file_list)) {

  # load data
  data <- file_list[[dataset]]

  # data prep
  source(file.path(data_scr, "01_data_preparation.R"))

  # HC
  source(file.path(data_scr, "02_HC_analysis.R"))

  # RF
  source(file.path(data_scr, "03_RF_variable_importance.R"))
}

gridExtra::grid.arrange(impo_perm$AUS,
  impo_perm$EU,
  impo_perm$NOA,
  impo_perm$NZ,
  nrow = 2,
  ncol = 2
)
ggsave(
  filename = file.path(data_out, "most_imp_traits_summary.png"),
  width = 20,
  height = 10.5,
  unit = "cm"
)

# _________________________________________________________
#### Subsets ####
# TODO: For subset analysis remove trait that is subsetted?
# Remarks: dev_holometabol -> too few species in NoA dataset
# _________________________________________________________

# # original aggregated data
# file_list <- load_data(pattern = "*.rds", path = data_in)
#
# # initialize list for clustered output data
# data_cluster <- list()
#
# # initialize list for predictions from RF on test data
# pred_on_test_data <- list()
#
# # initialize list to store most important variables inside
# most_important_variables <- list()
#
# # choose trait to subset
# trait_subsetted <- "dev_holometabol"
#
# for (i in seq_along(file_list)) {
#
#   # for subsets
#   data <- subset_trait_data(data = file_list[[i]],
#                             trait = trait_subsetted,
#                             trait_value = 1)
#
#   source(file.path(data_scr, "01_data_preparation.R"))
#
#   # name of loaded dataset
#   name_dataset <-
#     sub("(.+\\_)(.+)(\\_.+)", "\\2", names(file_list)[[i]]) %>%
#     paste0("_", trait_subsetted)
#
#   # data prep
#   source(file.path(data_scr, "01_data_preparation.R"))
#
#   # HC
#   source(file.path(data_scr, "02_HC_analysis.R"))
#
#   # RF
#   source(file.path(data_scr, "03_RF_variable_importance.R"))
# }
