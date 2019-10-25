# _________________________________________________________
#### Analysis without subsetting any dataset ####
# _________________________________________________________

# load data
file_list <- load_data(pattern = "*.rds", path = data_in)

# initialize list for clustered output data
data_cluster <- list()

# initialize list for predictions from RF on test data
pred_on_test_data <- list()

# initialize list to store most important variables inside
most_important_variables <- list()

for(i in seq_along(file_list)){
  # load data
  data <- file_list[[i]]
  
  # define name of loaded dataset
  name_dataset <- sub("(.+\\_)(.+)(\\_.+)", "\\2", names(file_list)[[i]])
  
  # data prep
  source(file.path(data_scr, "01_data_preparation.R"))
  
  # HC
  source(file.path(data_scr, "02_HC_analysis.R"))
  
  # RF 
  source(file.path(data_scr, "03_RF_variable_importance.R"))
  
}


# _________________________________________________________
#### Subsets ####
# TODO: For subset analysis remove trait that is subsetted?
# Remarks: dev_holometabol -> too few species in NoA dataset
# _________________________________________________________

# original aggregated data
file_list <- load_data(pattern = "*.rds", path = data_in)

# initialize list for clustered output data
data_cluster <- list()

# initialize list for predictions from RF on test data
pred_on_test_data <- list()

# initialize list to store most important variables inside
most_important_variables <- list()

# choose trait to subset
trait_subsetted <- "dev_holometabol"

for (i in seq_along(file_list)) {

  # for subsets
  data <- subset_trait_data(data = file_list[[i]],
                            trait = trait_subsetted,
                            trait_value = 1)
  
  source(file.path(data_scr, "01_data_preparation.R"))
  
  # name of loaded dataset
  name_dataset <-
    sub("(.+\\_)(.+)(\\_.+)", "\\2", names(file_list)[[i]]) %>%
    paste0("_", trait_subsetted)
  
  # data prep
  source(file.path(data_scr, "01_data_preparation.R"))
  
  # HC
  source(file.path(data_scr, "02_HC_analysis.R"))
  
  # RF
  source(file.path(data_scr, "03_RF_variable_importance.R"))
}





