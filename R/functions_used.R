# _________________________________________________________________________
#### Data handling ####
# _________________________________________________________________________

# loads RDS data from input directory, stores them in a list
# and assigns a name (filename)
load_data <- function(pattern, path){
  files <- list.files(path = path, pattern = pattern)
  data <- lapply(files, function(y) readRDS(file = file.path(data_in, y)))
  data <- setNames(data, files)
  data
}

# In case for subsetting (currently only for one variable can be subsetted)
subset_trait_data <- function(data, trait, trait_value){
  data  <- data[get(trait) == trait_value, ]
  data
}

# check if columnNR of elements in a list is the same
check_columnNr <- function(x){
  ck <- lapply(x, length) %>% unlist()
  if(abs(max(ck) - min(ck)) == 0){
    print("Files have the same number of columns")
  }
}

# check if colnames of data.frames stored in a list are the same 
# TODO: Improve output, but message when colnames differ
check_colNames <- function(x) {
  col_names_first_element <- lapply(x, names)[[1]]
  lapply(x, function(y) {
    all(names(y) %in% col_names_first_element)
  })
}

# _________________________________________________________________________
#### Statistical Analysis
# _________________________________________________________________________

# Clustering --------------------------------------------------------------
mycluster_hc <- function(x, k) {
  list(cluster = cutree(hclust(as.dist(x),
                               method = "ward.D2"),
                        k = k))
}

# RF Analysis -------------------------------------------------------------
# custom prediction function
custom_pred <- function(object, newdata) {
  pred <- predict(object, newdata)$predictions
  avg <- purrr::map_df(as.data.frame(pred), mean)
  return(avg)
}

