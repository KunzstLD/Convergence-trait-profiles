# ______________________________________________________________
# Data preparation 
# TODO: col_two_levels actually selects also for columns with less than
# two levels
# ______________________________________________________________

# read in data
trait_subsets <- readRDS(file.path(data_cache, "trait_subsets.rds"))

preproc_data <- list()
for (i in names(trait_subsets)) {
  data <- trait_subsets[[i]]
  
  # convert to data.frame -> data table does not support row.names
  setDF(data)
  
  # throw out all columns except traits and taxa col!
  rm <-
    grep("ER[0-9]{1,}|species|genus|family|order|tax.*",
         names(data))
  data <- data[, -rm]
  
  # data with only two levels/integer/factor variables have to be numeric
  col_int_fac <-
    names(Filter(function(y)
      is.integer(y) | is.factor(y), data))
  data[, col_int_fac] <- apply(data[, col_int_fac], 2, as.double)
  
  # add row.names
  row.names(data) <- data$ID_AQEM
  
  preproc_data[[i]] <- data
}
saveRDS(object = preproc_data,
        file = file.path(data_cache,
                         "preproc_data.rds"))
