# ______________________________________________________________
# Data preparation 
# TODO: col_two_levels actually selects also for columns with less than
# two levels
# ______________________________________________________________

# read in data
trait_subsets <- readRDS(file.path(data_cache, "trait_subsets.rds"))

# lists to keep the data
preproc_data <- list()
preproc_data_omit <- list()

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
  data$ID_AQEM <- NULL
  
  # try two approaches: 
  # 1) remove for every dataset incomplete observations!
  # That means only taxa with complete trait profiles remain
  data_omit <- na.omit(data)
  # 2) Substitute NA's with zeros
  # Need to adjust the Similarity function!  
  data[is.na(data)] <- 0
  
  preproc_data_omit[[i]] <- data_omit
  preproc_data[[i]] <- data
}
saveRDS(object = preproc_data_omit,
        file = file.path(data_cache,
                         "preproc_data_omit.rds"))
saveRDS(object = preproc_data,
        file = file.path(data_cache,
                         "preproc_data.rds"))