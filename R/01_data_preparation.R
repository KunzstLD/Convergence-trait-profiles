# ______________________________________________________________
# Data preparation 
# TODO: col_two_levels actually selects also for columns with less than
# two levels
# ______________________________________________________________

# read in data
# data <- readRDS(file.path(data_in, "Trait_AUS_agg.rds"))

# analyse base way -> data table does not support row.names
setDF(data)

# turn cols into ordered factors
# trait_cols <- names(data)[!names(data) %like% "order|family"]
# data[, trait_cols] <- lapply(data[, trait_cols], as.ordered)

# data with only two levels have to be numeric 
col_two_levels <- data[, -grep("order|family", names(data))] %>%
  lapply(levels) %>%
  lapply(., function(y)
    length(y)[length(y) <= 2]) %>%
  unlist %>%
  names

if(length(col_two_levels) == 1) {
  data[[col_two_levels]] <- sapply(data[[col_two_levels]], as.numeric)
}
if (length(col_two_levels) > 1) {
  data[, col_two_levels] <- lapply(data[, col_two_levels], as.numeric)
}

# add row.names
row.names(data) <- data$family