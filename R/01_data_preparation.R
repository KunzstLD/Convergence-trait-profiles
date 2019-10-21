# -------------------------------------------------------------
# Data preparation 
# -------------------------------------------------------------

# read in data
data <- readRDS(file.path(".", "Data", "Trait_EU_agg.rds"))

# change NAs to zero (NAs are true zeros)
for (j in names(data)) {
  data.table::set(data, which(is.na(data[[j]])), j, 0)
}

# analyse base way -> data table does not support row.names
setDF(data)

# turn cols into ordered factors
trait_cols <- names(data)[!names(data) %like% "order|family"]
data[, cols] <- lapply(data[, cols], as.ordered)

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