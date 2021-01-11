# ______________________________________________________________
# Data preparation
# TODO: col_two_levels actually selects also for columns with less than
# two levels
# ______________________________________________________________

# read in data
trait_subsets <- readRDS(file.path(data_cache, "trait_subsets.rds"))
trait_subsets_uq <-
  readRDS(file.path(data_cache, "trait_subsets_uq.rds"))

# choose which data to process
# for interactive use:
if (interactive()) {
  omitted <-
    readline("Use full data? \n Otherwise only taxa specific to one climate region are used \n")

  if (omitted == "T" | omitted == "TRUE" | omitted == "True") {
    data_type <- "full"
  } else {
    data_type <- "unique"
  }
  rm(omitted)
}
# else choose data_type manually
trait_data <- switch(data_type,
  "full" = trait_subsets,
  "unique" = trait_subsets_uq
)

# data preprocessing
preproc_data <- list()

for (region in c("arid", "temperate", "cold", "polar")) {
  data <- trait_data[climateregion == region, ]

  # convert to data.frame -> data table does not support row.names
  setDF(data)

  # throw out all columns except traits and taxa col!
  rm <-
    grep(
      "ER[0-9]{1,}|species|genus|family|order|tax.*|climateregion",
      names(data)
    )
  data <- data[, -rm]

  # add row.names
  row.names(data) <- data$ID_AQEM
  data$ID_AQEM <- NULL

  # data with only two levels/integer/factor variables have to be numeric
  col_int_fac <-
    names(Filter(function(y) {
      is.integer(y) | is.factor(y)
    }, data))
  data[, col_int_fac] <- apply(data[, col_int_fac], 2, as.double)

  # rm development columns for now
  data$dev_hemimetabol <- NULL
  data$dev_holometabol <- NULL

  preproc_data[[region]] <- data
}

# save
saveRDS(
  object = preproc_data,
  file = file.path(
    data_cache,
    paste0("preproc_data_", data_type, ".rds")
  )
)