# ______________________________________________________________
# Data preparation
# TODO: col_two_levels actually selects also for columns with less than
# two levels
# ______________________________________________________________

# read in data
trait_subsets <- readRDS(file.path(data_cache, "trait_subsets.rds"))

# omit data we (for now) do not consider (development)
trait_subsets[, c("dev_hemimetabol", "dev_holometabol") := NULL]

# data preprocessing
preproc_data <- list()
for (region in c("arid", "temperate", "cold", "polar")) {
  data <- trait_subsets[climateregion == region, ]

  # ----- Trait Aggregation --------------------------------------------------
  # columns not to consider for aggregation
  non_trait_cols <- grep(
      "ER[0-9]{1,}|species|genus|family|order|tax.*|climateregion|ID.*",
      names(data),
      value = TRUE
    )

  # trait aggregation to genus lvl (for comparability reasons with NOA data)
  data <- direct_agg(
    trait_data = data[!is.na(species), ],
    non_trait_cols = non_trait_cols,
    method = median,
    taxon_lvl = "genus",
    na.rm = TRUE
  )

  # ---- Data preparation for HC -----------------------------------------------
  # convert to data.frame -> data table does not support row.names
  setDF(data)

  # add row.names
  row.names(data) <- data$genus
  data$genus <- NULL

  # data with only two levels/integer/factor variables have to be numeric
  col_int_fac <-
    names(Filter(function(y) {
      is.integer(y) | is.factor(y)
    }, data))
  data[, col_int_fac] <- apply(data[, col_int_fac], 2, as.double)

  # save to list
  preproc_data[[region]] <- data
}

# save
saveRDS(
  object = preproc_data,
  file = file.path(
    data_cache,
    paste0("preproc_data_", "genus", ".rds")
  )
)