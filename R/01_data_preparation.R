# __________________________________________________________________________
# Data preparation
# TODO: col_two_levels actually selects also for columns with less than
# two levels
# ___________________________________________________________________________

# Load trait data
trait_data_ww <- load_data(pattern = "*.rds", path = data_in)

# Check if trait datasets have the same colnames
trait_data_ww %>% check_colNames()

# Rm dev for now
trait_data_ww <- lapply(
  trait_data_ww,
  function(y) y[, c("dev_holometabol", "dev_hemimetabol") := NULL]
)

# rename file list
names(trait_data_ww) <- sub(
  "([A-z]{5})(\\_)([A-z]{2,})(\\_)([A-z]{3})(.+)",
  "\\3",
  names(trait_data_ww)
)

# Arrange colnames according to traits
lapply(
  trait_data_ww,
  function(y) {
    setcolorder(y, c(
      "family",
      "feed_shredder",
      "feed_gatherer",
      "feed_filter",
      "feed_predator",
      "feed_herbivore",
      "feed_parasite",
      "resp_teg",
      "resp_gil",
      "resp_pls_spi",
      "volt_semi",
      "volt_uni",
      "volt_bi_multi",
      "locom_swim",
      "locom_crawl",
      "locom_burrow",
      "locom_sessil",
      "ovip_ter",
      "ovip_ovo",
      "ovip_aqu",
      "size_medium",
      "size_small",
      "size_large",
      "bf_streamlined",
      "bf_flattened",
      "bf_cylindrical",
      "bf_spherical"
    ))
  }
)

# ---- Data preparation for HC -----------------------------------------------

preproc_traits <- list()
for (region in c("AUS", "EU", "NOA", "NZ")) {
  data <- trait_data_ww[[region]]
  data[, order := NULL]

  # convert to data.frame -> data table does not support row.names
  setDF(data)

  # add row.names
  row.names(data) <- data$family
  data$family <- NULL

  # data with only two levels/integer/factor variables have to be numeric
  col_int_fac <-
    names(Filter(function(y) {
      is.integer(y) | is.factor(y)
    }, data))
  data[, col_int_fac] <- apply(data[, col_int_fac], 2, as.double)

  # Convert to ktab object
  vec <- sub("\\_.*", "\\1", names(data))
  blocks <- rle(vec)$lengths
  data <- prep.fuzzy(data, blocks)
  data <- ktab.list.df(list(data))

  # save to list
  preproc_traits[[region]] <- data
}

saveRDS(
  object = preproc_traits,
  file = file.path(data_cache, "preproc_traits.rds")
)