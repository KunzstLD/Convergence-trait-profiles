# __________________________________________________________________________
# Data preparation
# TODO: col_two_levels actually selects also for columns with less than
# two levels
# ___________________________________________________________________________

# Load trait data
trait_data_ww <- load_data(pattern = ".*agg\\.rds", path = data_in)

# Check if trait datasets have the same colnames
trait_data_ww %>% check_colNames()

# Rm dev for now
# Looking at the trait distributions revealed that in all databases no parasites were present 
trait_data_ww <- lapply(
  trait_data_ww,
  function(y) y[, c("dev_holometabol", 
                    "dev_hemimetabol", 
                    "feed_parasite") := NULL]
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

# ---- Overview over orders ------------------------------------------------------------------------

# Bind trait data
trait_data_bind <- rbindlist(trait_data_ww, idcol = "continent")

# Calculate nr of taxa & prop of orders
trait_data_bind[, nr_taxa := .N, by = "continent"]
trait_data_bind[, prop_order := .N/nr_taxa, by = .(continent, order)]

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(trait_data_bind, aes(x = as.factor(order),
                            y = prop_order * 100)) +
  geom_pointrange(aes(
    ymin = 0,
    ymax = prop_order * 100,
    color = as.factor(continent)
  ),
  position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = cbbPalette,
                      labels = c(
                        paste0("AUS, n =", unique(trait_data_bind[continent == "AUS", nr_taxa])),
                        paste0("EU, n =", unique(trait_data_bind[continent == "EU", nr_taxa])),
                        paste0("NOA, n =", unique(trait_data_bind[continent == "NOA", nr_taxa])),
                        paste0("NZ, n =", unique(trait_data_bind[continent == "NZ", nr_taxa]))
                      )) +
  labs(x = "Orders",
       y = "Percentage per trait dataset",
       color = "Continent") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(
      family = "Roboto Mono",
      size = 14,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(family = "Roboto Mono",
                               size = 14),
    legend.title = element_text(family = "Roboto Mono",
                                size = 16),
    legend.text = element_text(family = "Roboto Mono",
                               size = 14),
    strip.text = element_text(family = "Roboto Mono",
                              size = 14),
    panel.grid = element_blank()
  )
ggsave(
  filename = file.path(data_paper,
                       "Graphs",
                       "Perc_order_per_continent.png"),
  width = 35,
  height = 25,
  units = "cm"
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
