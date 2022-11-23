# ________________________________________________________________________
# Overview over aggregated trait datasets ----
# ________________________________________________________________________

## Continental datasets ----
trait_data_bind <- readRDS(file.path(data_cache,
                                     "trait_data_ww_bind.rds"))

# Calculate nr of taxa & prop of orders
trait_data_bind[, nr_taxa := .N, by = "continent"]
trait_data_bind[, nr_fam_per_order := .N, by = .(continent, order)]
trait_data_bind[, prop_order := nr_fam_per_order / nr_taxa]

# Plot (will rather not end up in publication)
cbbPalette <-
  c(
    "#000000",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"
  )
ggplot(trait_data_bind, aes(
  x = as.factor(order),
  y = prop_order * 100
)) +
  geom_pointrange(aes(
    ymin = 0,
    ymax = prop_order * 100,
    color = as.factor(continent)
  ),
  position = position_dodge(width = 0.4)
  ) +
  scale_colour_manual(
    values = cbbPalette,
    labels = c(
      paste0("AUS, n =", unique(trait_data_bind[continent == "AUS", nr_taxa])),
      paste0("EU, n =", unique(trait_data_bind[continent == "EU", nr_taxa])),
      paste0("NOA, n =", unique(trait_data_bind[continent == "NOA", nr_taxa])),
      paste0("NZ, n =", unique(trait_data_bind[continent == "NZ", nr_taxa]))
    )
  ) +
  labs(
    x = "Orders",
    y = "Percentage per trait dataset",
    color = "Continent"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(
      family = "Roboto Mono",
      size = 14,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    legend.title = element_text(
      family = "Roboto Mono",
      size = 16
    ),
    legend.text = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    strip.text = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    panel.grid = element_blank()
  )
ggsave(
  filename = file.path(
    data_paper,
    "Graphs",
    "Perc_order_per_continent.png"
  ),
  width = 35,
  height = 25,
  units = "cm"
)










