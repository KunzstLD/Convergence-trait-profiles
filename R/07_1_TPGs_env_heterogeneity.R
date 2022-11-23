# Load script 07_TPGs_env_heterogeneity.R
source(file.path(data_scr, "07_TPGs_env_heterogeneity.R"))

# TODO: Standardization by nr. of taxa within a family?
# Load and normalize mean trait profiles of TPGs
trait_CONT_mtpgs <- trait_CONT[, .(mean_affinity = mean(affinity)),
                               by = c("continent" , "group", "trait")]
trait_CONT_mtpgs <- trait_CONT_mtpgs[, .(continent,
                                         group,
                                         trait,
                                         mean_affinity)] %>%
  unique(.) %>%
  dcast(., ... ~ trait, value.var = "mean_affinity")
trait_CONT_mtpgs[, tpg_id := paste0(continent, "_", group)]
normalize_by_rowSum(trait_CONT_mtpgs,
                    non_trait_cols = c("continent",
                                       "group",
                                       "tpg_id"))
trait_CONT_mtpgs <- melt(
  trait_CONT_mtpgs,
  id.vars = c("continent", "group", "tpg_id"),
  value.name = "mean_affinity",
  variable.name = "trait"
)
trait_CONT_mtpgs[, grouping_feautre := sub("([a-z]{1,})(\\_)(.+)", "\\1", trait)]

# Trait data
# load script 07_TPGs_env_heterogeneity.R
trait_non_agg_lf <- list("EU" = Trait_EU_lf,
                         "NOA" = Trait_NOA_lf,
                         "AUS" = Trait_AUS_lf,
                         "NZ" = Trait_NZ_lf)

var_around_mean_tps <- list()
for(i in names(trait_non_agg_lf)){
  dat <- trait_non_agg_lf[[i]]
  
  # Quantify variation around mean trait profile per TPG
  # Merge back groups and mean trait profiles
  dat[trait_CONT[continent == i, ], group := i.group, on = "family"]
  dat[trait_CONT_mtpgs[continent == i,],
              mean_affinity := i.mean_affinity,
              on = c("group", "trait")]
  
  # Calculate distance to mean trait profile (i.e. absolute difference) for each trait
  # Then weight by the number of traits within each grouping feature  
  dat[, grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", trait)]
  dat[, weight := .N, by = c("taxa", "grouping_feature")]
  dat[!is.na(affinity),
              abs_diff_affinity := abs(affinity - mean_affinity)]
  dat[, dist_to_mean_tp := weighted.mean(abs_diff_affinity,
                                                 w = weight,
                                                 na.rm = TRUE),
              by = "taxa"]
  dat[, mean_dist_to_mean_tp_family := mean(dist_to_mean_tp, na.rm = TRUE), by = "family"]
  # Order families according to their distance to the mean tp of their tpg
  custom_order_families <- unique(dat[order(-mean_dist_to_mean_tp_family), family])
  dat[, family := factor(family, levels = custom_order_families)]
  
  # Cal. nr of taxa per family
  dat[, n_taxa_family := uniqueN(taxa), by = family]
  dat[is.na(affinity),
              n_taxa_used_family_grf := uniqueN(taxa),
              by = c("family", "grouping_feature")]
  
  # Plot distance of each taxa to the mean trait profile of their TPG per family 
  pl <- ggplot(dat,
         aes(x = family, y = dist_to_mean_tp)) +
    geom_boxplot(aes(fill = order)) +
    geom_point(
      size = 0.8,
      alpha = .5,
      position = position_jitter(seed = 1, width = .2),
      aes(color = order)
    ) +
    coord_flip() +
    scale_color_d3() +
    scale_fill_d3() +
    labs(x = "Family",
         y = "Distance to mean trait profiles of corresponding TPG",
         fill = "Order") +
    guides(color = "none") +
    theme_bw() +
    theme(
      axis.title = element_text(size = 16),
      axis.text.x = element_text(family = "Roboto Mono",
                                 size = 14),
      axis.text.y = element_text(family = "Roboto Mono",
                                 size = 14),
      legend.title = element_text(family = "Roboto Mono",
                                  size = 16),
      legend.text = element_text(family = "Roboto Mono",
                                 size = 14),
      strip.text = element_text(family = "Roboto Mono",
                                size = 14)
    )
  ggsave(
    plot = pl,
    filename = file.path(
      data_paper,
      "Graphs",
      paste0("Variation_mean_tp_family_", i, ".png")
    ),
    width = 35,
    height = 45,
    units = "cm"
  )
  
  # Save for further analysis & plots
  var_around_mean_tps[[i]] <- dat
}

# Create a data.table to work with
var_around_mean_tps_subset <- lapply(var_around_mean_tps, function(y)
  unique(y[, .(species,
               genus,
               family,
               order,
               taxa,
               dist_to_mean_tp,
               mean_dist_to_mean_tp_family)]))

var_around_mean_tps_subset <- rbindlist(var_around_mean_tps_subset, idcol = "continent")
var_around_mean_tps_subset[,continent := factor(continent, levels = c("NZ", "AUS", "EU", "NOA"))]

## Sd & ranges of distances to mean TPs ----
# Overall
var_around_mean_tps_subset[, sd(dist_to_mean_tp, na.rm = TRUE), by = "continent"]
var_around_mean_tps_subset[, range(dist_to_mean_tp, na.rm = TRUE), by = "continent"]

# By order
var_around_mean_tps_subset[, sd(dist_to_mean_tp, na.rm = TRUE),
                           by = c("continent", "order")]

# Distance of taxa to mean trait profile of their TPG per order 
ggplot(var_around_mean_tps_subset,
       aes(x = continent,
           y = dist_to_mean_tp)) +
  geom_boxplot(aes(fill = order)) +
  labs(x = "Continent or region",
       y = "Distance to mean trait profile of corresponding TPG",
       fill = "Order") +
  guides(color = "none") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(family = "Roboto Mono",
                               size = 14),
    axis.text.y = element_text(family = "Roboto Mono",
                               size = 14),
    legend.title = element_text(family = "Roboto Mono",
                                size = 16),
    legend.text = element_text(family = "Roboto Mono",
                               size = 14),
    strip.text = element_text(family = "Roboto Mono",
                              size = 14)
  ) +
  annotate(
    "segment",
    x = 0.5,
    y = -0.1,
    xend = 4.5,
    yend = -0.1,
    size = 5.5,
    linejoin = "mitre",
    arrow = arrow(type = "closed", length = unit(0.01, "npc"))
  ) +
  annotate(
    "text",
    x = 2.5,
    y = -0.1,
    label = "Environmental heterogeneity",
    color = "white",
    size = 4,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 1,
    y = -0.06,
    label = "KG-zones: 9",
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 2,
    y = -0.06,
    label = "KG-zones: 17",
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 3,
    y = -0.06,
    label = "KG-zones: 18",
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 4,
    y = -0.06,
    label = "KG-zones: 22",
    color = "black",
    size = 4,
    fontface = "bold"
  )
ggsave(
  filename = file.path(
    data_paper,
    "Graphs",
    "Variation_mean_tp_order.png"
  ),
  width = 35,
  height = 20,
  units = "cm"
)


# Distance of taxa to mean trait profile of their TPG per family 
pl_list <- list()
for (i in unique(var_around_mean_tps_subset$continent)) {
  pl_list[[i]] <-
    ggplot(var_around_mean_tps_subset[continent == i, ],
           aes(x = family, y = dist_to_mean_tp)) +
    geom_boxplot(aes(fill = order)) + # data =~ .x[n_taxa_family > 5, ],
    geom_point(
      size = 0.8,
      alpha = .5,
      position = position_jitter(seed = 1, width = .2),
      aes(color = order)
    ) +
    coord_flip() +
    scale_color_d3() +
    scale_fill_d3() +
    labs(x = "Family",
         y = "Distance to mean trait profiles of corresponding TPG",
         fill = "Order") +
    guides(color = "none") +
    theme_bw() +
    theme(
      axis.title = element_text(size = 16),
      axis.text.x = element_text(family = "Roboto Mono",
                                 size = 14),
      axis.text.y = element_text(family = "Roboto Mono",
                                 size = 14),
      legend.title = element_text(family = "Roboto Mono",
                                  size = 16),
      legend.text = element_text(family = "Roboto Mono",
                                 size = 14),
      strip.text = element_text(family = "Roboto Mono",
                                size = 14)
    )
}
lapply(names(pl_list), function(y)
  ggsave(
    filename = file.path(
      data_paper,
      "Graphs",
      paste0("Variation_mean_tp_family_", y, ".png"),
      width = 35,
      height = 20,
      units = "cm"
    )
  ))
pl_list$NZ


