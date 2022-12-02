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
trait_non_agg_lf <- list("EU" = Trait_EU_lf,
                         "NOA" = Trait_NOA_lf,
                         "AUS" = Trait_AUS_lf,
                         "NZ" = Trait_NZ_lf)
trait_non_agg_lf <- lapply(trait_non_agg_lf, function(y)
  y[, grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", trait)])

# Use only complete trait profiles
# Do not consider body form traits
trait_names <- unique(trait_CONT[grouping_feature %in% c("feed",
                                                         "resp",
                                                         "volt",
                                                         "size",
                                                         "locom"), trait]) #,"bf"
compl_taxa <- lapply(trait_non_agg_lf, function(y)
  y[grouping_feature %in% c("feed",
                            "resp",
                            "volt",
                            "size",
                            "locom"), ] %>%
    .[!is.na(affinity), .N, by = "taxa"] %>%
    .[N == length(trait_names), ]) 
trait_non_agg_lf_subset <- Map(function(x, y) x[taxa %in% y$taxa], trait_non_agg_lf,compl_taxa)
trait_non_agg_lf_subset <-
  lapply(trait_non_agg_lf_subset, function(y)
    y[grouping_feature %in% c("feed",
                              "resp",
                              "volt",
                              "size",
                              "locom"), ])

var_around_mean_tps <- list()
pl <- list()
for(i in names(trait_non_agg_lf_subset)) {
  dat <- trait_non_agg_lf_subset[[i]]
  
  # Quantify variation around mean trait profile per TPG
  # Merge back groups and mean trait profiles
  dat <- dat[family %in% trait_CONT[continent == i, family],]
  dat[trait_CONT[continent == i, ], group := i.group, on = "family"]
  dat[trait_CONT_mtpgs[continent == i,],
      mean_affinity := i.mean_affinity,
      on = c("group", "trait")]
  
  # Overlap
  # Then weight by the number of traits within each grouping feature
  dat[, grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", trait)]
  dat[, overlap_per_trait := fun_overlap_ind(p = affinity,
                                             q = mean_affinity),
      by = c("taxa", "grouping_feature")]
  dat[, dist_to_mean_tp := sqrt(mean(overlap_per_trait)), by = "taxa"]
  
  # Alternatively, absolute distance for each trait
  # dat[, weight := .N, by = c("taxa", "grouping_feature")]
  # dat[!is.na(affinity),
  #     abs_diff_affinity := abs(affinity - mean_affinity)]
  # dat[, dist_to_mean_tp := weighted.mean(abs_diff_affinity,
  #                                        w = weight,
  #                                        na.rm = TRUE),
  #     by = "taxa"]
  
  
  dat[, mean_dist_to_mean_tp_family := mean(dist_to_mean_tp, na.rm = TRUE), by = "family"]
  # Order families according to their distance to the mean tp of their tpg
  custom_order_families <-
    unique(dat[order(-mean_dist_to_mean_tp_family), family])
  dat[, family := factor(family, levels = custom_order_families)]
  
  # Cal. nr of taxa per family
  dat[, n_taxa_family := uniqueN(taxa), by = family]
  dat[is.na(affinity),
      n_taxa_used_family_grf := uniqueN(taxa),
      by = c("family", "grouping_feature")]
  
  # Plot distance of each taxa to the mean trait profile of their TPG per family
  pl[[i]] <- ggplot(dat,
                    aes(x = family, y = dist_to_mean_tp)) +
    geom_boxplot(aes(fill = order)) +
    geom_point(
      size = 0.8,
      alpha = .5,
      position = position_jitter(seed = 1, width = .2),
      aes(color = order)
    ) +
    ggtitle(fcase(i == "NOA", "NA",
                  i == "EU", "EUR",
                  default = i)) +
    coord_flip() +
    annotate(
      "text",
      x = fcase(i == "NZ", 62,
                i == "AUS", 99,
                i == "EU", 80,
                i == "NOA", 82),
      y = fcase(i == "NZ", 0.75,
                i == "AUS", 0.75,
                i == "EU", 0.7,
                i == "NOA", 0.8),
      label = paste(
        "KG-zones:",
        fcase(i == "NZ", 9,
              i == "AUS", 17,
              i == "EU", 18,
              i == "NOA", 22)
      ),
      color = "black",
      size = 4,
      fontface = "bold"
    ) +
    scale_color_d3() +
    scale_fill_d3() +
    labs(x = "Family",
         y = "Distance to mean trait profile of corresponding TPG",
         fill = "Order") +
    guides(color = "none") +
    theme_bw() +
    theme(
      axis.title = element_text(size = 14),
      axis.text.x = element_text(family = "Roboto Mono",
                                 size = 12),
      axis.text.y = element_text(family = "Roboto Mono",
                                 size = 6),
      legend.title = element_text(family = "Roboto Mono",
                                  size = 14),
      legend.text = element_text(family = "Roboto Mono",
                                 size = 12),
      strip.text = element_text(family = "Roboto Mono",
                                size = 14)
    )
  # ggsave(
  #   plot = pl,
  #   filename = file.path(
  #     data_paper,
  #     "Graphs",
  #     paste0("Variation_mean_tp_family_", i, ".png")
  #   ),
  #   width = 35,
  #   height = 45,
  #   units = "cm"
  # )
  
  # Save for further analysis & plots
  var_around_mean_tps[[i]] <- dat
}
# lapply(var_around_mean_tps, function(y) y[is.na(dist_to_mean_tp), ])

# Plotting
# ?Add information on how many taxa could be used
pl$NZ + pl$AUS + pl$EU + pl$NOA +
  plot_layout(guides = "collect")
ggsave(
  filename = file.path(
    data_paper,
    "Graphs",
    paste0("Variation_mean_tp_family_combined.png")
  ),
  width = 45,
  height = 35,
  units = "cm"
)

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
var_around_mean_tps_subset[, n_taxa := .N, by = "continent"]

## Sd & ranges of distances to mean TPs ----
# Overall
var_around_mean_tps_subset[, mean(dist_to_mean_tp, na.rm = TRUE), by = "continent"]
var_around_mean_tps_subset[, median(dist_to_mean_tp, na.rm = TRUE), by = "continent"]
var_around_mean_tps_subset[, sd(dist_to_mean_tp, na.rm = TRUE), by = "continent"]
var_around_mean_tps_subset[, range(dist_to_mean_tp, na.rm = TRUE), by = "continent"] %>% 
  .[, .(continent, round(V1, digits = 3))]

# By order
var_around_mean_tps_subset[, sd(dist_to_mean_tp, na.rm = TRUE),
                           by = c("continent", "order")]

# N families orig vs. after subsetting to complete trait profiles
trait_CONT[, uniqueN(family), by = "continent"]
var_around_mean_tps_subset[, uniqueN(family), by = "continent"]

data_av <- list()
for(i in names(trait_non_agg_lf)) {
  data_av[[i]] <-
    trait_non_agg_lf[[i]][family %in% var_around_mean_tps_subset[continent == i, family],
                          uniqueN(taxa), by = "family"]
}
data_av <- rbindlist(data_av, idcol = "continent")
setnames(data_av, "V1", "taxa_per_family_orignal")
data_av[var_around_mean_tps_subset[, uniqueN(taxa), by = c("continent", "family")],
        taxa_per_family_complete_tp := V1,
        on = c("continent", "family")]
data_av[, perc_taxa_used := round(100*(sum(taxa_per_family_complete_tp) / sum(taxa_per_family_orignal))),
        by = "continent"]

# Variation within TPGs summary plot
var_around_mean_tps_subset[continent == "EU", ] %>% 
  .[dist_to_mean_tp == min(dist_to_mean_tp), ]

set.seed(1234)
ggplot(var_around_mean_tps_subset,
       aes(x = continent,
           y = dist_to_mean_tp)) +
  geom_violin() +
  geom_jitter(width = 0.06,
              alpha = 0.5) +
  stat_summary(fun = "mean", color = "red") +
  labs(x = "Continent or region",
       y = "Distance taxon trait profile to mean trait profile of its family TPG") +
  scale_x_discrete(labels = c("NZ", "AUS", "EUR", "NA")) +
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
    y = -0.2,
    xend = 4.5,
    yend = -0.2,
    size = 5.5,
    linejoin = "mitre",
    arrow = arrow(type = "closed", length = unit(0.01, "npc"))
  ) +
  annotate(
    "text",
    x = 2.5,
    y = -0.2,
    label = "Environmental heterogeneity",
    color = "white",
    size = 4,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 1,
    y = -0.1,
    label = "KG-zones: 9",
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 2,
    y = -0.1,
    label = "KG-zones: 17",
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 3,
    y = -0.1,
    label = "KG-zones: 17",
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 4,
    y = -0.1,
    label = "KG-zones: 22",
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  annotate(
    "curve",
    x = 2.98,
    y = 0.0651,
    xend = 2.85,
    yend = 0.05,
    curvature = 0.2,
    arrow = arrow(length = unit(2, "mm"))
  ) +
  annotate("text",
           x = 2.75,
           y = 0.031,
           label = "Orthetrum (Libellulidae)")
ggsave(
  filename = file.path(
    data_paper,
    "Graphs",
    "Variation_within_tpgs.png"
  ),
  width = 35,
  height = 20,
  units = "cm"
)
