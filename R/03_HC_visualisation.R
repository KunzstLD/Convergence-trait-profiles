#___________________________________________________________________________________________________
#### Hierarchical clustering results & Visualization ####
# __________________________________________________________________________________________________

#### Loading data ####

# Results from HC
hc_output_ww <- readRDS(file.path(data_cache,
                                  "hc_output_ww.rds"))

# Load trait data (with order information)
trait_data_ww <- load_data(pattern = "agg*.rds", path = data_in)

# Rm dev stage 
trait_data_ww <- lapply(trait_data_ww, function(y) y[, c("dev_hemimetabol", 
                                                         "dev_holometabol", 
                                                         "feed_parasite") := NULL])

# Desired order of traits (for plotting purposes later on)
des_traits_order <- c(
  "bf_cylindrical",
  "bf_flattened",
  "bf_spherical",
  "bf_streamlined",
  "feed_filter",
  "feed_gatherer",
  "feed_herbivore",
  "feed_predator",
  "feed_shredder",
  "locom_burrow",
  "locom_crawl",
  "locom_sessil",
  "locom_swim",
  "ovip_aqu",
  "ovip_ter",
  "ovip_ovo",
  "resp_gil",
  "resp_pls_spi",
  "resp_teg",
  "size_small",
  "size_medium",
  "size_large",
  "volt_semi",
  "volt_uni",
  "volt_bi_multi"
)

# __________________________________________________________________________________________________
#### Grouping features ####
# __________________________________________________________________________________________________

# Extract Grouping features that most contributed to the global distance
lookup_gf <- data.table(
  code = c("F1",
           "F2",
           "F3",
           "F4",
           "F5",
           "F6",
           "F7"),
  gf = c(
    "feeding mode",
    "respiration",
    "voltinism",
    "locomotion",
    "oviposition",
    "size",
    "body form"
  )
)

global_dist <- list()
for (i in names(hc_output_ww)) {
  gf_intm <- data.table(
    hc_output_ww[[i]]$contribu_global_dist$glocor,
    "code" = rownames(hc_output_ww[[i]]$contribu_global_dist$glocor)
  )
  gf_intm[lookup_gf, grouping_feature := i.gf, on = "code"]
  global_dist[[i]] <- gf_intm[order(-`global distance`),]
}
saveRDS(object = global_dist,
        file = file.path(data_cache, "global_dist.rds"))

# __________________________________________________________________________________________________
#### Create plots ####
# __________________________________________________________________________________________________

# _________________________________________________________
#### Gap Statistic ####
# _________________________________________________________
for(region in names(hc_output_ww)) {
  png(
    file = file.path(
      data_paper,
      "Graphs",
      paste0("Gap_statistic_", region, ".png")
    ),
    width = 1300,
    height = 1100,
    res = 100
  )
  plot(
    hc_output_ww[[region]]$gap_statistic,
    main = paste(
      "GAP",
      region,
      "\n Optimal number of groups:",
      hc_output_ww[[region]]$optimal_nog
    )
  )
  dev.off()
}

# _______________________________________________________
#### Dendrogramm #### 
# TODO: Dendrogram with order as labels
# TODO: A few family names are too long, 
# needs visual improvement
# _______________________________________________________
dendrograms <- list()
for (i in names(hc_output_ww)) {
  plot <- fun_dendrog_pl(
    hc = hc_output_ww[[i]]$hc_element,
    optimal_nog = hc_output_ww[[i]]$optimal_nog,
    labels = hc_output_ww[[i]]$labels_dendrogram
  )
  dendrograms[[i]] <- plot
}

# save dendrogram plots
for (i in names(dendrograms)) {
  png(
    file = file.path(data_paper,
                     "Graphs",
                     paste0(
                       "Dendrogram_", i, ".png"
                     )),
    width = 1100,
    height = 1300,
    res = 100
  )
  plot(dendrograms[[i]], horiz = TRUE)
  dev.off()
}

# _______________________________________________________
#### Heatmap plots ####
# _______________________________________________________

# trait profile groups for all continents
trait_profile_groups <- list()
for(i in names(hc_output_ww)) {
  trait_profile_groups[[i]] <-
    data.table(
      family = names(
        cutree(
          dendrograms[[i]],
          k = hc_output_ww[[i]]$optimal_nog,
          order_clusters_as_data = FALSE
        )
      ),
      group = cutree(
        dendrograms[[i]],
        k = hc_output_ww[[i]]$optimal_nog,
        order_clusters_as_data = FALSE
      )
    )
}

#*****************************************************#
#--- AUS ----
#*****************************************************#

# Trait data
trait_AUS <- trait_data_ww$Trait_AUS_agg.rds
trait_AUS[trait_profile_groups$AUS, group := i.group,
     on = "family"]

# Change col order
setcolorder(trait_AUS, des_traits_order)

# Save for RF analyses
saveRDS(object = trait_AUS,
        file = file.path(data_cache, "trait_AUS_with_groups.rds"))

# Order families according to dendrogram
trait_AUS[, family := factor(family, levels = labels(dendrograms$AUS))]
trait_AUS_lf <- melt(
  trait_AUS,
  id.vars = c("family", "order", "group"),
  variable.name = "trait",
  value.name = "affinity"
)
trait_AUS_lf[, grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", trait)] 
trait_AUS_lf[, group := factor(group)]

# Create trait label for plotting
trait_AUS_lf[, trait_label := as.character(sub("([a-z]{1,})(\\_)(.+)", "\\3", trait))]
trait_AUS_lf[, trait_label := factor(trait_label, levels = unique(trait_label))]

# names for grouping features 
grouping_feature_names <- c(
  "feed" = "Feeding mode",
  "locom" = "Locomotion",
  "resp" = "Respiration",
  "size" = "Body size",
  "volt" = "Voltinism",
  "bf" = "Body form",
  "ovip" = "Oviposition",
  "1" = "TPG 1",
  "2" = "TPG 2",
  "3" = "TPG 3",
  "4" = "TPG 4",
  "5" = "TPG 5",
  "6" = "TPG 6",
  "7" = "TPG 7"
)

# plot
fun_heatmap_single_cont(data = trait_AUS_lf)+
  ggtitle("TPGs AUS")
ggplot2::ggsave(
  filename = file.path(data_paper,
                       "Graphs", 
                       "Heatmap_tpgs_AUS.png"),
  width = 45,
  height = 33,
  units = "cm"
)

# Orders (do we have phylogenetic signal?)
trait_AUS[, nr_families_gr := .N, by = group]
trait_AUS[, prop_order := .N/nr_families_gr, by = .(group, order)]

# Label function for plot (and following)
tpg_names <- c("1" = "TPG 1",
               "2" = "TPG 2",
               "3" = "TPG 3",
               "4" = "TPG 4",
               "5" = "TPG 5",
               "6" = "TPG 6",
               "7" = "TPG 7",
               "8" = "TPG 8", 
               "9" = "TGP 9", 
               "10" = "TPG 10")

# Annotate nr of families per TPG 
annotations_aus <- unique(trait_AUS[order(group), .(group, nr_families_gr)])
annotations_aus[, `:=`(
  x = rep(7.3, 7),
  y = c(7.1, 6.1, 5.1, 4.1, 3.1, 2.1, 1.5),
  label = paste0("n = ", nr_families_gr)
)]

ggplot(trait_AUS, aes(x = as.factor(order),
                      y = prop_order * 100)) +
  geom_pointrange(aes(ymin = 0,
                      ymax = prop_order * 100))  +
  geom_text(
    data = annotations_aus,
    aes(x = x,
        y = y,
        label = label),
    colour = "black",
    inherit.aes = FALSE,
    parse = FALSE,
    size = 5
  ) +
  facet_grid(group ~ ., labeller = as_labeller(tpg_names)) +
  labs(x = "Order", y = "Percentage per TPG") +
  ggtitle("AUS: Orders per TPGs") +
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
ggplot2::ggsave(
  filename = file.path(data_paper,
                       "Graphs", 
                       "Orders_per_TPG_AUS.png"),
  width = 35,
  height = 33,
  units = "cm"
)

#*****************************************************#
#---- EU ----
#*****************************************************#

# Trait data
trait_EU <- trait_data_ww$Trait_EU_agg.rds
trait_EU[trait_profile_groups$EU, group := i.group,
          on = "family"]

# Re-order cols 
setcolorder(trait_EU, des_traits_order)

# Save for RF analyses
saveRDS(object = trait_EU,
        file = file.path(data_cache, "trait_EU_with_groups.rds"))

# Order according to dendrogram
trait_EU[, family := factor(family, levels = labels(dendrograms$EU))]

trait_EU_lf <- melt(
  trait_EU,
  id.vars = c("family", "order", "group"),
  variable.name = "trait",
  value.name = "affinity"
)
trait_EU_lf[, grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", trait)] 
trait_EU_lf[, group := factor(group)]

# Create trait label for plotting
trait_EU_lf[, trait_label := as.character(sub("([a-z]{1,})(\\_)(.+)", "\\3", trait))]
trait_EU_lf[, trait_label := factor(trait_label, levels = unique(trait_label))]

# names for grouping features 
grouping_feature_names <- c(
  "feed" = "Feeding mode",
  "locom" = "Locomotion",
  "resp" = "Respiration",
  "size" = "Body size",
  "volt" = "Voltinism",
  "bf" = "Body form",
  "ovip" = "Oviposition",
  "1" = "TPG 1",
  "2" = "TPG 2",
  "3" = "TPG 3",
  "4" = "TPG 4",
  "5" = "TPG 5",
  "6" = "TPG 6",
  "7" = "TPG 7",
  "8" = "TPG 8",
  "9" = "TPG 9",
  "10" = "TPG 10"
)

# plot
fun_heatmap_single_cont(data = trait_EU_lf) +
  ggtitle("TPGs EU")
ggplot2::ggsave(
  filename = file.path(data_paper,
                       "Graphs", 
                       "Heatmap_tpgs_EU.png"),
  width = 47,
  height = 33,
  units = "cm"
)

# Orders (do we have phylogenetic signal?)
trait_EU[, nr_families_gr := .N, by = group]
trait_EU[, prop_order := .N/nr_families_gr, by = .(group, order)]

annotations_eu <- unique(trait_EU[order(group), .(group, nr_families_gr)])
annotations_eu[, `:=`(
  x = rep(7.3, 10),
  y = c(10.1, 9.1, 8.1, 7.1, 6.1, 5.1, 4.3, 3.3, 2.3, 1.5),
  label = paste0("n = ", nr_families_gr)
)]

ggplot(trait_EU, aes(x = as.factor(order),
                      y = prop_order * 100)) +
  geom_pointrange(aes(ymin = 0,
                      ymax = prop_order * 100))  +
  geom_text(
    data = annotations_eu,
    aes(x = x,
        y = y,
        label = label),
    colour = "black",
    inherit.aes = FALSE,
    parse = FALSE,
    size = 5
  ) +
  facet_grid(group ~ ., labeller = as_labeller(tpg_names)) +
  labs(x = "Order", y = "Percentage per TPG") +
  ggtitle("EU: Orders per TPGs") +
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
ggplot2::ggsave(
  filename = file.path(data_paper,
                       "Graphs", 
                       "Orders_per_TPG_EU.png"),
  width = 35,
  height = 33,
  units = "cm"
)

#*****************************************************#
#---- NOA ----
#*****************************************************#

# Trait data
trait_NOA <- trait_data_ww$Trait_NOA_agg.rds
trait_NOA[trait_profile_groups$NOA, group := i.group,
         on = "family"]

# Re-order columns
setcolorder(trait_NOA, des_traits_order)

# Save for RF analyses
saveRDS(object = trait_NOA,
        file = file.path(data_cache, "trait_NOA_with_groups.rds"))

# Order families according to dendrogram
trait_NOA[, family := factor(family, levels = labels(dendrograms$NOA))]

trait_NOA_lf <- melt(
  trait_NOA,
  id.vars = c("family", "order", "group"),
  variable.name = "trait",
  value.name = "affinity"
)
trait_NOA_lf[, grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", trait)] 
trait_NOA_lf[, group := factor(group)]

# Create trait label for plotting
trait_NOA_lf[, trait_label := as.character(sub("([a-z]{1,})(\\_)(.+)", "\\3", trait))]
trait_NOA_lf[, trait_label := factor(trait_label, levels = unique(trait_label))]

# Names for grouping features 
grouping_feature_names <- c(
  "feed" = "Feeding mode",
  "locom" = "Locomotion",
  "resp" = "Respiration",
  "size" = "Body size",
  "volt" = "Voltinism",
  "bf" = "Body form",
  "ovip" = "Oviposition",
  "1" = "TPG 1",
  "2" = "TPG 2",
  "3" = "TPG 3",
  "4" = "TPG 4",
  "5" = "TPG 5",
  "6" = "TPG 6"
)

# plot
fun_heatmap_single_cont(data = trait_NOA_lf) +
  ggtitle("TPGs NOA")
ggplot2::ggsave(
  filename = file.path(data_paper,
                       "Graphs",
                       "Heatmap_tpgs_NOA.png"),
  width = 50,
  height = 33,
  units = "cm"
)

# Orders (do we have phylogenetic signal?)
trait_NOA[, nr_families_gr := .N, by = group]
trait_NOA[, prop_order := .N/nr_families_gr, by = .(group, order)]

annotations_noa <- unique(trait_NOA[order(group), .(group, nr_families_gr)])
annotations_noa[, `:=`(
  x = rep(9.3, 6),
  y = c(6.1:1.1),
  label = paste0("n = ", nr_families_gr)
)]

ggplot(trait_NOA, aes(x = as.factor(order),
                     y = prop_order * 100)) +
  geom_pointrange(aes(ymin = 0,
                      ymax = prop_order * 100))  +
  geom_text(
    data = annotations_noa,
    aes(x = x,
        y = y,
        label = label),
    colour = "black",
    inherit.aes = FALSE,
    parse = FALSE,
    size = 5
  ) +
  facet_grid(group ~ ., labeller = as_labeller(tpg_names)) +
  labs(x = "Order", y = "Percentage per TPG") +
  ggtitle("NOA: Orders per TPGs") +
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
ggplot2::ggsave(
  filename = file.path(data_paper,
                       "Graphs", 
                       "Orders_per_TPG_NOA.png"),
  width = 35,
  height = 33,
  units = "cm"
)


#*****************************************************#
#---- NZ ----
#*****************************************************#

# Trait data
trait_NZ <- trait_data_ww$Trait_NZ_agg.rds
trait_NZ[trait_profile_groups$NZ, group := i.group,
          on = "family"]

# Re-order columns
setcolorder(trait_NZ, des_traits_order)

# Save for RF analyses
saveRDS(object = trait_NZ,
        file = file.path(data_cache, "trait_NZ_with_groups.rds"))

# Order according to dendrogram
trait_NZ[, family := factor(family, levels = labels(dendrograms$NZ))]

trait_NZ_lf <- melt(
  trait_NZ,
  id.vars = c("family", "order", "group"),
  variable.name = "trait",
  value.name = "affinity"
)
trait_NZ_lf[, grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", trait)] 
trait_NZ_lf[, group := factor(group)]

# Create trait label for plotting
trait_NZ_lf[, trait_label := as.character(sub("([a-z]{1,})(\\_)(.+)", "\\3", trait))]
trait_NZ_lf[, trait_label := factor(trait_label, levels = unique(trait_label))]

# names for grouping features 
grouping_feature_names <- c(
  "feed" = "Feeding mode",
  "locom" = "Locomotion",
  "resp" = "Respiration",
  "size" = "Body size",
  "volt" = "Voltinism",
  "bf" = "Body form",
  "ovip" = "Oviposition",
  "1" = "TPG 1",
  "2" = "TPG 2",
  "3" = "TPG 3",
  "4" = "TPG 4",
  "5" = "TPG 5",
  "6" = "TPG 6",
  "7" = "TPG 7",
  "8" = "TPG 8",
  "9" = "TPG 9",
  "10" = "TPG 10"
)

# plot
fun_heatmap_single_cont(data = trait_NZ_lf) +
  ggtitle("TPGs NZ")
ggplot2::ggsave(
  filename = file.path(data_paper,
                       "Graphs", 
                       "Heatmap_tpgs_NZ.png"),
  width = 45,
  height = 33,
  units = "cm"
)

# Orders (do we have phylogenetic signal?)
trait_NZ[, nr_families_gr := .N, by = group]
trait_NZ[, prop_order := .N/nr_families_gr, by = .(group, order)]

annotations_nz <- unique(trait_NZ[order(group), .(group, nr_families_gr)])
annotations_nz[, `:=`(
  x = rep(9.3, 9),
  y = c(9.1:1.1),
  label = paste0("n = ", nr_families_gr)
)]

ggplot(trait_NZ, aes(x = as.factor(order),
                      y = prop_order * 100)) +
  geom_pointrange(aes(ymin = 0,
                      ymax = prop_order * 100))  +
  geom_text(
    data = annotations_nz,
    aes(x = x,
        y = y,
        label = label),
    colour = "black",
    inherit.aes = FALSE,
    parse = FALSE,
    size = 5
  ) +
  facet_grid(group ~ ., labeller = as_labeller(tpg_names)) +
  labs(x = "Order", y = "Percentage per TPG") +
  ggtitle("NZ: Orders per TPGs") +
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
ggplot2::ggsave(
  filename = file.path(data_paper,
                       "Graphs", 
                       "Orders_per_TPG_NZ.png"),
  width = 35,
  height = 33,
  units = "cm"
)

# Save trait data with group assignments for further analyses
trait_CONT <- rbind(trait_AUS_lf,
                    trait_EU_lf,
                    trait_NOA_lf,
                    trait_NZ_lf,
                    idcol = "continent")
trait_CONT[, continent := factor(continent, labels = c("AUS", "EU", "NOA", "NZ"))]
saveRDS(object = trait_CONT,
        file = file.path(data_cache, "trait_dat_grp_assig.rds"))
