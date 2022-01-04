# ___________________________________________________________________________
# Hierarchical clustering results & Visualization ----
# ___________________________________________________________________________

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

# ___________________________________________________________________________
# Extract global distances ----
# ___________________________________________________________________________

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

# ___________________________________________________________________________
# Auxiliary plots ----
# ___________________________________________________________________________

## Gap Statistic ----
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

## Dendrogramm ----
# TODO: Dendrogram with order as labels
# TODO: A few family names are too long, 
# needs visual improvement
dendrograms <- list()
for (i in names(hc_output_ww)) {
  plot <- fun_dendrog_pl(
    hc = hc_output_ww[[i]]$hc_wardD,
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

## Ordination plots with superimposed single linkage clustering ----
for(i in names(hc_output_ww)) {
  png(
    file = file.path(
      data_paper,
      "Graphs",
      paste0("Ordination_superimposed_hc_", i, ".png")
    ),
    width = 1100,
    height = 1300,
    res = 100
  )
  s.label(hc_output_ww[[i]]$pcoa$li, clabel = 0.8)
  # plot(hc_output_ww[[i]]$pcoa$li)
  ordicluster(
    hc_output_ww[[i]]$pcoa$li,
    hc_output_ww[[i]]$hc_single,
    prune = 3,
    col = cutree(hc_output_ww[[i]]$hc_single,
                 hc_output_ww[[i]]$optimal_nog)
  )
  dev.off()
}

# _______________________________________________________
# Heatmap plots ----
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
## AUS ----
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
  "6" = "TPG 6"
)

# Plot TPGs and expressed traits
fun_heatmap_single_cont(data = trait_AUS_lf)+
  ggtitle("TPGs AUS")
ggplot2::ggsave(
  filename = file.path(data_paper,
                       "Graphs", 
                       "Heatmap_tpgs_AUS.png"),
  width = 50,
  height = 40,
  units = "cm"
)

#*****************************************************#
## EU ----
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
  "9" = "TPG 9"
)

# plot
fun_heatmap_single_cont(data = trait_EU_lf) +
  ggtitle("TPGs EUR")
ggplot2::ggsave(
  filename = file.path(data_paper,
                       "Graphs", 
                       "Heatmap_tpgs_EU.png"),
  width = 50,
  height = 40,
  units = "cm"
)

#*****************************************************#
## NOA ----
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
  "6" = "TPG 6",
  "7" = "TPG 7",
  "8" = "TPG 8",
  "9" = "TPG 9"
)

# plot
fun_heatmap_single_cont(data = trait_NOA_lf) +
  ggtitle("TPGs NA")
ggplot2::ggsave(
  filename = file.path(data_paper,
                       "Graphs",
                       "Heatmap_tpgs_NOA.png"),
  width = 50,
  height = 40,
  units = "cm"
)

#*****************************************************#
## NZ ----
# May merge cluster 6 to 7? Cluster 6 consists of only two families;
# similar defining traits
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
  width = 60,
  height = 35,
  units = "cm"
)

# Save trait data with group assignments for further analyses ----
trait_CONT <- rbind(trait_AUS_lf,
                    trait_EU_lf,
                    trait_NOA_lf,
                    trait_NZ_lf,
                    idcol = "continent")
trait_CONT[, continent := factor(continent, labels = c("AUS", "EU", "NOA", "NZ"))]
saveRDS(object = trait_CONT,
        file = file.path(data_cache, "trait_dat_grp_assig.rds"))

# _______________________________________________________
# Outliers ----
# Results from single linkage clustering
# _______________________________________________________
for(i in names(dendrograms)){
  png(
    file = file.path(data_paper,
                     "Graphs",
                     paste0(
                       "Dendrogram_single_linkage_", i, ".png"
                     )),
    width = 1100,
    height = 1300,
    res = 100
  )
  as.dendrogram(hc_output_ww[[i]]$hc_single) %>% 
    hang.dendrogram(., hang_height = 0.001) %>% 
    plot(., horiz = TRUE)
  dev.off()
}

# prepare outlier dataset for quick overview
outliers <- rbindlist(list(
  # AUS: "Heteroceridae", "Limnichidae"
  "AUS" = trait_AUS[family %in% c("Heteroceridae",
                                  "Limnichidae"), ],
  # EUR: "Chrysomelidae", "Ephemeridae"
  "EU" = trait_EU[family %in% c("Chrysomelidae",
                                "Ephemeridae",
                                "Baetidae",
                                "Corixidae"),],
  # NOA: "Dixidae", "Sisyridae"
  "NOA" = trait_NOA[family %in% c("Dixidae",
                                  "Sisyridae"),],
  # NZ
  "NZ" = trait_NZ[family %in% c("Stratiomyidae",
                                "Tanyderidae",
                                "Simuliidae",
                                "Hydrometridae"), ]
),
idcol = "continent")
saveRDS(outliers,
        file.path(data_cache, "outliers.rds"))

# long format
outliers_lf <- melt(
  outliers,
  id.vars = c("family", "order", "group", "continent"),
  variable.name = "trait",
  value.name = "affinity"
)
outliers_lf[, grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", trait)] 
outliers_lf[, group := factor(group)]
outliers_lf[order(continent), ] %>% 
  .[affinity >= 0.5, ] %>% 
  .[continent == "AUS", ]

p_list <- list()
for(cont in c("AUS", "EU", "NOA", "NZ")) {
  p_list[[cont]] <-
    fun_heatmap_single_cont(data = outliers_lf[continent == cont,])
}
(p_list[["AUS"]] + p_list[["EU"]])/ (p_list[["NOA"]] + p_list[["NZ"]])

# AUS: 
#   Heteroceridae:  streamlined, gatherers, burrowers, terrestrial ovipositon, respiration with plastron and spircale, small size, expressed all three voltinism traits
# Limnichidae: cylindrical, gatherers, swimmers, terrestrial ovipositon, respiration with plastron and spircale, small size, expressed all three voltinism traits (part of AUS_TPG5)
# EU: 
#   Chrysomelidae: cylindrical, shredder, sessil & burrow, ovip_aqu, resp_teg & resp_pls_spi, size_medium, volt_uni
# Ephemeridae: cylindrical, filterer, burrower, ovip_aqu, resp_teg & resp_gil, size_large, volt_semi
# Baetidae: bf_streamlined & cylindrical, locom_swim & locom_crawl, ovip_aqu, resp_teg & resp_gil, size_small, volt_bi_multi
# Corixidae: flattened & cylindrical, gatherer, swimmer, aquatic oviposition, resp_teg & resp_pls_spi, size small, volt_uni & volt_bi_multi 
# NOA: 
#   Dixidae: cylindrical, gatherer, swimmer, ovip_ter, resp_pls_spi, size_medium & size_small, volt_uni
# Sisyridae: streamlined, predator & herbivore, swimming & crawling, terrestrial oviposition, , resp_teg & resp_gil, size_small, volt_bi_multi & volt_semi (part of NA_TPG 3)
# NZ: 
#   Stratiomyidae: flattened & cylindrical, herbivore, crawling, ovip_ter & ovip_aqu, resp_pls_spi, size_medium, volt_semi 
# Hydrometridae: spherical & cylindrical, predator, crawling, ovip_aqu & ovip_ter, resp_pls_spi, size_medium & size_small, volt_bi_multi
# Simuliidae: cylindrical, herbivore & filterer, sessil, aquatic oviposition, resp_teg & resp_gil, size_small, volt_bi_multi
# Tanyderidae




