# __________________________________________________________________________________________________
#### Analysis of TPGs ####
# Identify traits that characterize the groups obtained through cluster analysis
# Clustering + Ordination  (principal coordinate analysis) together? 
# Consider correlations among traits (see Wilkes Paper!)
# __________________________________________________________________________________________________

# Read in trait data with group assignment
trait_CONT <- readRDS(file.path(data_cache, "trait_dat_grp_assig.rds"))

# --- TPG's ----------------------------------------------------------------------------------------

# Unique families per continent and group
trait_CONT[, nr_families := uniqueN(family), by = .(continent, group)]

# Calculate proportion of taxa per continent, TPG and trait. Traits that are expressed with 
# an affinity of more than 0.5 per TPG should presumably distinguish the TPG's (?)
trait_CONT[affinity > 0.5 , prop_taxa_high_aff := .N/nr_families, 
           by = .(continent,group, trait)]
disting_traits <-
  unique(trait_CONT[prop_taxa_high_aff > 0.5, ] %>%
           .[order(-prop_taxa_high_aff), .(continent, 
                                             group,
                                             trait,
                                             grouping_feature,
                                             prop_taxa_high_aff,
                                             nr_families)], by = c("continent", "group", "trait"))
disting_traits <- disting_traits[order(continent, group, -prop_taxa_high_aff), ] 

# Calculate nr of traits per tpg for upcoming comparisons
disting_traits[, nr_traits_group := .N, by = .(continent, group)]

# Are there similarities across continents?
# A) between the clustered groups
output_tpgs <- list()
for(cont in c("AUS", "EU", "NOA", "NZ")) {
  groups <- unique(disting_traits[continent == cont, group])
  similar_tpgs <- list()
  
  for (i in groups) {
    tpg <- disting_traits[continent == cont & group == i, trait]
    
    # length of subsetted tpg (equal to the nr_traits_group for this subset!)
    l_tpg <- length(tpg)
    
    similar_tpgs[[i]] <-
      disting_traits[trait %in% tpg,] %>%
      .[order(continent, group),] %>%
      .[, .(nr_comp = .N,
            trait,
            prop_taxa_high_aff,
            nr_families,
            nr_traits_group),
        by = .(continent, group)] %>%
      .[nr_comp == l_tpg,]
  }
  output_tpgs[[cont]] <- similar_tpgs
}

# 1) select only tpgs where nr_traits_group equals nr_comp, i.e. the whole tpgs are identical 
# across regions, not only parts of it
lapply(output_tpgs, function(y)
  lapply(y, function(z)
    z[length(unique(continent)) == 4,]))

# 1.1) Distinguish low and high affinities?
# 1.2) Pilliere found differences within orders (e.g. Ephemeroptera) 
# -> something worth to check out
disting_traits[trait %in% c("locom_crawl", "ovip_aqu", "bf_cylindrical", "resp_gil"), ] %>% 
  .[order(continent, group), ]

# Get tpgs and plot in single heatmaps
# Add Order information in separate barchart

#*****************************************************************
# TPG 1 that occurs in all continents
#*****************************************************************
tpg1_traits <- c(
  "locom_crawl",
  "ovip_aqu",
  "bf_cylindrical",
  "resp_gil",
  "size_medium",
  "feed_predator",
  "volt_uni"
)
tpg_acr_cont <-
  trait_CONT[trait %in% tpg1_traits & (
    (continent == "AUS" & group %in% c(1, 2)) |
      (continent == "EU" & group == 7) |
      (continent == "NOA" & group == 4) |
      (continent == "NZ" & group == 8)
  ),]
tpg_acr_cont[, continent := paste0(continent, "_", group)]

# change orders of trait and continent column for plotting
tpg_acr_cont[, trait := factor(
  trait,
  levels = c(
    "locom_crawl",
    "ovip_aqu",
    "bf_cylindrical",
    "resp_gil",
    "size_medium",
    "feed_predator",
    "volt_uni"
  )
)]
tpg_acr_cont$continent %>% unique
tpg_acr_cont[, continent := factor(continent,
                                   levels = c("AUS_1",
                                              "AUS_2",
                                              "EU_7",
                                              "NOA_4",
                                              "NZ_8"))]

# Heatmap
conv_gr1 <- fun_heatmap_tpg(data = tpg_acr_cont)
conv_gr1 <- conv_gr1 +
  annotate(
    xmin = 0.5,
    xmax = 4.5,
    ymin = -Inf,
    ymax = Inf,
    geom = "rect",
    alpha = 0.1,
    color = "darkorange",
    size = 2.1
  ) 
  # +
  # ggtitle("TPG A that occurred in all continents")

# Add order information
tpg_acr_cont_uq <- unique(tpg_acr_cont, by = c("family","continent"))
tpg_acr_cont_uq[, prop_order := .N/nr_families, by = .(continent, order)]

orders_conv_gr1 <- ggplot(tpg_acr_cont_uq, aes(x = as.factor(order),
                                               y = prop_order * 100)) +
  geom_pointrange(aes(ymin = 0,
                      ymax = prop_order * 100)) +
  facet_grid(continent ~ .) +
  labs(x = "Order", y = "Percentage per TPG") +
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
conv_gr1 + orders_conv_gr1 + plot_annotation(tag_levels = "I")
ggplot2::ggsave(
  filename = file.path(data_paper,
                       "Graphs", 
                       "Heatmap_TPG1_across_all_continents.png"),
  width = 35,
  height = 31,
  units = "cm"
)

#*****************************************************************
# TPG 2 that occurs in all continents
#*****************************************************************
tpg2_traits <- c(
  "ovip_aqu",
  "volt_uni",
  "feed_predator",
  "bf_cylindrical",
  "locom_crawl",
  "resp_gil",
  "size_medium",
  "locom_swim",
  "resp_pls_spi",
  "size_small"
)
tpg_acr_cont_2 <-
  trait_CONT[trait %in% tpg2_traits & (
    (continent == "AUS" & group == 1) |
      (continent == "EU" & group %in% c(6, 8, 9)) |
      (continent == "NOA" & group == 3) |
      (continent == "NZ" & group == 8)
  ), ]
tpg_acr_cont_2[, continent := paste0(continent, "_", group)]

# change orders of trait and continent column for plotting
tpg_acr_cont_2[, trait := factor(
  trait,
  levels = c(
    "ovip_aqu",
    "volt_uni",
    "feed_predator",
    "bf_cylindrical",
    "locom_crawl",
    "resp_gil",
    "size_medium",
    "locom_swim",
    "resp_pls_spi",
    "size_small"
  )
)]
tpg_acr_cont_2[, continent := factor(continent,
                                     levels = c("AUS_1",
                                                "EU_6",
                                                "EU_8",
                                                "EU_9",
                                                "NOA_3",
                                                "NZ_8"))]

# NZ group 6; AUS group 4 also occur in TPG1 
conv_gr2 <- fun_heatmap_tpg(data = tpg_acr_cont_2)
conv_gr2 <- conv_gr2 +
  annotate(
    xmin = 0.5,
    xmax = 3.5,
    ymin = -Inf,
    ymax = Inf,
    geom = "rect",
    alpha = 0.1,
    color = "darkorange",
    size = 2.1
  ) # +
# ggtitle("TPG B that occurred in all continents")

# Add order information
tpg_acr_cont_uq2 <- unique(tpg_acr_cont_2, by = c("family","continent"))
tpg_acr_cont_uq2[, prop_order := .N/nr_families, by = .(continent, order)]

orders_conv_gr2 <- ggplot(tpg_acr_cont_uq2, aes(x = as.factor(order),
                                               y = prop_order * 100)) +
  geom_pointrange(aes(ymin = 0,
                      ymax = prop_order * 100)) +
  facet_grid(continent ~ .) +
  labs(x = "Order", y = "Percentage per TPG") +
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
conv_gr2 + orders_conv_gr2 + plot_annotation(tag_levels = "I")
ggplot2::ggsave(
  filename = file.path(data_paper,
                       "Graphs", 
                       "Heatmap_TPG2_across_all_continents.png"),
  width = 35,
  height = 31,
  units = "cm"
)

#*****************************************************************
# TPG 3 that occurs in all continents
#*****************************************************************
tpg3_traits <- c(
  "size_small",
  "ovip_aqu",
  "bf_cylindrical",
  "locom_crawl",
  "feed_herbivore",
  "volt_uni",
  "resp_pls_spi"
)

tpg_acr_cont_3 <-
  trait_CONT[trait %in% tpg3_traits & (
    (continent == "AUS" & group == 6) |
      (continent == "EU" & group %in% c(2, 3)) |
      (continent == "NOA" & group == 1) |
      (continent == "NZ" & group %in% c(1, 2, 3))
  ),]
tpg_acr_cont_3[, continent := paste0(continent, "_", group)]

# change orders of trait and continent column for plotting
tpg_acr_cont_3[, trait := factor(
  trait,
  levels = c(
    "size_small",
    "ovip_aqu",
    "bf_cylindrical",
    "locom_crawl",
    "feed_herbivore",
    "volt_uni",
    "resp_pls_spi"
    )
)]
tpg_acr_cont_3[, continent := factor(continent,
                                     levels = c("AUS_6",
                                                "EU_2",
                                                "EU_3",
                                                "NOA_1",
                                                "NZ_1",
                                                "NZ_2",
                                                "NZ_3"))]

# NZ group 6; AUS group 4 also occur in TPG1 
conv_gr3 <- fun_heatmap_tpg(data = tpg_acr_cont_3)
conv_gr3 <- conv_gr3 +
  annotate(
    xmin = 0.5,
    xmax = 3.5,
    ymin = -Inf,
    ymax = Inf,
    geom = "rect",
    alpha = 0.1,
    color = "darkorange",
    size = 2.1
  ) # +
# ggtitle("TPG C that occurred in all continents")

# Add order information
tpg_acr_cont_uq3 <- unique(tpg_acr_cont_3, by = c("family","continent"))
tpg_acr_cont_uq3[, prop_order := .N/nr_families, by = .(continent, order)]

orders_conv_gr3 <- ggplot(tpg_acr_cont_uq3, aes(x = as.factor(order),
                                                y = prop_order * 100)) +
  geom_pointrange(aes(ymin = 0,
                      ymax = prop_order * 100)) +
  facet_grid(continent ~ .) +
  labs(x = "Order", y = "Percentage per TPG") +
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
conv_gr3 + orders_conv_gr3 + plot_annotation(tag_levels = "I")
ggplot2::ggsave(
  filename = file.path(data_paper,
                       "Graphs", 
                       "Heatmap_TPG3_across_all_continents.png"),
  width = 35,
  height = 31,
  units = "cm"
)

# B) In terms of grouping features contributing to the global distance
# What does the value exactly mean?
# From the docs of kdist.cor(): provides the correlations between the (squared) distances obtained
# for each trait and the global (squared) distances obtained by mixing all the traits
# (= contributions of traits to the global distances);
global_dist <- readRDS(file.path(data_cache, "global_dist.rds"))
global_dist <- rbindlist(global_dist, idcol = "continent")


# For word output flextable seems to be the best option
# library(flextable)
# flextable(global_dist[, .(continent, grouping_feature, `global distance`)]) %>%
#   theme_vanilla(.) %>%
#   hline(.) %>%
#   set_header_labels(
#     .,
#     continent = "Continent",
#     grouping_feature = "Grouping feature",
#     `global distance` = "Global distance"
#   ) %>%
#   autofit(.) %>% 
#   set_caption(
#     "Grouping features that contributed most to the global distance in the distance
#                 matrices in terms of correlation between the distances obtained for each grouping
#                 feature and the global distance obtained by mixing all the grouping features"
#   ) %>% 
#   save_as_docx(., path = file.path(data_paper, "Tables", "Tables_convergence_tp.docx"))

# Markdown output
# kable(global_dist[order(continent, -`global distance`), .(continent, 
#                                                           grouping_feature,
#                                                           `global distance`)],
#       format = "markdown")

