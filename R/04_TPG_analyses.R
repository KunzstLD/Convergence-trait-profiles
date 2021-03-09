# __________________________________________________________________________________________________
#### Analysis of TPGs ####
# Identify traits that characterize the groups obtained through cluster analysis
# Clustering + Ordination  (principal coordinate analysis) together? 
# Consider correlations among traits (see Wilkes Paper!)
# __________________________________________________________________________________________________

# Read in trait data with group assignment
trait_CONT <- readRDS(file.path(data_cache, "trait_dat_grp_assig.rds"))

# Rm dev traits for now
trait_CONT <- trait_CONT[!trait %in% c("dev_hemimetabol", "dev_holometabol"),] 


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

# Tables for results part
# library(knitr)
# library(flextable)
# trait_cols <- as.character(unique(disting_traits$trait))
# path_tables <- "./Paper/Tables/Tables_convergence_tp.docx"
# 
# disting_traits_tv <-
#   dcast(disting_traits[, .(continent, group, trait)],
#         continent + group ~ trait,
#         value.var = "trait")
# imp_traits_wf <- disting_traits_tv[, apply(.SD, 1, function(y)
#   paste(y[!is.na(y)], collapse = ", ")), .SDcols = trait_cols]
# 
# # Create table for word 
# dcast(cbind(disting_traits_tv[, .(continent, group)],
#       imp_traits_wf), 
#       group ~ continent, 
#       value.var = "imp_traits_wf") %>%
#   flextable(.) %>%
#   theme_vanilla(.) %>% 
#   hline(.) %>% 
#   set_caption("Traits expressed by the invertebrates clustered in the trait profile groups. Only 
#       traits have been selected that more than 50 % of the families within a TPG express with an trait 
#       affinity of at least 0.75.") %>% 
#   save_as_docx(., path = path_tables)
# 
# # For Appendix
# disting_traits[, perc_taxa_express_05 := ceiling(prop_taxa_high_aff*100)]
# kable(
#   dcast(
#     data = disting_traits,
#     formula = continent + group + nr_families  ~ trait,
#     value.var = "perc_taxa_express_05",
#     fill = NA
#   ),
#   format = "markdown"
# )

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
    (continent == "AUS" & group %in% c(3, 4)) |
      (continent == "EU" & group == 8) |
      (continent == "NOA" & group == 3) |
      (continent == "NZ" & group == 6)
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
                                   levels = c("AUS_3",
                                              "AUS_4",
                                              "EU_8",
                                              "NOA_3",
                                              "NZ_6"))]

# Heatmap
tpg1_plot <- fun_heatmap_tpg(data = tpg_acr_cont)
tpg1_plot +
  annotate(
    xmin = 0.5,
    xmax = 4.5,
    ymin = -Inf,
    ymax = Inf,
    geom = "rect",
    alpha = 0.1,
    color = "darkorange",
    size = 2.1
  ) # +
  # ggtitle("TPG A that occurred in all continents")
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
    (continent == "AUS" & group == 4) |
      (continent == "EU" & group %in% c(2, 5, 7)) |
      (continent == "NOA" & group == 5) |
      (continent == "NZ" & group == 6)
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
tpg_acr_cont_2$continent %>% unique
tpg_acr_cont_2[, continent := factor(continent,
                                     levels = c("AUS_4",
                                                "EU_2",
                                                "EU_5",
                                                "EU_7",
                                                "NOA_5",
                                                "NZ_6"))]

# NZ group 6; AUS group 4 also occur in TPG1 
tpg2_plot <- fun_heatmap_tpg(data = tpg_acr_cont_2)
tpg2_plot +
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
    (continent == "AUS" & group == 1) |
      (continent == "EU" & group %in% c(3, 4)) |
      (continent == "NOA" & group == 1) |
      (continent == "NZ" & group %in% c(3, 7, 9))
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
                                     levels = c("AUS_1",
                                                "EU_3",
                                                "EU_4",
                                                "NOA_1",
                                                "NZ_3",
                                                "NZ_7",
                                                "NZ_9"))]

# NZ group 6; AUS group 4 also occur in TPG1 
tpg3_plot <- fun_heatmap_tpg(data = tpg_acr_cont_3)
tpg3_plot +
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
ggplot2::ggsave(
  filename = file.path(data_paper,
                       "Graphs", 
                       "Heatmap_TPG3_across_all_continents.png"),
  width = 35,
  height = 31,
  units = "cm"
)

# TODO: Order information for the converged groups
# Barchart for information on orders
tpg_acr_cont[continent == "NOA", ] %>% 
  .[!duplicated(family), .N, by = order]

tpg_acr_cont_uq <- unique(tpg_acr_cont, by = c("family","continent"))
tpg_acr_cont_uq[, prop_order := .N/nr_families, by = .(continent, order)]

ggplot(tpg_acr_cont_uq, aes(x = as.factor(order),
                            y = prop_order * 100)) +
  geom_pointrange(aes(ymin = 0,
                      ymax = prop_order * 100)) +
  facet_grid(continent ~ .) +
  labs(x = "Order", y = "Percentage per TPG") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12),
    axis.text.x = element_text(
      family = "Roboto Mono",
      size = 11
    ),
    axis.text.y = element_text(family = "Roboto Mono",
                               size = 11),
    legend.title = element_text(family = "Roboto Mono",
                                size = 12),
    legend.text = element_text(family = "Roboto Mono",
                               size = 11),
    strip.text = element_text(family = "Roboto Mono",
                              size = 11),
    panel.grid = element_blank()
  )

# B) In terms of grouping features contributing to the global distance
# What does the value exactly mean?
# From the docs of kdist.cor(): provides the correlations between the (squared) distances obtained
# for each trait and the global (squared) distances obtained by mixing all the traits
# (= contributions of traits to the global distances);
global_dist <- readRDS(file.path(data_cache, "global_dist.rds"))
global_dist <- rbindlist(global_dist, idcol = "continent")

kable(global_dist[order(continent, -`global distance`), .(continent, 
                                                          grouping_feature,
                                                          `global distance`)],
      format = "markdown")

