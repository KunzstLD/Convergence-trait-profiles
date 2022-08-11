# __________________________________________________________________________________________________
# Analysis of TPGs ----
# Identify traits that characterize the groups obtained through cluster analysis
# Clustering + Ordination  (principal coordinate analysis) together?
# Consider correlations among traits (see Wilkes Paper!)
# __________________________________________________________________________________________________

# Read in trait data with group assignment
trait_CONT <- readRDS(file.path(data_cache, "trait_dat_grp_assig.rds"))

# __________________________________________________________________________________________________
# TPGs ----
# __________________________________________________________________________________________________

# Unique families per continent and group
trait_CONT[, nr_families := uniqueN(family), by = .(continent, group)]

# Delineated TPGs per continent
unique(trait_CONT[order(group), .(group), by = continent])
trait_CONT[, uniqueN(family), by = continent]

# Scaled to 100
trait_CONT[, uniqueN(family), by = continent] %>%
  .[, .(continent, 100 / V1 * c(8, 9, 8, 7, 10))]

## Do TPGs represent a taxonomic signal? ----
trait_CONT_wf <-
  dcast(trait_CONT[, .(continent, family, order, group, trait, affinity)],
    ... ~ trait,
    value.var = "affinity"
  )
trait_CONT_wf[, n_families_gr := .N, by = .(continent, group)]
trait_CONT_wf[, n_families_gr_order := .N, by = .(continent, group, order)]
trait_CONT_wf[, prop_order := n_families_gr_order / n_families_gr]

### Cluster sizes (nr. of families per TPG) ----
cluster_size <- unique(trait_CONT_wf[, .(continent, group, n_families_gr)])
cluster_size[order(continent, -n_families_gr), ]

### In how many clusters does a certain order occur? ----
cluster_occ <- unique(trait_CONT_wf[, .(group), by = .(continent, order)]) %>%
  dcast(., ... ~ group, fun.aggregate = fun_binary_length) %>%
  .[, occ_in_cluster := apply(.SD, 1, sum),
    .SDcols = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
  ] %>%
  .[, .(continent, order, occ_in_cluster)] 
# cluster_occ[order %in% c("Ephemeroptera", "Trichoptera", "Diptera", "Plecoptera"), ]
cluster_occ[order(-occ_in_cluster), ]

# Plotting
tpg_names <- c(
  "1" = "TPG 1",
  "2" = "TPG 2",
  "3" = "TPG 3",
  "4" = "TPG 4",
  "5" = "TPG 5",
  "6" = "TPG 6",
  "7" = "TPG 7",
  "8" = "TPG 8",
  "9" = "TPG 9",
  "10" = "TPG 10",
  "AUS" = "AUS",
  "EU" = "EUR",
  "NOA" = "NA",
  "NZ" = "NZ",
  "SA" = "SA"
)
ggplot(trait_CONT_wf,
       aes(x = as.factor(order),
           y = prop_order * 100)) +
  geom_pointrange(aes(
    ymin = 0,
    ymax = prop_order * 100,
    color = as.factor(continent)
  )) +
  geom_text(data = cluster_size,
            mapping = aes(
              x = 1,
              y = 80,
              label = paste0("n = ", n_families_gr)
            ), 
            size = 5) +
  facet_grid(as.factor(continent) ~ as.factor(group),
             labeller = as_labeller(tpg_names)) +
  coord_flip() +
  labs(x = "Order",
       y = "Proportion of families that belong to a certain order",
       color = "Region") +
  scale_color_d3() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(
      family = "Roboto Mono",
      size = 12,
      angle = 45,
      vjust = 0.6 # ,
      #    hjust = 0.7
    ),
    axis.text.y = element_text(family = "Roboto Mono",
                               size = 12),
    legend.title = element_text(family = "Roboto Mono",
                                size = 16),
    legend.text = element_text(family = "Roboto Mono",
                               size = 14),
    strip.text = element_text(family = "Roboto Mono",
                              size = 14),
    # panel.grid = element_blank(),
    legend.position = "none"
  )
ggplot2::ggsave(
    filename = file.path(data_paper,
                         "Graphs",
                         "Tax_signal_summary.png"),
    width = 42,
    height = 30,
    units = "cm",
    dpi = 600
  )

## Defining traits ----
# Criterion:
# Calculate proportion of taxa per continent, TPG and trait.
# Traits that are expressed with
# an affinity of more than 0.5 per TPG should presumably distinguish the TPG's
trait_CONT[affinity > 0.5,
           prop_taxa_high_aff := .N / nr_families,
           by = .(continent, group, trait)]

# Table for SI
# saveRDS(trait_CONT[affinity > 0.5, ],
#         file.path(data_cache, "def_traits_full.rds"))

# Get defining traits
disting_traits <-
  unique(trait_CONT[prop_taxa_high_aff >= 0.5, ] %>%
    .[order(-prop_taxa_high_aff), .(
      grouping_feature,
      prop_taxa_high_aff,
      nr_families
    ), by = c("continent", "group", "trait")])
disting_traits <-
  disting_traits[order(continent, group, -prop_taxa_high_aff), ]
# View(disting_traits)

# Dataset with traits and prop of taxa with high aff
disting_traits_compl <- copy(disting_traits)
disting_traits_compl[, tpg_id := paste0(continent, "_", group)]
disting_traits_compl[, .N, by = trait] %>% 
.[order(-N), ]

# Calculate nr of traits per tpg for upcoming comparisons
# Overall, TPGs have between 2 and 7 defining traits
disting_traits_compl[, nr_traits_group := .N, 
               by = .(continent, group)]
unique(disting_traits_compl[, nr_traits_group,
  by = c("continent", "group")
])
disting_traits_compl[, .(
  min_def_traits = min(nr_traits_group),
  max_def_traits = max(nr_traits_group)
), by = continent]

# Which traits did not define any TPG?
setdiff(trait_CONT$trait, disting_traits_compl$trait)

### Overview ----

# Dataset with defining trait combinations (TC)
disting_traits[, def_traits := paste(trait, collapse = ", "),
  by = c("continent", "group")
]

# Check for similar TPGs
lookup_traits <- data.table(
  trait = unique(disting_traits$trait),
  trait_label = c(
    "swimming",
    "small",
    "plastron and spiracle",
    "predator",
    "cylindrical",
    "large",
    "univoltinism",
    "crawling",
    "gills",
    "flattened",
    "herbivore",
    "medium",
    "tegument",
    "burrower",
    "streamlined",
    "semivoltinism",
    "shredder",
    "gatherer",
    "bi/multivoltinism",
    "filterer",
    "sessil"
  ))
# saveRDS(lookup_traits,
#   file = "/home/kunzst/Dokumente/Projects/Trait_DB/Convergence-trait-profiles/Cache/lookup_traits.rds"
# )
disting_traits[lookup_traits, 
           trait_label := i.trait_label, 
           on = "trait"]
disting_traits <- disting_traits %>%
  arrange(factor(
    grouping_feature,
    levels = c(
      "locom",
      "bf",
      "ovip",
      "size",
      "feed",
      "resp",
      "volt"
    )
  )) %>%
  .[, paste(trait_label, collapse = ", "),
    by = c("continent", "group")
  ]
setnames(disting_traits,
         "V1",
         "defining_tc")

# Add id
disting_traits[, tpg_id := paste0(continent, "_", group)]

# Save
# saveRDS(disting_traits,
#         file = file.path(data_cache, "def_traits.rds"))

# Exact duplicates 
disting_traits[duplicated(defining_tc) | duplicated(defining_tc, fromLast = TRUE), ]

# Fuzzy matching
tab_fmatch_tc <- disting_traits[, .(tidy_comb_all(defining_tc), continent)] %>%
  .[, .(continent, V1, V2, lv_dist = stringdist(V1, V2, method = "lv"))]
setDT(tab_fmatch_tc)

# Go through fuzzy matching table manually and select potential similar tc's
# tab_fmatch_tc[order(lv_dist), head(V2, n = 5), by = V1] %>% View()

# How many "similar" matches are there?
tab_fmatch_tc[, .N, by = V1] %>% 
.[order(-N), ]
tab_fmatch_tc[, .N, by = c("V1", "continent")] %>% 
.[order(-N)] 

# On how many continents does a certain TC
# and "similar" TCs occur? 
tab_fmatch_tc[, .(occur_continent = uniqueN(continent)),
  by = V1
] %>% 
.[order(-occur_continent), ]

## Explore different trait combinations ----

### Trait combinations that occur on all continents & regions ---- 
# - Crawling, clyindrical, univoltnism  
# occurs as combination on all continents
# other traits that occur often with these TPGs: gills, medium and small size, predator  
# 13 out of 42
disting_traits[defining_tc %like% "crawling" & defining_tc %like% "cylindrical" & defining_tc %like% "univoltinism", ] 
global_pat_id <- disting_traits[defining_tc %like% "crawling" & defining_tc %like% "cylindrical" & defining_tc %like% "univoltinism", tpg_id]
disting_traits_compl[tpg_id %in% global_pat_id, .(continent, .N),
  by = trait
] %>% 
.[order(-N), ]

disting_traits[(defining_tc %like% "crawling" & defining_tc %like% "cylindrical" & 
defining_tc %like% "univoltinism")| (defining_tc %like% "gills" | defining_tc %like% "predator"), ]

# - predator, univoltinism, crawling
# five of the former TPGs included
# Present on all continents and regions
disting_traits[defining_tc %like% "gill" & defining_tc %like% "univoltinism" & defining_tc %like% "crawling", ] %>% 
.[tpg_id %in% global_pat_id,]

disting_traits[defining_tc %like% "predator" & defining_tc %like% "univoltinism" & defining_tc %like% "crawling", ] %>% 
  .[tpg_id %in% global_pat_id,]


# - Crawling, small, univoltinism
# More or less a subset of the former (6 out of the 10 TPGs belong to the former global TC)
disting_traits[defining_tc %like% "crawling" & defining_tc %like% "small" & defining_tc %like% "univoltinism", ] %>% 
.[!tpg_id %in% global_pat_id, ]

# Small & plastron with spiracle 
disting_traits[defining_tc %like% "small" & defining_tc %like% "plastron.*", ] %>%
  .[order(continent), ]
global_pat_id2 <- disting_traits[defining_tc %like% "small" & defining_tc %like% "plastron.*", ] %>%
  .[order(continent), tpg_id]
setdiff(global_pat_id2, global_pat_id)

disting_traits[defining_tc %like% "small" & defining_tc %like% "plastron.*", ] %>%
  .[defining_tc %like% "predator", ]
disting_traits[defining_tc %like% "small" & defining_tc %like% "plastron.*", ] %>%
  .[defining_tc %like% "herbivore", ]
disting_traits[defining_tc %like% "small" & defining_tc %like% "plastron.*", ] %>%
  .[order(continent), ] %>%
  .[!tpg_id %in% global_pat_id, ]


### On 4 continents/regions ----
# - cylindrical, small, plastron and spiracle
# Most are univoltine, some predators, some herbivores
# (only 3 are part of the global TC)
# Subset of small & plastron with spiracle
disting_traits[defining_tc %like% "cylindrical" & defining_tc %like% "small" & defining_tc %like% "plastron.*", ] 
disting_traits[defining_tc %like% "cylindrical" & defining_tc %like% "small" & defining_tc %like% "plastron.*", ] %>% 
.[!tpg_id %in% global_pat_id, ]


# - Flattened, crawling, small
# With gills on three continents and regions
# Only one TPG in global_pat_id
disting_traits[defining_tc %like% "flattened", ] %>% 
.[!tpg_id %in% c(global_pat_id, global_pat_id2), ]


### On 3 regions/continents ----

# - Crawling, Herbivore, univoltinism
# Most TPGs are part of the other TPGs above
disting_traits[defining_tc %like% "crawling" & defining_tc %like% "herbivore" & defining_tc %like% "univoltinism", ]
disting_traits[defining_tc %like% "crawling" & defining_tc %like% "herbivore" & defining_tc %like% "univoltinism", ] %>% 
.[!tpg_id %in% c(global_pat_id, global_pat_id2), ]

# - Herbivore, small, gills
disting_traits[defining_tc %like% "herbivore" & defining_tc %like% "small" & defining_tc %like% "gills", ] %>% 
.[!tpg_id %in% c(global_pat_id, global_pat_id2), ]

# - crawling, medium, gills
disting_traits[defining_tc %like% "crawling" & defining_tc %like% "medium" & defining_tc %like% "gills", ] %>% 
.[!tpg_id %in% c(global_pat_id, global_pat_id2), ]

# Crawling & large, mostly predators
disting_traits[defining_tc %like% "large", ]
trait_CONT[continent == "SA" & group == 5 & trait == "feed_predator",]

### Between Australia, South Afria, and New Zealand ----
# Specific TPGs similar in AUS, NZ & SA (more recent geological history)
disting_traits[continent %in% c("AUS", "NZ", "SA"), ] %>% 
  .[order(continent), ]


### On 2 regions/continents ----

# ?Tegument, univoltnism, cylindrical, small
disting_traits[defining_tc %like% "tegument", ] %>%
  .[!tpg_id %in% c(global_pat_id, global_pat_id2), ]

# Gatherer
# gills, crawling, small
# NOA & SA
disting_traits[defining_tc %like% "gatherer", ] %>%
  .[!tpg_id %in% c(global_pat_id, global_pat_id2), ]

# - Swimming
# On two continents, but totally different in the other defining traits 
disting_traits[defining_tc %like% "swimming", ] %>%
  .[!(tpg_id %in% global_pat_id | tpg_id %in% global_pat_id2), ]


### Traits that only occur once in a TPG/special TPGs ----

# Shredder only in NOA
# with crawling, cylindrical, small, univoltinism
disting_traits[defining_tc %like% "shredder", ] %>% 
.[!(tpg_id %in% global_pat_id | tpg_id %in% global_pat_id2), ]

# Streamlined only in NOA
# with burrower, streamlined, medium, gills, univoltinism
disting_traits[defining_tc %like% "streamlined", ] %>%
  .[!(tpg_id %in% global_pat_id | tpg_id %in% global_pat_id2), ]

# Filterer only in SA
# with crawling, flattened, cylindrical, small, gatherer, gills
disting_traits[defining_tc %like% "filterer", ] %>%
  .[!(tpg_id %in% global_pat_id | tpg_id %in% global_pat_id2), ]

# Sessil only in SA
# with cylindrical, small, herbivore, gills, univoltinism
disting_traits[defining_tc %like% "sessil", ] %>%
  .[!(tpg_id %in% global_pat_id | tpg_id %in% global_pat_id2), ]

# Burrower
disting_traits[defining_tc %like% "burrow", ] %>%
  .[!(tpg_id %in% global_pat_id | tpg_id %in% global_pat_id2), ]

# Bi/multivoltinism
# NOA and NZ
# with plastron and spiracle, cylindrical body form
disting_traits[defining_tc %like% "bi/multivoltinism", ] %>% 
  .[!(tpg_id %in% global_pat_id | tpg_id %in% global_pat_id2), ]

# semivoltinism
disting_traits[defining_tc %like% "semivoltinism", ] %>%
  .[!(tpg_id %in% global_pat_id | tpg_id %in% global_pat_id2), ]



#*****************************************************************
## Crawling, clyindrical body form & univoltnism ----
# Basically this trait combination with
# variations occurs on all continents
# TODO: adapt code
#*****************************************************************
# tpg1_traits <- c(
#   "bf_cylindrical",
#   "ovip_aqu",
#   "locom_crawl", 
#   "size_small",
#   "size_medium",
#   "volt_uni"
# )
# tpg_acr_cont <-
#   trait_CONT[trait %in% tpg1_traits & (
#     (continent == "AUS" & group %in% c(1, 2, 3)) |
#       (continent == "EU" & group %in% c(3, 4, 5, 6, 7)) |
#       (continent == "NOA" & group %in% c(3, 4, 5, 8)) |
#       (continent == "NZ" & group %in% c(3, 5, 6, 10))
#   ), ]
# # tpg_acr_cont$continent %>% unique
# 
# # change order of continent column for plotting
# tpg_acr_cont[, continent := paste0(continent, "_", group)]
# tpg_acr_cont$continent %>% unique()
# tpg_acr_cont[, continent := factor(continent,
#   levels = c(
#     "AUS_1",
#     "AUS_2",
#     "AUS_3",
#     "EU_3",
#     "EU_4",
#     "EU_5",
#     "EU_6",
#     "EU_7",
#     "NOA_3",
#     "NOA_4", 
#     "NOA_5",
#     "NOA_8",
#     "NZ_3",
#     "NZ_5",
#     "NZ_6",
#     "NZ_10"
#   )
# )]
# 
# # Heatmap
# tpg_names1 <- c(
#   "AUS_1" = "AUS_TPG1",
#   "AUS_2" = "AUS_TPG2",
#   "AUS_3" = "AUS_TPG3",
#   "EU_3" = "EU_TPG3",
#   "EU_4" = "EU_TPG4",
#   "EU_5" = "EU5_TPG",
#   "EU_6" = "EU6_TPG",
#   "EU_7" = "EU7_TPG",
#   "NOA_3" = "NA_TPG3",
#   "NOA_4" = "NA_TPG4",
#   "NOA_5" = "NA_TPG5",
#   "NOA_8" = "NA_TPG8",
#   "NZ_3" = "NZ_TPG3",
#   "NZ_5" = "NZ_TPG5",
#   "NZ_6" = "NZ_TPG6",
#   "NZ_10" = "NZ_TPG10"
# )
# plot_list <- list()
# for(cont in c("AUS", "EU", "NOA", "NZ")) {
#   plot_list[[cont]] <- fun_heatmap_tpg(data = tpg_acr_cont[continent %like% cont,],
#                                        facet_names = tpg_names1)
# }
# (plot_list[["AUS"]] + plot_list[["EU"]]) / (plot_list[["NOA"]] + plot_list[["NZ"]])
# ggsave(
#   filename = file.path(
#     data_paper,
#     "Graphs",
#     "Heatmap_TPG_across_all_continents_50.png"
#   ),
#   width = 60,
#   height = 30,
#   units = "cm"
# )
# 
# #*****************************************************************
# ## Predators, small size, respiration with plastron and spiracle ----
# #*****************************************************************
# tpg2_traits <- c(
#   "feed_predator",
#   "size_small",
#   "resp_pls_spi",
#   "ovip_aqu",
#   "volt_uni",
#   "locom_swim"
# )
# 
# # Similar TPGs:
# # AUS_TPG5, EU_TPG1, EU_TPG2,  NA_TPG2,  NZ_TPG5
# tpg2_acr_cont <-
#   trait_CONT[trait %in% tpg2_traits & (
#     (continent == "AUS" & group == 5) |
#       (continent == "EU" & group %in% c(1, 2)) |
#       (continent == "NOA" & group == 2) |
#       (continent == "NZ" & group %in% 6)
#   ), ]
# 
# # change order of continent column for plotting
# tpg2_acr_cont[, continent := paste0(continent, "_", group)]
# tpg2_acr_cont$continent %>% unique()
# tpg2_acr_cont[, continent := factor(continent,
#                                     levels = c("AUS_5",
#                                                "EU_1",
#                                                "EU_2",
#                                                "NOA_2",
#                                                "NZ_6"))]
# tpg2_acr_cont[, trait := factor(
#   trait,
#   levels = c(
#     "volt_uni",
#     "ovip_aqu",
#     "locom_swim",
#     "feed_predator",
#     "size_small",
#     "resp_pls_spi"
#   )
# )]
# 
# # Heatmap
# tpg2_names <- c(
#   "AUS_5" = "AUS_TPG5",
#   "EU_1" = "EU_TPG1",
#   "EU_2" = "EU_TPG2",
#   "NOA_2" = "NA_TPG2",
#   "NZ_6" = "NZ_TPG6")
# plot_list <- list()
# for(cont in c("AUS", "EU", "NOA", "NZ")) {
#   plot_list[[cont]] <- fun_heatmap_tpg(data = tpg2_acr_cont[continent %like% cont,],
#                                        facet_names = tpg2_names)
# }
# (plot_list[["AUS"]] + plot_list[["EU"]]) / (plot_list[["NOA"]] + plot_list[["NZ"]])
# ggsave(
#   filename = file.path(
#     data_paper,
#     "Graphs",
#     "Heatmap_TPG_across_all_continents_40.png"
#   ),
#   width = 55,
#   height = 35,
#   units = "cm"
# )