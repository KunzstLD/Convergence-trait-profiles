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
trait_CONT[, .(uniqueN(family), uniqueN(group)), by = continent] %>%
  .[, .(continent, 100 / V1 * V2)]

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
cluster_occ <-
  unique(trait_CONT_wf[, .(group), by = .(continent, order)]) %>%
  dcast(., ... ~ group, fun.aggregate = fun_binary_length) %>%
  .[, occ_in_cluster := apply(.SD, 1, sum),
    .SDcols = c("1",
                "2",
                "3",
                "4",
                "5",
                "6",
                "7",
                "8",
                "9",
                "10")] %>%
  .[, .(continent, order, occ_in_cluster)] 
# cluster_occ[order %in% c("Ephemeroptera", "Trichoptera", "Diptera", "Plecoptera"), ]
cluster_occ[order(-occ_in_cluster), ]

# Plotting
# tpg_names <- c(
#   "1" = "TPG 1",
#   "2" = "TPG 2",
#   "3" = "TPG 3",
#   "4" = "TPG 4",
#   "5" = "TPG 5",
#   "6" = "TPG 6",
#   "7" = "TPG 7",
#   "8" = "TPG 8",
#   "9" = "TPG 9",
#   "10" = "TPG 10",
#   "AUS" = "AUS",
#   "EU" = "EUR",
#   "NOA" = "NA",
#   "NZ" = "NZ",
#   "SA" = "SA"
# )
# ggplot(trait_CONT_wf,
#        aes(x = as.factor(group),
#            y = prop_order * 100)) +
#   geom_col() +
#   # geom_pointrange(aes(
#   #   ymin = 0,
#   #   ymax = prop_order * 100,
#   #   color = as.factor(continent)
#   # )) +
#   geom_text(data = cluster_size,
#             mapping = aes(
#               x = 1,
#               y = 80,
#               label = paste0("n = ", n_families_gr)
#             ), 
#             size = 5) +
#   facet_grid(order ~ as.factor(continent),
#              labeller = as_labeller(tpg_names),
#              scales = "free") +
#   coord_flip() +
#   labs(x = "Order",
#        y = "Proportion of families that belong to a certain order",
#        color = "Region") +
#   scale_color_d3() +
#   theme_bw() +
#   theme(
#     axis.title = element_text(size = 16),
#     axis.text.x = element_text(
#       family = "Roboto Mono",
#       size = 12,
#       angle = 45,
#       vjust = 0.6 # ,
#       #    hjust = 0.7
#     ),
#     axis.text.y = element_text(family = "Roboto Mono",
#                                size = 12),
#     legend.title = element_text(family = "Roboto Mono",
#                                 size = 16),
#     legend.text = element_text(family = "Roboto Mono",
#                                size = 14),
#     strip.text = element_text(family = "Roboto Mono",
#                               size = 14),
#     # panel.grid = element_blank(),
#     legend.position = "none"
#   )
# ggplot2::ggsave(
#     filename = file.path(data_paper,
#                          "Graphs",
#                          "Tax_signal_summary.png"),
#     width = 42,
#     height = 30,
#     units = "cm",
#     dpi = 600
#   )

# New attempt, difficult to visualize since orders are on the x-axis
# better in a table?
# unique(trait_CONT_wf[, .(continent, group, order, prop_order)]) %>%
#  .[group %in% c(1,2) & continent %in% c("AUS", "EU"), ] %>% 
#   ggplot(.,
#          aes(x = as.factor(order),
#              y = prop_order*100)) +
#   geom_col() +
#   facet_grid(as.factor(continent) ~ group)


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
# volt_semi, locom_burrow, bf_spherical
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
    "large",
    "crawling",
    "predator",
    "cylindrical",
    "gills",
    "univoltinism",
    "small",
    "herbivore",
    "medium", 
    "plastron and spiracle",
    "swimming",
    "bi/multivoltinism",
    "tegument",
    "gatherer",
    "flattened",
    "shredder",
    "streamlined",
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
# 1) crawling, cylindrical, univoltinism, ...& related
# 2) ?
tab_fmatch_tc[, .(occur_continent = uniqueN(continent)),
  by = V1
] %>% 
.[order(-occur_continent), ]

# Consider distances between defining traits (in characters)
tab_fmatch_tc[order(lv_dist), ]

tab_fmatch_tc[, .(occur_continent = uniqueN(continent), lv_dist),
              by = V1
] %>% 
  .[order(-occur_continent, lv_dist), ] %>%
  View(.)


## Explore different trait combinations ----

### Trait combinations that occur on all continents & regions ---- 
# - Crawling, clyindrical, gills, univoltinism  
# occurs as combination on all continents
# other traits that occur often with these TPGs: medium and small size, predator  
# 10 times
disting_traits[defining_tc %like% "crawling" &
                 defining_tc %like% "cylindrical" &
                 defining_tc %like% "gills",] 

# 13 times (out of 42)
disting_traits[defining_tc %like% "crawling" &
                 defining_tc %like% "cylindrical" &
                 defining_tc %like% "univoltinism",] 

# Also on every continent/region
disting_traits[defining_tc %like% "crawling" &
                 defining_tc %like% "cylindrical" &
                 defining_tc %like% "gills" & 
                 defining_tc %like% "univoltinism",] 

global_pat_id <-
  disting_traits[defining_tc %like% "crawling" &
                   defining_tc %like% "cylindrical" &
                   defining_tc %like% "univoltinism", tpg_id]

# traits that frequently occur with the combination crawling+cylindrical+univoltinism
disting_traits_compl[tpg_id %in% global_pat_id, .(continent, .N),
                     by = trait] %>%
  .[order(-N),]
disting_traits[defining_tc %like% "crawling" &
                 defining_tc %like% "cylindrical" &
                 defining_tc %like% "small",] 

# Crawling & cylindrical occur in more than half of the TPGs as defining trait combination (23)
disting_traits[(defining_tc %like% "crawling" & defining_tc %like% "cylindrical"),]

# - predator & plastron with spiracle

# predator & plastron and spiracle
# frequently associated small, cylindrical, and swimming
disting_traits[defining_tc %like% "predator" &
                 defining_tc %like% "plastron.*", ] %>%
  .[order(continent), ] %>% 
  .[!tpg_id %in% global_pat_id, ]


global_pat_id2 <-
  disting_traits[defining_tc %like% "predator" &
                   defining_tc %like% "plastron.*",] %>%
  .[order(continent), tpg_id]
# Completely different TPGs compared to the first group 
# setdiff(global_pat_id2, global_pat_id)

### On 4 continents/regions ----

# - herbivore, cylindrical, small
disting_traits[defining_tc %like% "herbivore" &
                 defining_tc %like% "cylindrical" & 
                 defining_tc %like% "small", ]

# # - crawling, medium, gills
disting_traits[defining_tc %like% "crawling" & defining_tc %like% "medium" & 
                 defining_tc %like% "gills", ] # %>% 
  # .[!tpg_id %in% c(global_pat_id, global_pat_id2), ]

### On 3 regions/continents ----

# crawling, cylindrical, predator, tegument, univoltinism
disting_traits[defining_tc %like% "tegument", ]

# Bi/multivoltinism
# 3 TPGs AUS, NOA, NZ
# with plastron and spiracle, small
disting_traits[defining_tc %like% "bi/multivoltinism", ] # %>%
#   .[!(tpg_id %in% global_pat_id | tpg_id %in% global_pat_id2), ]

# Tegument, univoltnism, cylindrical, small
# EU, NZ, NOA
disting_traits[defining_tc %like% "tegument", ] 


# 
# # Crawling & large, mostly predators
# disting_traits[defining_tc %like% "large", ]
# trait_CONT[continent == "SA" & group == 5 & trait == "feed_predator",]
# 
# ### Between Australia, South Afria, and New Zealand ----
# Specific TPGs similar in AUS, NZ & SA (more recent geological history)
# Not really something specific
# disting_traits[continent %in% c("AUS", "NZ", "SA"), ] %>%
#   .[order(continent), ]


### On 2 regions/continents ----
disting_traits[defining_tc %like% "shredder",] %>%
  .[!tpg_id %in% global_pat_id,]

# Gatherer
## gills, crawling, small
## EU & SA
disting_traits[defining_tc %like% "gatherer", ] %>%
  .[!tpg_id %in% c(global_pat_id, global_pat_id2), ]
 
# - Shredder in NOA & SA
# # with crawling, cylindrical, small, univoltinism
disting_traits[defining_tc %like% "shredder", ]

# Filterer in NOA & NZ
# # with crawling, flattened, cylindrical, small, gatherer, gills
disting_traits[defining_tc %like% "filterer", ] #%>%
#   .[!(tpg_id %in% global_pat_id | tpg_id %in% global_pat_id2), ]

# Sessil in SA & NZ
# with cylindircal, apart from that TPG are quite different
disting_traits[defining_tc %like% "sessil", ] # %>%
#   .[!(tpg_id %in% global_pat_id | tpg_id %in% global_pat_id2), ]

# - large
disting_traits[defining_tc %like% "large", ]
 
#### Traits that only occur once in a TPG/special TPGs ----
# streamlined and ?  

# Streamlined only in NOA
# with filterer, gills, univoltinism
disting_traits[defining_tc %like% "streamlined", ] #%>%
#   .[!(tpg_id %in% global_pat_id | tpg_id %in% global_pat_id2), ]


# TPGs with similar defining trait combinations ----
# NZ 3 & 4 Actually identical
trait_CONT[continent == "NZ" & group %in% c(3,4), ] %>%
  .[prop_taxa_high_aff >= 0.2, ] %>% 
  .[order(group),]
  .[order(group), uniqueN(family), by = c("group", "trait")]