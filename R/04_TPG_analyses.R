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
                "10",
                "11",
                "12")] %>%
  .[, .(continent, order, occ_in_cluster)] 
# cluster_occ[order %in% c("Ephemeroptera", "Trichoptera", "Diptera", "Plecoptera"), ]
cluster_occ[order(-occ_in_cluster), ]

# New attempt, difficult to visualize since orders are on the x-axis
# better in a table?
# unique(trait_CONT_wf[, .(continent, group, order, prop_order)]) %>%
#  .[group %in% c(1,2) & continent %in% c("AUS", "EU"), ] %>%
#   ggplot(.,
#          aes(x = as.factor(order),
#              y = prop_order*100)) +
#   geom_col() +
#   facet_grid(as.factor(continent) ~ group)
# ggplot2::ggsave(
#     filename = file.path(data_paper,
#                          "Graphs",
#                          "Tax_signal_summary.png"),
#     width = 42,
#     height = 30,
#     units = "cm",
#     dpi = 600
#   )
prop_order_tpg <- unique(trait_CONT_wf[, .(continent, group, order, prop_order, n_families_gr)]) %>%
  .[, .(continent, group, order, prop_order = round(prop_order, digits = 3)*100, n_families_gr)] %>% 
  dcast(., ... ~ order, value.var = "prop_order")
setnames(prop_order_tpg,
         c("continent", "group", "n_families_gr"),
         c("Continent", "TPG", "Nr. families TPG"))
# fwrite(prop_order_tpg,
#        file.path(data_paper, "Tables", "prop_order_tpg.csv"),
#        sep = ";")

# Distribution of orders
# unique(trait_CONT_wf[, .(continent, group, order, prop_order, n_families_gr)]) %>%
#   .[, .(continent, group, order, prop_order = round(prop_order, digits = 3)*100, n_families_gr)] %>% 
#   .[prop_order > 0, .N, by = c("continent", "order")] %>% 
#   .[order(continent, -N), ]


## Defining traits ----
# Criterion:
# Calculate proportion of taxa per continent, TPG and trait.
# Traits that are expressed with
# an affinity of more than 0.5 per TPG should presumably distinguish the TPG's
trait_CONT[affinity > 0.5,
           prop_taxa_high_aff := .N / nr_families,
           by = .(continent, group, trait)]
# saveRDS(trait_CONT, file.path(data_cache, "trait_CONT_prop_taxa_high.rds"))
# Table for SI
# saveRDS(trait_CONT[affinity > 0.5, ],
#         file.path(data_cache, "def_traits_full.rds"))

# Example families
# unique(trait_CONT_wf[, .(continent, group, family, order, prop_order, n_families_gr)]) %>% 
#   .[order(continent, group, -prop_order, order), ] %>% 
#   .[continent == "SA" & group %in% c(7,8,9,10), ]

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
# saveRDS(disting_traits_compl, file.path(data_cache, "disting_traits_compl.rds"))

# Which traits did not define any TPG?
# locom_burrow, bf_spherical
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
    "shredder",
    "tegument",
    "sessil",
    "semivoltinism",
    "flattened",
    "streamlined",
    "filterer",
    "gatherer"
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

# Save, table in paper
# saveRDS(disting_traits,
#         file = file.path(data_cache, "def_traits.rds"))

# Example families
trait_CONT[, .SD, by = c("continent", "group")] %>% 
  .[order(continent, group), ]

# Exact duplicates 
disting_traits[duplicated(defining_tc) | duplicated(defining_tc, fromLast = TRUE), ]
# trait_CONT[continent == "NZ" & group %in% c(3, 4) & affinity > 0.5, ] %>% 
#   .[order(group), ]

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
# tab_fmatch_tc[, .(occur_continent = uniqueN(continent), lv_dist),
#               by = V1
# ] %>% 
#   .[order(-occur_continent, lv_dist), ] %>%
#   View(.)


## Explore different trait combinations ----

### Trait combinations that occur on all continents & regions ---- 
# - Crawling, clyindrical, gills, univoltinism  
# occurs as combination on all continents
# other traits that occur often with these TPGs: medium and small size, predator  

# 11 times (out of 39)
disting_traits[defining_tc %like% "crawling" &
                 defining_tc %like% "cylindrical" &
                 defining_tc %like% "univoltinism", ] 

# 10 times
disting_traits[defining_tc %like% "crawling" &
                 defining_tc %like% "cylindrical" &
                 defining_tc %like% "small",] 

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

disting_traits[defining_tc %like% "small" &
                 defining_tc %like% "plastron.*", ]

global_pat_id2 <-
  disting_traits[defining_tc %like% "predator" &
                   defining_tc %like% "plastron.*",] %>%
  .[order(continent), tpg_id]
# Completely different TPGs compared to the first group 
# setdiff(global_pat_id2, global_pat_id)

### On 4 continents/regions ----

# - crawling, cylindrical, gills, univoltinism
disting_traits[defining_tc %like% "crawling" &
                 defining_tc %like% "cylindrical" &
                 defining_tc %like% "gills" & 
                 defining_tc %like% "univoltinism",] 

# - herbivore, cylindrical, small
disting_traits[defining_tc %like% "herbivore" &
                 defining_tc %like% "cylindrical" & 
                 defining_tc %like% "small", ]

### On 3 regions/continents ----

# crawling, cylindrical, predator, tegument, univoltinism
disting_traits[defining_tc %like% "tegument", ]

# - crawling, medium, gills
disting_traits[defining_tc %like% "crawling" & defining_tc %like% "medium" & 
                 defining_tc %like% "gills", ] # %>% 
# .[!tpg_id %in% c(global_pat_id, global_pat_id2), ]

# Bi/multivoltinism
# 3 TPGs AUS, NOA, NZ
# with plastron and spiracle, small
disting_traits[defining_tc %like% "bi/multivoltinism", ] # %>%
#   .[!(tpg_id %in% global_pat_id | tpg_id %in% global_pat_id2), ]
trait_CONT[continent == "NOA" & group == 2, ] %>% 
  .[affinity > 0.5, ]

# # Crawling & large, mostly predators
# disting_traits[defining_tc %like% "large", ]
# trait_CONT[continent == "SA" & group == 5 & trait == "feed_predator",]
disting_traits[defining_tc %like% "large", ] 

# ### Between Australia, South Africa, and New Zealand ----
# Specific TPGs similar in AUS, NZ & SA (more recent geological history)
# Not really something specific
# disting_traits[continent %in% c("AUS", "NZ", "SA"), ] %>%
#   .[order(continent), ]


### On 2 regions/continents ----

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


#### Traits that only occur in TPG in one continent/region & special TPGs ----
## SA
disting_traits[defining_tc %like% "gatherer", ] %>%
  .[!tpg_id %in% c(global_pat_id, global_pat_id2), ]

disting_traits_compl[, .N, by = c("trait", "continent")] %>% 
  .[order(trait, N), ]

# Streamlined only in NOA
# with filterer, gills, univoltinism
disting_traits[defining_tc %like% "streamlined", ] #%>%
#   .[!(tpg_id %in% global_pat_id | tpg_id %in% global_pat_id2), ]

# Semivoltinism only in EUR
disting_traits[defining_tc %like% "semivoltinism", ]


