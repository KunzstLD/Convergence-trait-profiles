# __________________________________________________________________________________________________
# Analysis of TPGs ----
# Identify traits that characterize the groups obtained through cluster analysis
# Clustering + Ordination  (principal coordinate analysis) together?
# Consider correlations among traits (see Wilkes Paper!)
# __________________________________________________________________________________________________

# Read in trait data with group assignment
trait_CONT <- readRDS(file.path(data_cache, "trait_dat_grp_assig.rds"))

# __________________________________________________________________________________________________
# TPG's ----
# __________________________________________________________________________________________________

# Unique families per continent and group
trait_CONT[, nr_families := uniqueN(family), by = .(continent, group)]

# Delineated TPGs per continent
unique(trait_CONT[order(group), .(group), by = continent])
trait_CONT[, uniqueN(family), by = continent]

# Scaled to 100
trait_CONT[, uniqueN(family), by = continent][, .(continent, 100 / V1 * c(8, 9, 11, 7, 10))]

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
    .SDcols = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
  ] %>%
  .[, .(continent, order, occ_in_cluster)] 
# cluster_occ[order %in% c("Ephemeroptera", "Trichoptera", "Diptera", "Plecoptera"), ]

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
  "11" = "TPG 11",
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

# Calculate nr of traits per tpg for upcoming comparisons
# Overall, TPGs have between 2 and 7 defining traits
disting_traits[, nr_traits_group := .N, 
               by = .(continent, group)]
unique(disting_traits[, nr_traits_group, by = c("continent", "group")])
disting_traits[, .(
  min_def_traits = min(nr_traits_group),
  max_def_traits = max(nr_traits_group)
), by = continent]

# saveRDS(disting_traits,
#          file = file.path(data_cache, "def_traits.rds"))

### Overview ----
# for our initial criterion
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
    "shredder",
    "gatherer",
    "burrower",
    "streamlined",
    "bi/multivoltinism",
    "filterer",
    "sessil"
  ))
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

# Exact duplicates 
disting_traits[duplicated(defining_tc) | duplicated(defining_tc, fromLast = TRUE), ]

# Fuzzy matching
library(tidystringdist)
# library(tidyverse)

tab_fmatch_tc <- disting_traits[, tidy_comb_all(defining_tc)] %>%
  tidy_stringdist() 
setDT(tab_fmatch_tc)
# Go through fuzzy matching table manually and select
# potential similar tc's
tab_fmatch_tc[order(lv), head(V2, n = 5), by = V1] %>% View()


# 17 out of 39 TPGs are characterised by cylindrical, aquatic ovipositon & crawling
unique(disting_traits[, .(continent, group, def_traits)]) %>% 
  .[grepl("(?=.*cylindrical)(?=.*aqu)(?=.*crawl)", def_traits, perl = TRUE), ]

# Remaining 
unique(disting_traits[, .(continent, group, def_traits)]) %>% 
  .[!grepl("(?=.*cylindrical)(?=.*aqu)(?=.*crawl)", def_traits, perl = TRUE), ]

# "unique" TPGs
unique(disting_traits[, .(continent, group, def_traits)]) %>% 
  .[grepl("(?=.*ter)|(?=.*gatherer)|(?=.*filter)|(?=.*streamlined)", def_traits, perl = TRUE), ]
unique(disting_traits[, .(continent, group, def_traits)]) %>% 
  .[grepl("(?=.*flattened)", def_traits, perl = TRUE), ]

# add id
disting_traits[, tpg_id := paste0(continent, "_", group)]

# Defining traits for each TPG long format
disting_traits_lf <-
  dcast(disting_traits[, .(continent, group, tpg_id, trait)],
    ... ~ trait,
    fun.aggregate = fun_binary_length
  )
# disting_traits_lf[, apply(.SD, 1, function(y) sum(y)), 
#                   .SDcols = patterns("bf|feed|locom|ovip|resp|size|volt")]

### Detailed look ----
# Returns the proportion to which TPGs share the same defining traits
test <- disting_traits_lf[, apply(.SD, 1, function(y)
  ifelse(y == 1, names(y), NA)), # select the defining traits
  .SDcols = patterns("bf|feed|locom|ovip|resp|size|volt")]
test <- as.data.frame(test)

# assign id, makes life easier for later 
names(test) <-  disting_traits_lf$tpg_id

identical_new <- function(x, y) {
  x <- na.omit(x)
  y <- na.omit(y)
  if (length(x) >= length(y)) {
    equiv <- unique(x) %in% unique(y)
  } else{
    equiv <- unique(y) %in% unique(x)
  }
  sum(equiv) / length(equiv)
}

# get those tpgs that share at least 50 % of defining traits
comp_tpgs <- combn(names(test), 2)

out <- list()
for(i in 1:ncol(comp_tpgs)){
  var1 <- comp_tpgs[, i][[1]]
  var2 <- comp_tpgs[, i][[2]]
  out[[paste0(var1, "_", var2)]] <- identical_new(test[, var1], test[, var2])
}
out <- unlist(out)
share_traits <- as.data.table(out, keep.rownames = TRUE)
setnames(share_traits,
         "out",
         "amount_shared_traits")
share_traits[, `:=`(
  tpg = sub("(.+\\_[0-9]{1,})(\\_)(.+\\_[0-9]{1,})", "\\1", rn),
  tpg_match = sub("(.+\\_[0-9]{1,})(\\_)(.+\\_[0-9]{1,})", "\\3", rn)
)]
share_traits[, continent_match := sub("(.+)(\\_)(.+)", "\\1", tpg_match)]

# tpgs that occur on three or four continents and share at least 37.5 % (i.e, at least 3 of 8)
# of their defining traits
tpgs_shared <- share_traits[amount_shared_traits >= 0.375, ] %>%
  .[, paste0(unique(continent_match), collapse = ","), by = tpg] %>%
  .[V1 == "AUS,EU,NOA,NZ,SA" | V1 == "EU,NOA,NZ,SA", tpg]
tpg_across_cont <-
  share_traits[tpg %in% tpgs_shared &
                 amount_shared_traits >= 0.375, tpg_match]

## Defining traits that occur most often ----
# this is the basis for the selection afterwards
disting_traits[tpg_id %in% tpg_across_cont, .N, by = trait] %>% 
  .[order(-N), ]

## TPG similarity across all continents ----
# Various combinations of bf_cylindrical, locom_crawl, ovip_aqu
# additionally: , feed_predator,resp_teg, resp_gil, size_medium, feed_herbivore
disting_traits[tpg_id %in% tpg_across_cont &
                 grepl("(?=.*bf\\_cylindrical)(?=.*ovip\\_aqu)(?=.*locom\\_crawl)",
                       def_traits,
                       perl = TRUE),]

## Associated traits that occurred on all continents with these traits ---- 
# - size small (mostly with locomotion crawling, EU1 is the exception;
# most of the TPGs have bf_cylindrical, exceptions are NZ9, NOA6 and NOA7)
# - volt_uni (some don't have locom_crawl, e.g. NZ1)
# - size_medium
# - resp_gil
# - feed_herbivore
disting_traits[tpg_id %in% tpg_across_cont &
                 grepl("(?=.*bf\\_cylindrical)(?=.*ovip\\_aqu)(?=.*size\\_small)",
                       def_traits,
                       perl = TRUE),]

disting_traits[tpg_id %in% tpg_across_cont &
                 grepl("(?=.*ovip_aqu)(?=.*locom_crawl)(?=.*size\\_small)",
                       def_traits,
                       perl = TRUE),] %>%
  .[grepl("bf_flattened", def_traits), ]

disting_traits[tpg_id %in% tpg_across_cont &
                 grepl("(?=.*bf\\_cylindrical)(?=.*ovip\\_aqu)(?=.*volt\\_uni)",
                       def_traits,
                       perl = TRUE),]

disting_traits[tpg_id %in% tpg_across_cont &
                 grepl("(?=.*bf\\_cylindrical)(?=.*ovip\\_aqu)(?=.*size_medium)",
                       def_traits,
                       perl = TRUE),]

disting_traits[tpg_id %in% tpg_across_cont &
                 grepl("(?=.*locom\\_crawl)(?=.*ovip\\_aqu)(?=.*resp_gil)",
                       def_traits,
                       perl = TRUE),]

disting_traits[tpg_id %in% tpg_across_cont &
                 grepl("(?=.*bf\\_cylindrical)(?=.*ovip\\_aqu)(?=.*feed_herbivore)",
                       def_traits,
                       perl = TRUE),]

# Feed predator is somehow on the edge (i.e. TPGs differ from each other relativly strongly)
# occurs in TPGsin EU, NOA & NZ
# in AUS only with ovip_ter & size large (different to many other TPGs that have predator)
disting_traits[tpg_id %in% tpg_across_cont &
                 grepl("(?=.*bf\\_cylindrical)(?=.*ovip\\_aqu)(?=.*feed_predator)",
                       def_traits,
                       perl = TRUE),]

disting_traits[tpg_id %in% tpg_across_cont &
                 grepl("(?=.*feed_predator)",
                       def_traits,
                       perl = TRUE),] %>% 
  .[grepl("(?=.*locom_crawl)(?=.*bf\\_cylindrical)",
          def_traits,
          perl = TRUE), ]

## Three region similarity ----

# resp_teg in TPGs in EU, NOA & NZ
disting_traits[tpg_id %in% tpg_across_cont &
                 grepl("(?=.*bf\\_cylindrical)(?=.*ovip\\_aqu)(?=.*resp_teg)",
                       def_traits,
                       perl = TRUE),]

# bf_flattened (EU, NOA, SA)
disting_traits[tpg_id %in% tpg_across_cont &
                 grepl("(?=.*bf_flattened)",
                       def_traits,
                       perl = TRUE), ]

# Two region similarity

# feed_shredder -> part of the combination of defining traits that occurs on all continents
disting_traits[tpg_id %in% tpg_across_cont &
                 grepl("(?=.*feed\\_shredder)",
                       def_traits,
                       perl = TRUE), ]

# size_large -> not really similar TPGs?
disting_traits[tpg_id %in% tpg_across_cont &
                 grepl("(?=.*size\\_large)",
                       def_traits,
                       perl = TRUE), ]

## Unique TPGs ----
# Except for ovip_ter, those others are part of the already mentioned TPGs
# ovip ter
disting_traits[tpg_id %in% tpg_across_cont &
                 grepl("(?=.*ovip\\_ter)",
                       def_traits,
                       perl = TRUE), ]

# feed filter
disting_traits[tpg_id %in% tpg_across_cont &
                 grepl("(?=.*feed\\_filter)",
                       def_traits,
                       perl = TRUE), ]

# locom_sessil
disting_traits[tpg_id %in% tpg_across_cont &
                 grepl("(?=.*locom\\_sessil)",
                       def_traits,
                       perl = TRUE), ]

# bf_spherical
disting_traits[tpg_id %in% tpg_across_cont &
                 grepl("(?=.*bf\\_spherical)",
                       def_traits,
                       perl = TRUE), ]

# feed_gatherer -> largely similar to EU_TPG8
disting_traits[tpg_id %in% tpg_across_cont &
                 grepl("(?=.*feed\\_gatherer)",
                       def_traits,
                       perl = TRUE) | (tpg_id == "EU_8"),]

# look at those TPGs that don't really share traits with other TPGs (the ones
# we looked at so far share at least 37.5 percent)
disting_traits[, .N, by = trait] %>% 
  .[order(-N), ]

# volt bi multi -> part of the aformentioned TPGs
disting_traits[grepl("(?=.*volt\\_bi\\_multi)",
                       def_traits,
                       perl = TRUE), ]

# bf_streamlined -> part of the aformentioned TPGs
disting_traits[grepl("(?=.*bf\\_streamlined)",
                       def_traits,
                       perl = TRUE), ]

## Non-defining traits for TPGs ----
unique(trait_CONT$trait[!trait_CONT$trait %in% disting_traits$trait])


### Are there similarities across continents? ----
# output_tpg <- list()
# for(cont in c("AUS", "EU", "NOA", "NZ")) { # "AUS",
#   groups <- unique(disting_traits[continent == cont, group])
#   similar_tpg <- list()
#   for (i in groups) {
#     trait_comb <-
#       unique(disting_traits[continent == cont & group == i, trait])
#     
#     # use binary coded disting_traits dataset
#     # 1 indicates that the trait is a defining trait based on the previously established criterion
#     tpg <- melt(disting_traits_lf, id.vars = c("continent", "group")) %>%
#       .[value == 1,] %>%
#       .[, n_defTraits := .N, by = .(continent, group)] %>% # Number of defining traits of each TPG
#       .[variable %in% trait_comb,] %>%
#       .[order(continent, group), ] %>%
#       .[, n_subsTraits := .N, by = c("continent", "group")] %>% 
#       .[n_subsTraits == length(trait_comb), ] %>% # Number of the defining traits that should occur 
#       # in other TPGs (i.e. the other TPGs can have additional defining traits)
#       #.[n_defTraits == length(trait_comb), ] %>%  # Number of defining traits of each TPG should be equal to 
#         # the number of the def. traits that are currently searched for, i.e. we obtain the whole set 
#         # of defining traits for a given TPG (i.e. identical TPG across regions)
#       .[uniqueN(continent) == 4, ]  
#     if(nrow(tpg) > 0){
#       similar_tpg[[i]] <- tpg 
#     } 
#   }
#   output_tpg[[cont]] <- similar_tpg
# }
# 
# #### Postprocessing ----
# # Get those TPGs that are similar across all continents
# output_tpg4 <- output_tpg
# simtpg_across_cont <- lapply(output_tpg4, rbindlist) %>% 
#   .[["AUS"]] %>% 
#   .[, .(continent, group)] %>% 
#   unique(.) %>% 
#   .[, paste0(continent, "_", group)]

# filter only for those groups that are similar between two or three regions 
# (not just subsets of those TPGs that are similar across all continents)
# Should be noted that the TPGs are often not entirely similar, so "new" similarities can
# be found (i.e. with additional traits) between two or three regions
# output_tpg3 <- output_tpg
# lapply(output_tpg3, function(y)
#   lapply(y, function(z)
#     z[, id := paste0(continent, "_", group)]))
# 
# lapply(output_tpg3, function(y)
#   lapply(y, function(z)
#     z[!id %in% simtpg_across_cont,]))
# 
# output_tpg2 <- output_tpg
# lapply(output_tpg2, function(y)
#   lapply(y, function(z)
#     z[, id := paste0(continent, "_", group)]))
# 
# lapply(output_tpg2, function(y)
#   lapply(y, function(z)
#     z[!id %in% simtpg_across_cont,]))

# ___________________________________________________________________________________________________
#  TPGs that occur in all continents ----
# Plot in single heatmaps?
# ___________________________________________________________________________________________________

#*****************************************************************
## Crawling, aquatic oviposition & cylindrical body form ----
# Basically this trait combination with slight
# variations (size_medium & size small sometimes)
#*****************************************************************
tpg1_traits <- c(
  "bf_cylindrical",
  "ovip_aqu",
  "locom_crawl", 
  "size_small",
  "size_medium",
  "volt_uni"
)
tpg_acr_cont <-
  trait_CONT[trait %in% tpg1_traits & (
    (continent == "AUS" & group %in% c(1, 2, 3)) |
      (continent == "EU" & group %in% c(3, 4, 5, 6, 7)) |
      (continent == "NOA" & group %in% c(3, 4, 5, 8)) |
      (continent == "NZ" & group %in% c(3, 5, 6, 10))
  ), ]
# tpg_acr_cont$continent %>% unique

# change order of continent column for plotting
tpg_acr_cont[, continent := paste0(continent, "_", group)]
tpg_acr_cont$continent %>% unique()
tpg_acr_cont[, continent := factor(continent,
  levels = c(
    "AUS_1",
    "AUS_2",
    "AUS_3",
    "EU_3",
    "EU_4",
    "EU_5",
    "EU_6",
    "EU_7",
    "NOA_3",
    "NOA_4", 
    "NOA_5",
    "NOA_8",
    "NZ_3",
    "NZ_5",
    "NZ_6",
    "NZ_10"
  )
)]

# Heatmap
tpg_names1 <- c(
  "AUS_1" = "AUS_TPG1",
  "AUS_2" = "AUS_TPG2",
  "AUS_3" = "AUS_TPG3",
  "EU_3" = "EU_TPG3",
  "EU_4" = "EU_TPG4",
  "EU_5" = "EU5_TPG",
  "EU_6" = "EU6_TPG",
  "EU_7" = "EU7_TPG",
  "NOA_3" = "NA_TPG3",
  "NOA_4" = "NA_TPG4",
  "NOA_5" = "NA_TPG5",
  "NOA_8" = "NA_TPG8",
  "NZ_3" = "NZ_TPG3",
  "NZ_5" = "NZ_TPG5",
  "NZ_6" = "NZ_TPG6",
  "NZ_10" = "NZ_TPG10"
)
plot_list <- list()
for(cont in c("AUS", "EU", "NOA", "NZ")) {
  plot_list[[cont]] <- fun_heatmap_tpg(data = tpg_acr_cont[continent %like% cont,],
                                       facet_names = tpg_names1)
}
(plot_list[["AUS"]] + plot_list[["EU"]]) / (plot_list[["NOA"]] + plot_list[["NZ"]])
ggsave(
  filename = file.path(
    data_paper,
    "Graphs",
    "Heatmap_TPG_across_all_continents_50.png"
  ),
  width = 60,
  height = 30,
  units = "cm"
)

#*****************************************************************
## Predators, small size, respiration with plastron and spiracle ----
#*****************************************************************
tpg2_traits <- c(
  "feed_predator",
  "size_small",
  "resp_pls_spi",
  "ovip_aqu",
  "volt_uni",
  "locom_swim"
)

# Similar TPGs:
# AUS_TPG5, EU_TPG1, EU_TPG2,  NA_TPG2,  NZ_TPG5
tpg2_acr_cont <-
  trait_CONT[trait %in% tpg2_traits & (
    (continent == "AUS" & group == 5) |
      (continent == "EU" & group %in% c(1, 2)) |
      (continent == "NOA" & group == 2) |
      (continent == "NZ" & group %in% 6)
  ), ]

# change order of continent column for plotting
tpg2_acr_cont[, continent := paste0(continent, "_", group)]
tpg2_acr_cont$continent %>% unique()
tpg2_acr_cont[, continent := factor(continent,
                                    levels = c("AUS_5",
                                               "EU_1",
                                               "EU_2",
                                               "NOA_2",
                                               "NZ_6"))]
tpg2_acr_cont[, trait := factor(
  trait,
  levels = c(
    "volt_uni",
    "ovip_aqu",
    "locom_swim",
    "feed_predator",
    "size_small",
    "resp_pls_spi"
  )
)]

# Heatmap
tpg2_names <- c(
  "AUS_5" = "AUS_TPG5",
  "EU_1" = "EU_TPG1",
  "EU_2" = "EU_TPG2",
  "NOA_2" = "NA_TPG2",
  "NZ_6" = "NZ_TPG6")
plot_list <- list()
for(cont in c("AUS", "EU", "NOA", "NZ")) {
  plot_list[[cont]] <- fun_heatmap_tpg(data = tpg2_acr_cont[continent %like% cont,],
                                       facet_names = tpg2_names)
}
(plot_list[["AUS"]] + plot_list[["EU"]]) / (plot_list[["NOA"]] + plot_list[["NZ"]])
ggsave(
  filename = file.path(
    data_paper,
    "Graphs",
    "Heatmap_TPG_across_all_continents_40.png"
  ),
  width = 55,
  height = 35,
  units = "cm"
)

# ___________________________________________________________________________________________________
# TPGs that occur in specific regions ----
# ___________________________________________________________________________________________________

#*****************************************************************
## Groups with locomotion swimming ----
# EUR TPG 1 (also ovip_aqu, bf_cylindrical)
# AUS TPG 5 
# NOA TPG 3 (also ovip_aqu ,bf_clyindrical)
#*****************************************************************
# locom_swim, feed_predator, resp_pls_spi, size small
# trait_CONT[continent == "AUS" & group == 5 & prop_taxa_high_aff >= 0.4, ]
# trait_CONT[continent == "EU" & group == 1 & prop_taxa_high_aff >= 0.4, ]
# unique(trait_CONT[(continent == "EU" &
#                      group == 1 & prop_taxa_high_aff >= 0.4) |
#                     (continent == "AUS" &
#                        group == 5 & prop_taxa_high_aff >= 0.5) |
#                     (continent == "NOA" &
#                        group == 3 &
#                        prop_taxa_high_aff >= 0.4),]) %>%
#   .[trait %in% c("locom_swim", "feed_predator", "resp_pls_spi", "size_small"),]
# 
# trait_CONT[continent == "AUS" &
#              group == 5 & prop_taxa_high_aff >= 0.4, ]
# 
# # .(continent, group, trait)
# 
# #*****************************************************************
# ## Groups with feed shredder ----
# # AUS: 3
# # EU: 3, partly in 4 and 6
# # NoA: 7
# #*****************************************************************
# 
# 
# #*****************************************************************
# ## Groups with volt semi ----
# # AUS: 3
# # EU: 3, partly in 4 and 6
# # NoA: 7
# #*****************************************************************
# 
# 
# 
# #*****************************************************************
# ## Groups with terrestrial oviposition ----
# # Actually only AUS TPG5
# # EU_TPG4 has some taxa, also NOA_TPG2, and NZ_TPG4/NZ_TPG5
# #*****************************************************************
# 
# groups_ter_ovip <- rbind(
#   unique(trait_CONT[(continent == "EU" &
#     group == 4 &
#     prop_taxa_high_aff >= 0.5), .(trait, continent, group, prop_taxa_high_aff)]),
#   unique(trait_CONT[(continent == "AUS" &
#     group %in% c(5, 6) &
#     prop_taxa_high_aff >= 0.35), .(trait, continent, group, prop_taxa_high_aff)]),
#   unique(trait_CONT[(continent == "NOA" &
#     group == 2 &
#     prop_taxa_high_aff >= 0.25), .(trait, continent, group, prop_taxa_high_aff)]),
#   unique(trait_CONT[continent == "NZ" &
#     group == 5, .(trait, continent, group, prop_taxa_high_aff)])
# )
# groups_ter_ovip[continent == "AUS", ] %>% 
#   .[order(group), ]
# 
# #*****************************************************************
# ## Groups with bf flattened ----
# # AUS: 4 (partly)
# # EU: 8, 4 (partly)
# # NoA: 1, 2, and 6
# # NZ: 2, 7, 8 (all partly)
# # EU & NoA: bf_flattened, feed_predator, locom_crawl, ovip_aqu volt_uni 
# # 2 groups NoA: size small & size large, feed_gatherer & feed_predator 
# #*****************************************************************
# unique(trait_CONT[(continent == "NOA" & group %in% c(2, 6) &
#               prop_taxa_high_aff >= 0.5) |
#              (continent == "EU" & group == 7 &
#                 prop_taxa_high_aff >= 0.5) |
#              (continent == "AUS" & group == 4 &
#                 prop_taxa_high_aff >= 0.4),
#            .(continent, group, trait)]
# )









