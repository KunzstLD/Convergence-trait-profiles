# ___________________________________________________________________________
# Post analysis ----
# mainly for the discussion part of the paper
# ___________________________________________________________________________

## Trait space analysis ----

### AUS trait syndromes that deviates from other trait spaces  ---- 
Trait_AUS_agg <- readRDS(file.path(data_in, "Trait_AUS_agg.rds"))
AUS_outside_taxa_lf <- melt(Trait_AUS_agg[family %in% c(
  "Veliidae",
  "Pleidae",
  "Hydrometridae",
  "Notonectidae",
  "Dytiscidae",
  "Gerridae",
  "Hygrobiidae",
  "Eustheniidae",
  "Synthemistidae",
  "Libellulidae",
  "Corduliidae"
), ], id.vars = c("family", "order"), variable.name = "traits")
AUS_outside_taxa_lf[, grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", traits)]

ggplot(AUS_outside_taxa_lf, aes(
  x = factor(traits),
  y = value
)) +
  geom_jitter(width = 0.05) +
  geom_violin(alpha = 0.2) +
  coord_flip() +
  facet_wrap(grouping_feature ~ .,
    scales = "free"
  ) + # labeller = as_labeller(facet_names)
  labs(
    x = "Traits",
    y = "Affinity"
  ) +
  ggtitle(paste("Trait affinity distribution")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    axis.text.y = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    legend.title = element_text(
      family = "Roboto Mono",
      size = 16
    ),
    legend.text = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    strip.text = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    plot.title = element_text(
      family = "Roboto Mono",
      size = 16
    ),
  )

# To which TPG do these taxa belong? 
# AUS TPG 5 & 6
trait_CONT <- readRDS(file.path(data_cache, "trait_dat_grp_assig.rds"))
unique(trait_CONT[continent == "AUS" & family %in% aus_scores_syndrom$family, group])


### Taxa in high density area ----
# Characterize trait spaces in these areas
contour_25 <- readRDS(file.path(data_cache, "contour_25_TS.rds"))
contour_50 <- readRDS(file.path(data_cache, "contour_50_TS.rds"))
pcoa_scores[, family := sub("(.+)(\\_)(.+)", "\\3", id)]

# N per continent (include total N)
pcoa_scores[, total_N := .N, by = continent]
unique(pcoa_scores[contour_25, .N/total_N, by = continent])
unique(pcoa_scores[contour_50, .N/total_N, by = continent])

# continue with 50 % contour (i.e. smallest region with 50 % probability mass)
trait_CONT <- readRDS(file.path(data_cache, "trait_dat_grp_assig.rds"))
trait_CONT[, id := paste0(continent, "_", family)]

# Search and summarize trait profiles for the taxa in the 50 % probability occurrence area
trait_CONT_50 <- trait_CONT[id %in% pcoa_scores[contour_50, id],]

# Which groups are present per continent? 
group_50 <- unique(trait_CONT_50[, group, by = continent]) %>% 
  .[order(continent, group), ]

unique(trait_CONT_50[, .(group, family), by = continent]) %>%
  .[, paste(family, collapse = ", "), by = .(continent, group)] %>%
  .[order(continent, group), ] %>% 
  fwrite(.,
         file = file.path(data_out, "TPGs_in_50_contour.csv"),
         sep = ";")

# Then check defining traits for each group 
trait_CONT_50[, nr_families := uniqueN(family), by = .(continent, group)]
trait_CONT_50[affinity > 0.5,
              prop_taxa_high_aff := .N / nr_families,
              by = .(continent, group, trait)]

# Traits by which these families are mostly defined 
trait_CONT_50[prop_taxa_high_aff >= 0.5, .N, by = trait] %>% 
  .[order(-N), ] 
trait_CONT_50[trait %in% c("ovip_aqu",
                           "volt_uni",
                           "bf_cylindrical",
                           "locom_crawl",
                           "resp_gil",
                           "size_small"), .N, by = .(continent, group)] %>% 
  .[order(continent, group)]

# ___________________________________________________________________________
## Global pattern in TPGs ----
# What are the dominant grouping features for every dataset?
# ___________________________________________________________________________
trait_CONT <- readRDS(file.path(data_cache, "trait_dat_grp_assig.rds"))

# life-history, mobility, morphology
# voltinism, locomotion, body form
trait_CONT[grouping_feature %in% c("bf", "size", "volt", "locom"), ]

disting_traits <-
  readRDS(file = file.path(data_cache, "def_traits.rds"))

# check how certain grouping features are distributed across continents
disting_traits[grouping_feature %in% c("bf", "locom", "volt", "size"), ] %>%
  .[, .(trait, prop_taxa_high_aff, .N), by = c("continent", "group")] %>%
  .[N >= 4, ]

# AUS: gf according to most important traits: bf, locom, ovip, size
disting_traits[grouping_feature %in% c("bf", "locom", "ovip", "size"), ] %>% 
  .[, .(trait, prop_taxa_high_aff, .N), by = c("continent", "group")] %>%
  .[N >= 4, ]

# NZ: feed, ovip, resp
disting_traits[grouping_feature %in% c("feed", "ovip", "resp"),] %>%
  .[, .(trait, prop_taxa_high_aff, .N), by = c("continent", "group")] %>%
  .[N >= 3,]


# ___________________________________________________________________________
## Patterns in the aggregated trait datasets ----
# ___________________________________________________________________________
trait_data_bind <- readRDS(file.path(data_cache,
                                     "trait_data_ww_bind.rds"))

# Calculate nr of taxa & prop of orders
trait_data_bind[, nr_taxa := .N, by = "continent"]
trait_data_bind[, nr_fam_per_order := .N, by = .(continent, order)]
trait_data_bind[, prop_order := nr_fam_per_order / nr_taxa]

## Terrestrial oviposition ----
trait_data_bind[ovip_ter >= 0.5, .N, by = order]
trait_data_bind[ovip_ter >= 0.5, .N, by = continent]
trait_data_bind[ovip_ter >= 0.5, .N, by = c("continent", "order")]
trait_data_bind[ovip_ter >= 0.5 & continent == "EU", ]

# Australia 
# also a lot of Odonates and Hemipterans
trait_data_bind[ovip_ter >= 0.5 & continent == "AUS", .N, by = order]

# Although there are mainly dipterans & coleopterans in the NZ dataset
# only few dipterans exhibit mainly terrestrial oviposition (i.e. > 0.5) 
trait_data_bind[ovip_ter > 0 & continent == "NZ", ] %>% 
  .[order %in% c("Coleoptera", "Diptera"), .N, by = order]
trait_data_bind[ovip_ter >= 0.5 & continent == "NZ", ]

# ___________________________________________________________________________
## Variability in traits ----
# Most important:
# AUS: bf_cylindrical|locom_crawl|ovip_ter|size_medium|size_small
# EU: feed predator|locom_crawl|locom_swim|size_small|volt_semi
# NOA: bf_cylindrical|bf_flattened|size_large|size_medium|size_small
# NZ: feed_herbivore|feed_predator|ovip_aqu|ovip_ter|resp_teg
# ___________________________________________________________________________

# five most important
trait_data_bind[, lapply(.SD, sd), 
                .SDcols = patterns("bf|feed|resp|ovip|volt|size|locom"), 
                by = continent] %>% 
  melt(., id.var = "continent") %>% 
  .[order(continent, -value), head(variable, n = 10), by = continent]

# Variation in body size
trait_data_bind[, lapply(.SD, mean), .SDcols = patterns("size"),
                by = "continent"]

trait_data_bind[, lapply(.SD, sd), .SDcols = patterns("size"),
                by = "continent"]

ggplot(trait_data_bind) +
  geom_point(aes(x = family, y = size_large)) +
  facet_wrap( ~ continent) +
  theme_bw()


# ___________________________________________________________________________
## Comparison to non-aggregated data ----
# ___________________________________________________________________________
# go back to the non-aggregate trait datasets to investigate certain patterns further

trait_non_agg <- load_data(pattern = ".*\\.rds",
                           path = file.path(data_in, "Not_aggregated"))
trait_non_agg <- rbindlist(trait_non_agg, idcol = "dataset", fill = TRUE)

# Remove orders Neuroptera and Megaloptera
trait_non_agg <- trait_non_agg[!order %in% c("Neuroptera", "Megaloptera"), ]

# Load the aggregated datasets as well
trait_data_bind <- readRDS(file.path(data_cache,
                                     "trait_data_ww_bind.rds"))

# filter for aquatic insects
aq_insects <- c(
  "Ephemeroptera",
  "Hemiptera",
  "Odonata",
  "Trichoptera",
  "Coleoptera",
  "Plecoptera",
  "Diptera",
  "Megaloptera",
  "Neuroptera"
)

### Ovoviviparity ----
# 223 taxa that exhibit ovoviviparity to a certain degree
sum(trait_non_agg[ovip_ovo > 0, .N, by = dataset]$N)

# Only 36 aquatic insects
trait_non_agg[ovip_ovo > 0 & order %in% aq_insects, .N, by = dataset]
ovo_families <-
  unique(trait_non_agg[ovip_ovo > 0 & order %in% aq_insects, family])

trait_non_agg[order %in% aq_insects &
                !is.na(family), .SD,
              .SDcols = patterns("ovip|species|family|order|dataset")] %>%
  .[family %in% ovo_families, .(
    median(ovip_aqu, na.rm = TRUE),
    median(ovip_ovo, na.rm = TRUE),
    median(ovip_ter, na.rm = TRUE)
  ),
  by = c("dataset", "family")]


### Terrestrial oviposition ----

# Overall: families with taxa with ovip_ter >= 0.5
# families with taxa ovip_ter >= 0.5: AUS 43, NZ 11, EU, 17, NOA 12
trait_non_agg[ovip_ter >= 0.5, uniqueN(family), by = dataset]

# Families of Diptera & Coleoptera order with taxa with ovip_ter >= 0.5
trait_non_agg[ovip_ter >= 0.5 &
                order %in% c("Diptera", "Coleoptera"),
              uniqueN(family), by = "dataset"]

# Overall number of dipterans and coleoptera families in each dataset
trait_non_agg[order %in% aq_insects, uniqueN(family),
              by = c("dataset", "order")] %>%
  .[order(dataset, order), ] %>%
  data.table::dcast(., ... ~ dataset) %>%
  fwrite(.,
         file = file.path(data_out, "summary_non_aggregated.csv"),
         sep = ";")

## Comparison to aggregated datasets ----

### Nr. of families ----

# Overall
# Some orders not represent in aggregated trait datasets
trait_data_bind[, uniqueN(family), by = continent]
trait_non_agg[order %in% aq_insects, uniqueN(family), by = dataset]

# Per order
nfamilies_agg <- trait_data_bind[, uniqueN(family), by = c("continent", "order")]
nfamilies_non_agg <-
  trait_non_agg[order %in% aq_insects, uniqueN(family), by = c("dataset", "order")]
nfamilies_non_agg[, continent := fcase(
  dataset == "Trait_AUS_harmonized.rds", "AUS",
  dataset == "Trait_freshecol_2020_pp_harmonized.rds", "EU",
  dataset == "Traits_US_LauraT_pp_harmonized.rds", "NOA",
  dataset == "Trait_NZ_pp_harmonized.rds", "NZ",
  dataset == "Trait_SA_pp_harmonised.rds", "SA"
)]
nfamilies_agg[nfamilies_non_agg,
              `:=`(occur_ratio = V1 / i.V1,
                   nr_families_non_agg = i.V1),
              on = c("continent", "order")]
setnames(nfamilies_agg,
         c("continent", "order", "V1"),
         c("Continent", "Order", "nr_families_agg"))
setcolorder(
  nfamilies_agg,
  c(
    "Continent",
    "Order",
    "nr_families_agg",
    "nr_families_non_agg",
    "occur_ratio"
  )
)
nfamilies_agg[, occur_ratio := round(occur_ratio, digits = 2)]
fwrite(
  nfamilies_agg[order(Continent, Order),],
  file.path(data_paper, "Tables", " n_families_agg_non_agg_June_9.csv"),
  sep = ";"
)

nfamilies_agg[, .(
  sum_non_agg = sum(nr_families_non_agg),
  sum_agg = sum(nr_families_agg)
),
by = Continent] %>% 
  .[, .(Continent, sum_agg, sum_non_agg, sum_agg/sum_non_agg)]

# ___________________________________________________________________________
## Clustering with complete NZ dataset ----
# take most complete trait dataset and apply cluster analysis
# Do we find different strategies?
# ___________________________________________________________________________
trait_non_agg <- load_data(pattern = ".*\\.rds",
                           path = file.path(data_in, "Not_aggregated"))
Trait_NZ <- trait_non_agg$Trait_NZ_pp_harmonized.rds
Trait_NZ <- Trait_NZ[order %in% c(
  "Ephemeroptera",
  "Hemiptera",
  "Odonata",
  "Trichoptera",
  "Coleoptera",
  "Plecoptera",
  "Diptera",
  "Megaloptera",
  "Neuroptera"
), ]

# test how complete trait sets are 
completeness_trait_data(
  x = Trait_NZ,
  non_trait_cols = c("order",
                     "family",
                     "genus",
                     "species",
                     "unique_id")
)
Trait_NZ[, taxa := coalesce(species, genus, family, order)]

# change column order
setcolorder(
  Trait_NZ,
  c(
    "feed_shredder",
    "feed_gatherer",
    "feed_filter",
    "feed_predator",
    "feed_herbivore",
    "resp_teg",
    "resp_gil",
    "resp_pls_spi",
    "volt_semi",
    "volt_uni",
    "volt_bi_multi",
    "locom_swim",
    "locom_crawl",
    "locom_burrow",
    "locom_sessil",
    "ovip_ter",
    "ovip_ovo",
    "ovip_aqu",
    "size_medium",
    "size_small",
    "size_large",
    "bf_streamlined",
    "bf_flattened",
    "bf_cylindrical",
    "bf_spherical"
  )
)
Trait_NZ_cp <- copy(Trait_NZ)
Trait_NZ_cp[, c(
  "species",
  "genus",
  "family",
  "order",
  "unique_id",
  "dev_hemimetabol",
  "dev_holometabol"
) := NULL]

### Cluster analysis NZ ----
setDF(Trait_NZ_cp)

# add row.names
row.names(Trait_NZ_cp) <- Trait_NZ_cp$taxa
Trait_NZ_cp$taxa <- NULL

# Convert to ktab object
vec <- sub("\\_.*", "\\1", names(Trait_NZ_cp))
blocks <- rle(vec)$lengths
Trait_NZ_cp <- prep.fuzzy(Trait_NZ_cp, blocks)
Trait_NZ_cp <- ktab.list.df(list(Trait_NZ_cp))

# distance matrix and clustering
dist_mat_NZ <- dist.ktab(Trait_NZ_cp, type = "F") 
hc_NZ <- hclust(dist_mat_NZ, method = "ward.D")

# Get labels of dendrogram
dend_label_NZ <- hc_NZ %>%
  as.dendrogram() %>%
  labels()

# Optimal number of groups using the gap statistic
# Gap statistic always increases a bit
# Decide number of groups rather graphically
# gap <- clusGap(
#   x = as.matrix(dist_mat_NZ),
#   FUN = mycluster_hc,
#   K.max = 30,
#   B = 500
# )
# optimal_nog_NZ <- maxSE(gap$Tab[, "gap"],
#                      gap$Tab[, "SE.sim"],
#                      method = "Tibs2001SEmax")
# plot(gap)

# What kind of strategies could be found in the NZ trait dataset?
# dendrogram plot
dendro_NZ <- fun_dendrog_pl(
  hc = hc_NZ,
  optimal_nog = 1,
  labels = dend_label_NZ
)
plot(dendro_NZ, horiz = TRUE)

# trait profile groups for all continents
# cut at 5 revealed 14 groups
tpg_NZ <- data.table(
  taxa = names(cutree(
    dendro_NZ,
    h = 5,
    order_clusters_as_data = FALSE
  )),
  group = cutree(dendro_NZ,
                 h = 5,
                 order_clusters_as_data = FALSE)
)

# merge back tpgs to NZ trait data
Trait_NZ[tpg_NZ, group := i.group,
          on = "taxa"]
Trait_NZ[, c("dev_hemimetabol",
             "dev_holometabol") := NULL]

# Check which TPGs and defining trait combinations exist
Trait_NZ[, nr_families := .N, by = group]
Trait_NZ_lf <- melt(
  Trait_NZ,
  id.vars = c(
    "taxa",
    "order",
    "family",
    "genus",
    "species",
    "unique_id",
    "nr_families",
    "group"
  ),
  variable.name = "trait",
  value.name = "affinity"
)
Trait_NZ_lf[affinity > 0.5,
           prop_taxa_high_aff := .N / nr_families,
           by = .(group, trait)]
def_traits_NZ <- unique(Trait_NZ_lf[prop_taxa_high_aff >= 0.5,] %>%
                          .[order(-prop_taxa_high_aff), .(prop_taxa_high_aff,
                                                          nr_families),
                            by = .(group, trait)]) %>%
  .[order(group),]

# re-arrange for better representation
def_traits_NZ[, grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", trait)]
saveRDS(def_traits_NZ,
        file = file.path(data_cache, "NZ_complete_clustering.rds"))
def_traits_NZ %>% arrange(factor(
  grouping_feature,
  levels = c(
    "locom",
    "bf",
    "ovip",
    "size",
    "feed",
    "resp",
    "volt"
  ))) %>%
  .[, paste(trait, collapse = ", "),
    by = group] %>% 
  .[order(group), ]