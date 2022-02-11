# ___________________________________________________________________________
# Post analysis 
# mainly for the discussion part of the paper
# ___________________________________________________________________________

# AUS trait syndromes that deviates from other trait spaces  ----
# Identify families in the right upper corner of the Australia FS 
# which are not overlapping with the other taxa
# (could be probably solved easier)

hull <- readRDS(file.path(data_cache, "hull_pcoa.rds"))
setDT(hull)
# right and above the highest NZ coord (upper right) & the lowest righthand NZ coord

# highest NZ coord
coord1 <- hull[continent == "NZ",] %>%
  .[order(A2, decreasing = TRUE),] %>%
  .[1, .(A1, A2)]

# lowest righthand NZ coord (which is also the coord with the highest score on axis1)
coord2 <- hull[continent == "NZ",] %>%
  .[order(A1, decreasing = TRUE),] %>%
  .[1, .(A1, A2)]

# calculate slope
slope <- (coord2$A2 - coord1$A2)/(coord2$A1 - coord1$A1)
# y2-y1 = m * (x2-x1)
# if x1 = 0 then we can find the intercept 
# -> y2-y1 = m*(x2-0)
# -> y2 = m*(x2) + y1

# calculate intercept
# y = m * x + b
intercept <- coord2$A2 - slope*coord2$A1

# get values for Australia (right from highest NZ coord)
aus_scores_syndrom <-
  pcoa_scores[A1 > coord1$A1 & continent == "AUS", ]

# check for each A1 value (i.e. X-Coordinate) of AUS subset if it's not enclosed by the NZ space
# x = (y - b)/m
aus_scores_syndrom[, x_right_side_fsNZ := (A2 - intercept)/slope]
aus_scores_syndrom <- aus_scores_syndrom[A1 >= x_right_side_fsNZ, ]

# AUS taxa "above" the other trait spaces in the PCoA plot
# last two already captured! 
aus_scores_syndrom <- rbind(aus_scores_syndrom,
                            pcoa_scores[A2 > coord1$A2 & continent == "AUS", ][1,],
                            fill = TRUE)
aus_scores_syndrom[, family := sub("(.+)(\\_)(.+)", "\\3", id)]

# Check traits for these 
Trait_AUS_agg <- readRDS(file.path(data_in, "Trait_AUS_agg.rds"))

Trait_AUS_agg[family %in% aus_scores_syndrom$family, ] %>% 
  melt(., id.vars = c("family", "order"), variable.name = "traits") %>% 
  .[, grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", traits)] %>% 
  ggplot(., aes(x = factor(traits),
                y = value)) +
  geom_jitter(width = 0.05) +
  geom_violin(alpha = 0.2) +
  facet_wrap(grouping_feature ~ .,
             scales = "free"
  ) + # labeller = as_labeller(facet_names)
  labs(x = "Traits",
       y = "Affinity") +
  ggtitle(paste("Trait affinity distribution")) +
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
                              size = 14),
    plot.title = element_text(family = "Roboto Mono",
                              size = 16),
  )

# To which TPG do these taxa belong? 
# AUS TPG 5 & 6
trait_CONT <- readRDS(file.path(data_cache, "trait_dat_grp_assig.rds"))
unique(trait_CONT[continent == "AUS" & family %in% aus_scores_syndrom$family, group])

# ___________________________________________________________________________
# Global pattern in TPGs ----
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
# Patterns in the aggregated trait datasets ----
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
trait_data_bind[ovip_ter > 0 & continent == "NZ", ]
trait_data_bind[ovip_ter >= 0.5 & continent == "NZ", ]

# ___________________________________________________________________________
# Variability in traits ----
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

most_imp_vars[continents == "EU" & !is.na(five_most_imp), ]

# ___________________________________________________________________________
# Comparison to non-aggregated data ----
# ___________________________________________________________________________
# go back to the non-aggregate trait datasets to investigate certain patterns further

## Ovoviviparity ----
trait_non_agg <- load_data(pattern = ".*\\.rds",
                           path = file.path(data_in, "Not_aggregated"))

trait_non_agg <- rbindlist(trait_non_agg, idcol = "dataset", fill = TRUE)

# 223 taxa that exhibit ovoviviparity to a certain degree
sum(trait_non_agg[ovip_ovo > 0, .N, by = dataset]$N)

# filter for aquatic insects
aq_insects <- c("Ephemeroptera",
                "Hemiptera",
                "Odonata",
                "Trichoptera",
                "Coleoptera",
                "Plecoptera",
                "Diptera",
                "Megaloptera",
                "Neuroptera")

# Only 36 aquatic insects
trait_non_agg[ovip_ovo > 0 & order %in% aq_insects, .N, by = dataset]
ovo_families <- unique(trait_non_agg[ovip_ovo > 0 & order %in% aq_insects, family])

trait_non_agg[order %in% aq_insects &
                !is.na(family), .SD,
              .SDcols = patterns("ovip|species|family|order|dataset")] %>% 
  .[family %in% ovo_families, .(median(ovip_aqu, na.rm = TRUE), 
                                median(ovip_ovo, na.rm = TRUE),
                                median(ovip_ter, na.rm = TRUE)),
    by = c("dataset", "family")]


## Terrestrial oviposition ----
trait_non_agg[ovip_ter >= 0.5, .N, by = c("dataset", "order", "family")] %>% 
  .[order(dataset, order), sum(.N), by = c("dataset")]
# families with taxa ovip_ter >= 0.5: AUS 43, NZ 11, EU, 20, NOA 12

# overall number of dipterans and coleoptera families in each dataset
trait_non_agg[order %in% c("Diptera", "Coleoptera"), 
              .N, by = c("dataset", "order", "family")] %>% 
  .[, sum(.N), by = c("dataset", "order")]

