# __________________________________________________________________________________________________
# Comparing the complete trait profiles with the number of Köppen Geiger zones ----
# Calc. a distance between trait profiles of families and their distance to the 
# mean trait profile of their tpgs 
# TODO weighted by family number
# __________________________________________________________________________________________________

# Non aggregated trait data ----
trait_non_agg <- load_data(pattern = ".*\\.rds",
                           path = file.path(data_in, "Not_aggregated"))

# TPGs and their mean trait profiles
trait_CONT <- readRDS(file.path(data_cache, "trait_dat_grp_assig.rds"))

## EU ----
Trait_EU <- trait_non_agg$Trait_freshecol_2020_pp_harmonized.rds
Trait_EU <- Trait_EU[, .SD,
                     .SDcols = names(Trait_EU) %like% "locom|feed|resp|volt|size|bf|species|genus|family|order"] %>%
  .[order %in% c(
    "Ephemeroptera",
    "Hemiptera",
    "Odonata",
    "Trichoptera",
    "Coleoptera",
    "Plecoptera",
    "Diptera"
  ),]
Trait_EU[, taxa := coalesce(species, genus, family, order)]
Trait_EU[, feed_parasite := NULL]

# Search those families that could be used for TPG delineation
Trait_EU <- Trait_EU[family %in% unique(trait_CONT[continent == "EU", family]), ]
Trait_EU_lf <- melt(
  Trait_EU,
  id.vars = c("species",
              "genus",
              "family",
              "order",
              "taxa"),
  value.name = "affinity",
  variable.name = "trait"
)

## NOA ----
Trait_NOA <- trait_non_agg$Traits_US_LauraT_pp_harmonized.rds
# Trait_NOA[is.na(species) & is.na(genus), ]
Trait_NOA <- Trait_NOA[, .SD,
                       .SDcols = names(Trait_NOA) %like% "locom|feed|resp|volt|size|bf|species|genus|family|order"] %>%
  .[order %in% c(
    "Ephemeroptera",
    "Hemiptera",
    "Odonata",
    "Trichoptera",
    "Coleoptera",
    "Plecoptera",
    "Diptera"
  ), ]
Trait_NOA[, taxa := coalesce(species, genus, family, order)]
Trait_NOA[, c("feed_parasite",
              "unique_id") := NULL]

# Subset to those families that could be used for TPG delineation
Trait_NOA <- Trait_NOA[family %in% unique(trait_CONT[continent == "NOA", family]), ]
Trait_NOA_lf <- melt(
  Trait_NOA,
  id.vars = c("species",
              "genus",
              "family",
              "order",
              "taxa"),
  value.name = "affinity",
  variable.name = "trait"
)

## AUS ----
Trait_AUS <- trait_non_agg$Trait_AUS_harmonized.rds
Trait_AUS <- Trait_AUS[, .SD,
                       .SDcols = names(Trait_AUS) %like% "locom|feed|resp|volt|size|bf|species|genus|family|order"] %>%
  .[order %in% c(
    "Ephemeroptera",
    "Hemiptera",
    "Odonata",
    "Trichoptera",
    "Coleoptera",
    "Plecoptera",
    "Diptera"
  ), ]
Trait_AUS[, feed_parasite := NULL]
Trait_AUS[, taxa := coalesce(species, genus, family, order)]

# Add BF data
BF_AUS <- fread(file.path(
  data_in,
  "AUS_BF_missing_final.csv"
)) %>%
  na.omit(.)
BF_AUS[, bf_streamlined := as.double(bf_streamlined)]
Trait_AUS[BF_AUS,
          `:=`(
            bf_spherical = i.bf_spherical,
            bf_flattened = i.bf_flattened,
            bf_cylindrical = i.bf_cylindrical,
            bf_streamlined = i.bf_streamlined
          ),
          on = c(taxa = "family")]

# Add newly assigned traits
trait_assignment_LM_BK <- fread(file.path(data_in, "AUS_added_LM_BK.csv"))
trait_assignment_LM_BK <- trait_assignment_LM_BK[, .SD,
                                                 .SDcols = names(trait_assignment_LM_BK) %like%
                                                   "locom|feed|resp|volt|size|bf|species|genus|family|order"]
trait_assignment_LM_BK <- trait_assignment_LM_BK[family != "Coloburiscidae", ]
trait_assignment_LM_BK <- na.omit(trait_assignment_LM_BK)
trait_assignment_LM_BK[, taxa := family]
trait_assignment_LM_BK[,feed_parasite := NULL]
Trait_AUS <-
  rbind(Trait_AUS[!(is.na(species) &
                      is.na(genus) &
                      family %in% trait_assignment_LM_BK$family),],
        trait_assignment_LM_BK, fill = TRUE)
Trait_AUS <- Trait_AUS[family %in% unique(trait_CONT[continent == "AUS", family]), ]
Trait_AUS_lf <- melt(
  Trait_AUS,
  id.vars = c("species",
              "genus",
              "family",
              "order",
              "taxa"),
  value.name = "affinity",
  variable.name = "trait"
)


## NZ ----
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
  "Diptera"
), ]
Trait_NZ[, c("ovip_aqu",
             "ovip_ter",
             "ovip_ovo",
             "dev_holometabol",
             "dev_hemimetabol",
             "unique_id") := NULL]
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
    "size_medium",
    "size_small",
    "size_large",
    "bf_streamlined",
    "bf_flattened",
    "bf_cylindrical",
    "bf_spherical"
  )
)
Trait_NZ_lf <- melt(
  Trait_NZ,
  id.vars = c("species",
              "genus",
              "family",
              "order",
              "taxa"),
  value.name = "affinity",
  variable.name = "trait"
)





# test code with family dytiscidae (subset Trait_NZ_lf to dytisicidae before)
# Take family from NZ
# Trait_NZ[family == "Dytiscidae", ]
# # Search corrsponding TPG
# trait_CONT[continent == "NZ" & family == "Dytiscidae", group]
# Get mean trait profile for this group
# trait_CONT_mtpgs[continent == "NZ" & group == 6, ]
# test_dytiscidae[trait_CONT_mtpgs[tpg_id == "NZ_6", ], 
#                 mean_affinity := i.mean_affinity, 
#                 on = "trait"]
# test_dytiscidae[, grouping_feautre := sub("([a-z]{1,})(\\_)(.+)", "\\1", trait)]
# ggplot(test_dytiscidae) +
#  geom_boxplot(aes(x = trait, y = affinity)) +
#  geom_point(aes(x = trait, y = me3an_affinity), color = "red")

# Variability between traits, should be grouped 
# test_dytiscidae[, weight := .N, by = c("taxa", "grouping_feautre")]
# test_dytiscidae[, abs_diff_affinity := abs(affinity - mean_affinity)] 
# test_dytiscidae[, dist_to_mean_tp := weighted.mean(abs_diff_affinity, w = weight), 
#                 by = "taxa"]
# 
# ggplot(test_dytiscidae, aes(x = family, y = dist_to_mean_tp))+
#          geom_boxplot()


       
       
       