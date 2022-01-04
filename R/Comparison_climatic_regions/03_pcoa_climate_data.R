# __________________________________________________________________________________________________
# Comparison between climatic regions EU ----
# _________________________________________________________________________________________________

# _____________________________________________________________________
## EU ----
# total nr. of genera described (including data on species-lvl; for the traits we consider):
# 1457
# total nr. of genera with distributional information: 
# 353 (2422 taxa, mostly species-lvl)
# total nr. of genera restricted to cold and temperate: 
# 167 (cold 122, temperate 45)
# total nr. genera information available (cold and temperate): 
# 111 (cold 57, temperate 54)
# total nr. genera with complete information (for traits 
# size, feed, resp, lovom, volt; cold and temperate):
# 63 (cold 21, temperate 42)
# _____________________________________________________________________

# Read in trait data for genera from temperate and cold regions in Europe
traits_eu_clim <- readRDS(file = file.path(
  data_cache,
  "Cache_comparison_climatic_reg",
  "eu_glvl_murria.rds"
))

# Nr taxa per climateregion
traits_eu_clim[, .N, by = climateregion]

# Check completeness and subset only for complete data (temp preference not considered)
# completeness_trait_data(x = traits_eu_clim[, .SD, 
#                                            .SDcols = patterns("size|feed|resp|locom|volt|temp|ovip")])
traits_eu_clim <-
  na.omit(traits_eu_clim[, .SD, 
  .SDcols = patterns("feed|resp|volt|size|locom|ovip|climateregion|genus|family|order")])
str(traits_eu_clim)
traits_eu_clim[, .N, by = climateregion]

# Add label to climateregion column
traits_eu_clim[, continent := "EU"]
traits_eu_clim[, climateregion := paste0(continent, "_", climateregion)]
saveRDS(
  object = traits_eu_clim,
  file = file.path(
    data_cache,
    "Cache_comparison_climatic_reg",
    "eu_glvl_murria_pp.rds"
  )
)

### PCoA EU ----
cr_eu <- as.factor(traits_eu_clim$climateregion)
traits_eu_clim[, c("continent", "climateregion", "family", "order") := NULL]

# convert to df (dt does not support row names)
setDF(traits_eu_clim)
row.names(traits_eu_clim) <- traits_eu_clim$genus 
traits_eu_clim$genus <- NULL

# Convert to ktab object
vec <- sub("\\_.*", "\\1", names(traits_eu_clim))
blocks <- rle(vec)$lengths
traits_eu_clim <- prep.fuzzy(traits_eu_clim, blocks)
traits_eu_clim <- ktab.list.df(list(traits_eu_clim))
eu_clim_dist <- dist.ktab(traits_eu_clim, type = "F") 

# PcOA
eu_clim_pcoa <- dudi.pco(eu_clim_dist, scannf = FALSE)
summary(eu_clim_pcoa)

eu_clim_pcoa_scores <- eu_clim_pcoa$li[, 1:2]
eu_clim_pcoa_scores$genus <- rownames(eu_clim_pcoa_scores)
setDT(eu_clim_pcoa_scores)
eu_clim_pcoa_scores[, climateregion := cr_eu]
eu_clim_pcoa_scores[, climateregion := as.factor(climateregion)]

hull_eu <- eu_clim_pcoa_scores %>%
  group_by(climateregion) %>%
  slice(chull(A1, A2))
# hull_eu$climateregion %>% unique()

cbb_palette <-
  c(
    "#000000",
    "#f04a18",
    "#34a4e6",
    "#ad1669"
  )
ggplot(eu_clim_pcoa_scores, aes(x = A1, y = A2)) +
  geom_point(alpha = 0.2) +
  geom_polygon(
    data = hull_eu, 
    alpha = 0.4,
    aes(fill = climateregion)
  ) +
  scale_fill_manual(
    values = cbb_palette,
    name = "Climate region"
  ) +
  labs(
    x = paste0(
      "Axis 1",
      " (",
      round(
        eu_clim_pcoa$eig[1] / sum(eu_clim_pcoa$eig) * 100,
        digits = 2
      ),
      "%)"
    ),
    y = paste0("Axis 2", " (", round(
      eu_clim_pcoa$eig[2] / sum(eu_clim_pcoa$eig) * 100,
      digits = 2
    ), "%)")
  ) +
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
    panel.grid = element_blank()
  )

# "Quality" of trait space
# Q <- coranking(Xi = clim_dist, 
#                X = clim_pcoa$li[, 1:2]) 
# NX <- coRanking::R_NX(Q)
# coRanking::AUC_ln_K(NX)

### Permanova EU ----
eu_perm_clim <- adonis2(eu_clim_dist ~ cr_eu)

# Find no evidence that regions differ in distances based on their trait profiles
# 1.9 of the total variance (in the distance matrix) is explained
# by the factor climateregion
eu_perm_clim

# Test assumption of multivariate normal dispersion
# looks relatively similar
eu_bd <- betadisper(eu_clim_dist,
                    cr_eu)

# Permutation test
# Find no evidence that multivariate normal dispersions of cold and temperate regions
# is different -> assumptions is met
permutest(eu_bd)


# __________________________________________________________________________________________________
## NoA ----
# total nr. of genera described (including data on species-lvl, for the traits we consider):
# 1330
# total nr. of genera with distributional information: 
# 626
# total nr. of genera restricted to cold and temperate: 
# 167 (cold 122, temperate 45)
# total nr. genera information available (cold and temperate): 
# 134 (cold 111, temperate 23)
# total nr. genera with complete information (for traits 
# size, feed, resp, lovom, volt; cold and temperate):
# 91 (cold 84, temperate 7)
# __________________________________________________________________________________________________
traits_noa_clim <- readRDS(file = file.path(data_cache,
                                            "Cache_comparison_climatic_reg",
                                            "noa_glvl.rds"))

# Much more taxa in cold climate region 
traits_noa_clim[, .N, by = climateregion]
traits_noa_clim <- traits_noa_clim[!climateregion %in% c("polar", "arid"), ]

# Check completeness and subset only for complete data (temp preference not considered)
# completeness_trait_data(x = traits_noa_clim[, .SD,
#                                             .SDcols = patterns("size|feed|resp|locom|volt|ovip")])
traits_noa_clim <-
  na.omit(traits_noa_clim[, .SD,
                          .SDcols = patterns("feed|size|locom|resp|volt|climateregion|genus|family|order")])
traits_noa_clim[, .N, by = climateregion]

# _________________________________________________________________________________________
### Get those genera where trait information is missing ----
# Missing information
noa_genera_temp <- readRDS(file.path(data_cache, "noa_genera_temp.rds"))
genera_missing <- noa_genera_temp[!Genus %in% traits_noa_clim$genus, ]

# traits_noa_clim <- readRDS(file = file.path(data_cache,
#                                             "Cache_comparison_climatic_reg",
#                                             "noa_glvl.rds"))


trait_NOA <- readRDS(file.path(data_in, "Traits_US_LauraT_pp_harmonized.rds"))
genera_missing <- merge(genera_missing,
                        trait_NOA[is.na(species), ],
                        by.x = "Genus",
                        by.y = "genus",
                        all.x = TRUE)
genera_missing[, c("unique_id", "species") := NULL]

# Add taxonomic information
library(taxize)
genera_ids <- get_ids(sci_com = genera_missing[is.na(family), Genus], 
        db = "gbif",
        rows = 1)
genera_gbif <- data.table(id_genus = genera_ids$gbif, 
                          genus = names(genera_ids$gbif))
genera_gbif <- cbind(classification(id = genera_gbif[!is.na(id_genus), id_genus],
                                    db = "gbif"))[c("genus", "family", "order")]
setDT(genera_gbif)
genera_missing[genera_gbif,
               `:=`(family = i.family,
                    order = i.order),
               on = c(Genus = "genus")]
setnames(genera_missing, "Genus", "genus")
cols <-
  grep("genus|family|order|feed.*|locom.*|size.*|resp.*|volt.*",
       names(genera_missing), 
       value = TRUE)
setcolorder(genera_missing,
            c("genus",
              "family",
              "order",
              sort(cols[!cols %in% c("genus", "family", "order")])))
fwrite(genera_missing[, ..cols],
       file.path(data_out, "missing_genera_temperate.csv"))
# _________________________________________________________________________________________

# Add label to climateregion column
traits_noa_clim[, continent := "NOA"]
traits_noa_clim[, climateregion := paste0(continent, "_", climateregion)]
saveRDS(
  object = traits_noa_clim,
  file = file.path(
    data_cache,
    "Cache_comparison_climatic_reg",
    "noa_glvl_pp.rds"
  )
)


### PCoA NOA ----
cr_noa <- as.factor(traits_noa_clim$climateregion)
traits_noa_clim[, c("continent", "climateregion", "family", "order") := NULL]

# convert to df (dt does not support row names)
setDF(traits_noa_clim)
row.names(traits_noa_clim) <- traits_noa_clim$genus 
traits_noa_clim$genus <- NULL

# Convert to ktab object
vec <- sub("\\_.*", "\\1", names(traits_noa_clim))
blocks <- rle(vec)$lengths
traits_noa_clim <- prep.fuzzy(traits_noa_clim, blocks)
traits_noa_clim <- ktab.list.df(list(traits_noa_clim))
noa_clim_dist <- dist.ktab(traits_noa_clim, type = "F") 

# PcOA
noa_clim_pcoa <- dudi.pco(noa_clim_dist, scannf = FALSE)
summary(noa_clim_pcoa)


### Permanova NoA ----
noa_perm_clim <- adonis2(noa_clim_dist ~ cr_noa)

# Find no evidence that regions differ in distances based on their trait profiles
# 1.1 of the total variance (in the distance matrix) is explained
# by the factor climateregion
noa_perm_clim

# Test assumption of multivariate normal dispersion
# looks relatively similar
noa_bd <- betadisper(noa_clim_dist,
                     cr_noa)

# Permutation test
# Find no evidence that multivariate normal dispersions of cold and temperate regions
# is different -> assumptions is met
permutest(noa_bd)


# _____________________________________________________________________
# Combined data from NOA and EU ----
# _____________________________________________________________________

# Prep
# no ovip in NOA, thus rm in EU data
traits_clim <- rbind(traits_eu_clim[, -c("ovip_aqu", "ovip_ovo", "ovip_ter")],
                     traits_noa_clim)

# genera not unique
# traits_clim[duplicated(genus), ]
traits_clim[, genus := paste0(continent, "_", genus)]

# extra dataset with cr and continent factory
cr_continents <- traits_clim[, .(genus, climateregion, continent)]
traits_clim[, c("continent", "climateregion", "family", "order") := NULL]

# convert to df (dt does not support row names)
setDF(traits_clim)
row.names(traits_clim) <- traits_clim$genus 
traits_clim$genus <- NULL

# Convert to ktab object
vec <- sub("\\_.*", "\\1", names(traits_clim))
blocks <- rle(vec)$lengths
traits_clim <- prep.fuzzy(traits_clim, blocks)
traits_clim <- ktab.list.df(list(traits_clim))
clim_dist <- dist.ktab(traits_clim, type = "F") 


## PCoA ----
clim_pcoa <- dudi.pco(clim_dist, scannf = FALSE)
summary(clim_pcoa)

# Factor map with colors for the different climateregions
# s.class(clim_pcoa$li,
#         climateregions,
#         col = TRUE)

# Plot scores and hulls for the different climate regions
clim_pcoa_scores <- clim_pcoa$li[, 1:2]
clim_pcoa_scores$genus <- rownames(clim_pcoa_scores)
setDT(clim_pcoa_scores)
clim_pcoa_scores[cr_continents,
                 `:=`(climateregion = i.climateregion,
                      continent = i.continent),
                 on = "genus"]
clim_pcoa_scores[, `:=`(continent = as.factor(continent),
                        climateregion = as.factor(climateregion))]

# Hulls for each cr/continent combination
hull <- clim_pcoa_scores %>%
  group_by(climateregion) %>%
  slice(chull(A1, A2))
hull$climateregion %>% unique()

cbb_palette <-
  c(
    "#000000",
    "#f04a18",
    "#34a4e6",
    "#ad1669"
  )
ggplot(clim_pcoa_scores, aes(x = A1, y = A2)) +
  geom_point(alpha = 0.2) +
  geom_polygon(
    data = hull, 
    alpha = 0.4,
    aes(fill = climateregion)
  ) +
  scale_fill_manual(
    values = cbb_palette,
    name = "Climate region"
  ) +
  labs(
    x = paste0(
      "Axis 1",
      " (",
      round(
        clim_pcoa$eig[1] / sum(clim_pcoa$eig) * 100,
        digits = 2
      ),
      "%)"
    ),
    y = paste0("Axis 2", " (", round(
      clim_pcoa$eig[2] / sum(clim_pcoa$eig) * 100,
      digits = 2
    ), "%)")
  ) +
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
    panel.grid = element_blank()
  )
ggsave(
  filename = file.path(
    data_paper,
    "Graphs",
    "PCOA_EU_climateregion.png"
  ),
  width = 35,
  height = 20,
  units = "cm"
)

# "Quality" of trait space
Q <- coranking(Xi = clim_dist, 
          X = clim_pcoa$li[, 1:2]) 
NX <- coRanking::R_NX(Q)
coRanking::AUC_ln_K(NX)


## Combined ----

# TODO: Returns error for some reason, need to fix this!
# (have there been changes in vegan?)
cr_continents[, `:=`(
  climateregion = as.factor(climateregion),
  continent = as.factor(continent)
)]

# Shuffeling design with continent and climateregion
# perm <- how(
#   within = Within(type = "free"),
#   plots = Plots(strata = cr_continents$climateregion),
#   blocks = cr_continents$continent
# )

# Shuffeling just within continents
perm <- how(
   within = Within(type = "free"),
   blocks = cr_continents$continent
 )
# sequential test
perm_clim <-
  adonis2(clim_dist ~ continent * climateregion,
          data = cr_continents,
          permutations = perm)
perm_clim

# test interaction continent:climateregion
# if there's some interaction reject convergence 
# of within-region variability
perm_clim_margin <-
  adonis2(clim_dist ~ continent:climateregion,
          data = cr_continents,
          permutations = perm,
          by = "margin")

# Assumption of multivariate normal dispersion not met
bd <- betadisper(clim_dist,
                 cr_continents$climateregion)
# Permutation test
permutest(bd)



