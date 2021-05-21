# _____________________________________________________________________
#### Comparison between climatic regions EU ####
# _____________________________________________________________________

# _______
## EU ##
# _______

# Read in trait data for genera from temperate and cold regions in Europe
traits_eu_clim <- readRDS(file = file.path(
  data_cache,
  "Cache_comparison_climatic_reg",
  "eu_glvl_murria.rds"
))

# Nr observations per climateregion
traits_eu_clim[, .N, by = climateregion]

# Check completeness and subset only for complete data (temp preference not considered)
completeness_trait_data(x = traits_eu_clim[, .SD, .SDcols = patterns("size|feed|resp|locom|volt|temp|ovip")])
traits_eu_clim <-
  na.omit(traits_eu_clim[, .SD, 
  .SDcols = patterns("feed|resp|volt|size|locom|ovip|climateregion|genus|family|order")])
str(traits_eu_clim)
traits_eu_clim[, .N, by = climateregion]

# _______
## NoA ##
# _______



# _____________________________________________________________________
# PCoA ----
# _____________________________________________________________________
climateregions <- traits_eu_clim$climateregion
climateregions <- as.factor(traits_eu_clim$climateregion)
traits_eu_clim[, c("family", "order", "climateregion") := NULL]
setDF(traits_eu_clim)
row.names(traits_eu_clim) <- traits_eu_clim$genus 
traits_eu_clim$genus <- NULL

# Convert to ktab object
vec <- sub("\\_.*", "\\1", names(traits_eu_clim))
blocks <- rle(vec)$lengths
traits_eu_clim <- prep.fuzzy(traits_eu_clim, blocks)
traits_eu_clim <- ktab.list.df(list(traits_eu_clim))
clim_dist <- dist.ktab(traits_eu_clim, type = "F") 

# PcOA
clim_pcoa <- dudi.pco(clim_dist, scannf = FALSE, nf = 21)
summary(clim_pcoa)

# Factor map with colors for the different climateregions
# s.class(clim_pcoa$li,
#         climateregions,
#         col = TRUE)

# Plot scores and hulls for the different climate regions
clim_pcoa_scores <- clim_pcoa$li[, 1:2]
clim_pcoa_scores$genus <- rownames(clim_pcoa_scores)
setDT(clim_pcoa_scores)
clim_pcoa_scores[, cr := climateregions]

# Hulls for each continent
hull <- clim_pcoa_scores %>%
  group_by(cr) %>%
  slice(chull(A1, A2))

cbb_palette <-
  c(
    "#000000",
    "#f04a18"
  )
ggplot(clim_pcoa_scores, aes(x = A1, y = A2, fill = cr)) +
  geom_point(alpha = 0.2) +
  geom_polygon(
    data = hull, alpha = 0.4,
    aes(fill = cr)
  ) +
  scale_fill_manual(
    values = cbb_palette,
    name = "Climate region",
    labels = c("Cold", "Temperate")
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

# _____________________________________________________________________
# Permanova
# _____________________________________________________________________
# TODO: Returns error for some reason, need to fix this!
# (have there been changes in vegan?)
test <- as.matrix(clim_dist)
perm_clim <- adonis(test ~ climateregions)

# Find no evidence that regions differ in distances based
# on their trait profiles
# 1.7 of the total variance (in the distance matrix) is explained
# by the factor climateregion
perm_clim
 
# Test assumption of multivariate normal dispersion
# looks relatively similar
bd <- betadisper(clim_dist,
           climateregions)

# Permutation test
# Find no evidence that multivariate normal dispersions of cold and temperate regions
# is different -> assumptions is met
permutest(bd)

