# __________________________________________________________________________________________________
#### Comparison between climatic regions ####
# __________________________________________________________________________________________________

# Read in trait data for genera from temperate and cold regions in Europe
traits_eu_clim <- readRDS(file = file.path(
  data_cache,
  "Cache_comparison_climatic_reg",
  "eu_glvl_murria.rds"
))

# Nr observations per climateregion
traits_eu_clim[, .N, by = climateregion]

completeness_trait_data(x = traits_eu_clim[, .SD, .SDcols = patterns("size|feed|resp|locom|volt|temp|ovip")])
# Subset only for complete data (temp preference not considered)
traits_eu_clim <-
  na.omit(traits_eu_clim[, .SD, .SDcols = patterns("feed|resp|volt|size|locom|ovip|climateregion|genus|family|order")])

# __________________________________________________________________________________________________
# PCoA ----
# __________________________________________________________________________________________________
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
clim_pcoa <- dudi.pco(clim_dist, scannf = FALSE, nf = 3)
# add color with groups
scatter(clim_pcoa)
summary(clim_pcoa)

# TODO: After between class or discriminant analysis only 1 axis remains. Why?
# Between class analysis
bet_climregions <- bca(clim_pcoa,
                       climateregions,
                       scannf = FALSE)
summary(bet_climregions)
library(adegraphics)
plot(bet_climregions, 
     row.pellipses.col = adegpar()$ppalette$quali(2))

# 
(rt_betclimregions <- randtest(bet_climregions))

dis_climregions <- discrimin(clim_pcoa, 
                             climateregions,
                             scannf = FALSE)
plot(dis_climregions, row.pellipses.col = adegpar()$ppalette$quali(2))

