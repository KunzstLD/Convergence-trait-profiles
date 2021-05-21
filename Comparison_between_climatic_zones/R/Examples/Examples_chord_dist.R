library(dplyr)
library(vegan)

# Example data
ex_trait <- data.frame(
  feed_shredder = c(0, NA, 0.7, 0),
  feed_gatherer = c(0.1, NA, 0.3, 0),
  feed_predator = c(0, NA, 0, 0),
  feed_parasite = c(0, NA, 0, 0),
  feed_filter = c(0, NA, 0, 0),
  feed_herbivore = c(0.9, NA, 0, 1),
  resp_teg = c(NA, 0.5, 0, 0),
  resp_gil = c(NA, 0.5, 0, 0.5 ), 
  
  resp_pls_spi = c(NA, 0, 1, 0.5), 
  
  volt_semi = c(NA, 0, 0, 0),
  volt_uni = c(NA, 0, 1, 1),
  volt_bi_multi = c(NA, 1, 0, 0))

# convert NAs to zero
ex_trait_zero <- vapply(ex_trait, function(y)
  ifelse(is.na(y), 0, y), FUN.VALUE = c(
    "1" = 0,
    "2" = 0,
    "3" = 0,
    "4" = 0
  ))

# Orloci's chord distance
dist_mat <- decostand(ex_trait, "norm", na.rm = TRUE) %>%
  vegdist(., "euclidean", na.rm = TRUE) %>%
  as.matrix()

# Example calculation by hand Orloci chord distance
sqrt(
  (0/sqrt(0.1^2+0.9^2)-0/sqrt(1^2+0.5^2+0.5^2+1^2))^2 +
  (0.1/sqrt(0.1^2+0.9^2)-0/sqrt(1^2+0.5^2+0.5^2+1^2))^2 +
  (0/sqrt(0.1^2+0.9^2)-0/sqrt(1^2+0.5^2+0.5^2+1^2))^2 +
  (0/sqrt(0.1^2+0.9^2)-0/sqrt(1^2+0.5^2+0.5^2+1^2))^2 +
  (0/sqrt(0.1^2+0.9^2)-0/sqrt(1^2+0.5^2+0.5^2+1^2))^2 +
  (0.9/sqrt(0.1^2+0.9^2)-1/sqrt(1^2+0.5^2+0.5^2+1^2))^2
)/sqrt(2)

# Standardized to range [0, 1]
dist_mat <- dist_mat/sqrt(2)

# With zero data
dist_mat_zero <- decostand(ex_trait_zero, "norm", na.rm = TRUE) %>%
  vegdist(., "euclidean", na.rm = TRUE) %>%
  as.matrix()
dist_mat_zero/sqrt(2)

# Using ade4 dist functions!
library(ade4)

vec <- sub("\\_.*", "\\1", names(ex_trait))
blocks <- rle(vec)$lengths
ex_trait <- prep.fuzzy(ex_trait, blocks)
ktab_ex_trait <- ktab.list.df(list(ex_trait))

dist_by_ktab <- dist.ktab(ktab_ex_trait, type = "F")
as.matrix(dist_by_ktab)

# glocor gives correlations between distances obtained for each trait and the global distances
# obtained by mixing all the traits (= contribution of traits to the global distances)
kdist.cor(ktab_ex_trait, type = "F")




data(bsetal97)

w <- prep.fuzzy(bsetal97$biol, bsetal97$biol.blo)
w[1:6, 1:10]
ktab1 <- ktab.list.df(list(w))

dis <- dist.ktab(ktab1, type = "F")
kdist.cor(ktab1, type = "F")
str(w)