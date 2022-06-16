# __________________________________________________________________________________________________
# PcOA ----
# TODO: check scripts for function to find suitable nr. of axis
# https://github.com/tanogc/overarching_functional_space
# __________________________________________________________________________________________________
trait_data_bind <- readRDS(file.path(
  data_cache,
  "trait_data_ww_bind.rds"
))

# Preproc
trait_data_bind[, id := paste0(continent, "_", family)]
continent_names <- factor(trait_data_bind$continent)
trait_data_bind_cp <- copy(trait_data_bind)
trait_data_bind[, c("family", "order", "continent") := NULL]
setDF(trait_data_bind)

# Add row.names
taxa_reg_names <- trait_data_bind$id
row.names(trait_data_bind) <- taxa_reg_names
trait_data_bind$id <- NULL

# Convert to ktab object
vec <- sub("\\_.*", "\\1", names(trait_data_bind))
blocks <- rle(vec)$lengths
trait_data_bind <- prep.fuzzy(trait_data_bind, blocks)
trait_data_bind <- ktab.list.df(list(trait_data_bind))
comb_dist <- dist.ktab(trait_data_bind, type = "F")
is.euclid(comb_dist)

# PcOA
comb_pcoa <- dudi.pco(comb_dist, scannf = FALSE, nf = 25)
summary(comb_pcoa)

# Explained variance
# 10 axis explain more than 80 %
# 5 axis around 57 %
sum(comb_pcoa$eig[1:10])/sum(comb_pcoa$eig)
sum(comb_pcoa$eig[1:5])/sum(comb_pcoa$eig)

# Factor map with colors for the different continents
# s.label(comb_pcoa$li, clabel = 0.7) # Fig. 2 in the main text
# s.class(comb_pcoa$li,
#   fac = continent_names,
#   col = TRUE
# )

## Quality of trait space ----
# with Guido Kraemer package
Q <- coranking(Xi = comb_dist, 
               X = comb_pcoa$li[, 1:2]) 
NX <- coRanking::R_NX(Q)
coRanking::AUC_ln_K(NX)


## Plotting ----
# Add loads? Eigenvector * sq.root eigenvalue
# ? comb_pcoa$l1
# comb_pcoa$tab

# Plot scores and hulls for the different continents
pcoa_scores <- comb_pcoa$li[1:2]
pcoa_scores$id <- rownames(pcoa_scores)
setDT(pcoa_scores)
pcoa_scores[, continent := factor(sub("(\\w)(\\_)(.+)", "\\1", id))]
saveRDS(pcoa_scores, file = file.path(data_cache, "pcoa_scores.rds"))

# Hulls for each continent
hull <- pcoa_scores %>%
  group_by(continent) %>%
  slice(chull(A1, A2))
saveRDS(hull, file.path(data_cache, "hull_pcoa.rds"))

### Tait space plot ----
pcoa_plot <- ggplot(pcoa_scores, aes(x = A1, y = A2)) +
  geom_point(alpha = 0.15) +
  geom_polygon(data = hull, alpha = 0.05,
               aes(color = continent)) +
  scale_color_d3(name = "Continent",
                 labels = c("AUS", "EUR", "NA", "NZ", "SA")) +
  labs(x = paste0("Axis 1",
                  " (",
                  round(
                    comb_pcoa$eig[1] / sum(comb_pcoa$eig) * 100,
                    digits = 2
                  ),
                  "%)"),
       y = paste0("Axis 2", " (", round(
         comb_pcoa$eig[2] / sum(comb_pcoa$eig) * 100,
         digits = 2
       ), "%)")) +
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
    panel.grid = element_blank()
  )
ggsave(
  filename = file.path(
    data_paper,
    "Graphs",
    "PCOA_continent.png"
  ),
  width = 35,
  height = 20,
  units = "cm"
)

# Plot those taxa from Australia that do not overlap with the other trait spaces (on the first two PCOA axis) 
pcoa_plot + geom_text_repel(
  data = function(x) {
    x[continent == "AUS", ]
  },
  aes(label = id)
)
pcoa_scores[, family := sub("(.+)\\_", "", id)]

pcoa_plot + geom_text_repel(
  data = function(x) {
    x[continent == "AUS" & family %in% c(
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
    ), ]
  },
  aes(label = family)
)
ggsave(
  filename = file.path(
    data_paper,
    "Graphs",
    "PCOA_continent_AUS_outliers.png"
  ),
  width = 35,
  height = 20,
  units = "cm"
)

### Contour/density plot ----
# Calculate with kde to get exact contour lines
kde_traits <- kde(pcoa_scores[, .(A1, A2)])
plot(kde_traits, cont = c(25, 50, 75))

# Get level for 75 % and 50 % (i.e. 25 % contour line and 50 % contour line)
level_75 <- kde_traits$cont[names(kde_traits$cont) == "75%"]
level_50 <- kde_traits$cont[names(kde_traits$cont) == "50%"]
kde_traits$cont[names(kde_traits$cont) == "25%"]

# Identify high density regions/points
estimate <- kde(pcoa_scores[, .(A1, A2)], eval.points = pcoa_scores[, .(A1, A2)])$estimate

# ~ 2.8 is the 25 % contour line (i.e. smallest region that contains 25 % probability mass) in the plot
contour_25 <- as.numeric(names(estimate[estimate >= level_75]))
contour_50 <- as.numeric(names(estimate[estimate >= level_50]))
saveRDS(contour_25,
        file.path(data_cache, "contour_25_TS.rds"))
saveRDS(contour_50,
        file.path(data_cache, "contour_50_TS.rds"))

pl_contour <- ggplot(pcoa_scores, aes(x = A1, y = A2)) +
  stat_density_2d(aes(color = ..level..),
                  #bins = 4,
                  breaks = c(1.72, 2.3, 2.86)) +
  geom_point(alpha = 0.15) +
  scale_color_viridis_c() +
  labs(x = paste0("Axis 1",
                  " (",
                  round(
                    comb_pcoa$eig[1] / sum(comb_pcoa$eig) * 100,
                    digits = 2
                  ),
                  "%)"),
       y = paste0("Axis 2", " (", round(
         comb_pcoa$eig[2] / sum(comb_pcoa$eig) * 100,
         digits = 2
       ), "%)"),
       color = "Level") +
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
    panel.grid = element_blank()
  )
# Interactive view with plotly
ggplotly(pl_contour)
ggsave(
  filename = file.path(
    data_paper,
    "Graphs",
    "Contour.png"
  ),
  width = 35,
  height = 20,
  units = "cm"
)

# __________________________________________________________________________________________________
# PERMANOVA ----
# __________________________________________________________________________________________________

# Is strata/permutations needed?
# perm <- how(
#    within = Within(type = "free"),
#    plots = Plots(strata = continent_names)
# #   blocks = continent_names
#  )
comb_perm <- adonis2(comb_dist ~ continent_names)
comb_perm
# 4.95 % of the total variance (in the distance matrix)
# is explained by the factor continent

# Assumptions of PERMANOVA do not hold
bd_comb <- betadisper(comb_dist,
  group = continent_names
)
boxplot(bd_comb)
permutest(bd_comb)
saveRDS(bd_comb,
        file.path(data_cache, "betadisper_continental.rds"))

# __________________________________________________________________________________________________
# Multivariate Welch Test ----
# __________________________________________________________________________________________________

# Wd statistic
WdS(dm = comb_dist, f = continent_names)
# permutation test
WdS.test(dm = comb_dist, f = continent_names)


