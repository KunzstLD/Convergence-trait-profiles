# __________________________________________________________________________________________________
### Variation between mean trait profiles of TPGs ----
# __________________________________________________________________________________________________
trait_CONT <-
  readRDS(file.path(data_cache, "trait_dat_grp_assig.rds"))
trait_CONT_mtpgs <- trait_CONT[, .(continent,
                                   mean_affinity = mean(affinity)),
                               by = c("continent" , "group", "trait")]
trait_CONT_mtpgs <- trait_CONT_mtpgs[, .(continent,
                                         group,
                                         trait,
                                         mean_affinity)] %>%
  unique(.) %>%
  dcast(., ... ~ trait, value.var = "mean_affinity")
trait_CONT_mtpgs[, tpg_id := paste0(continent, "_", group)]
normalize_by_rowSum(trait_CONT_mtpgs,
                    non_trait_cols = c("continent",
                                       "group",
                                       "tpg_id"))
dist_between_tpgs_aff <- list()
for (i in unique(trait_CONT_mtpgs$continent)) {
  dat <- trait_CONT_mtpgs[continent == i,]
  setDF(dat)
  rownames(dat) <- dat$tpg_id
  dat$tpg_id <- NULL
  dat$group <- NULL
  dat$continent <- NULL
  
  vec <- sub("\\_.*", "\\1", names(dat))
  blocks <- rle(vec)$lengths
  dat <- prep.fuzzy(dat, blocks)
  dat <- ktab.list.df(list(dat))
  dist_between_tpgs_aff[[i]] <-
    as.vector(dist.ktab(dat, type = "F"))
}
dist_between_tpgs_aff <-
  lapply(dist_between_tpgs_aff, as.data.table) %>%
  rbindlist(., idcol = "continent")
setnames(dist_between_tpgs_aff, "V1", "dist")
dist_between_tpgs_aff[, `:=`(mean = mean(dist),
                             median = median(dist),
                             sd = sd(dist)), by = "continent"]

# Range
dist_between_tpgs_aff[, range(dist), by = "continent"] %>% 
  .[, diff(V1), by = "continent"]

# Add kg zones
# KG Numbers: NZ 9, SA 14, AUS 17, EUR 18, NA 22
kg <- data.table(
  continent = c("AUS",
                "EU",
                "NOA",
                "NZ",
                "SA"),
  kg_zones = c(17,
               18,
               20,
               9,
               14)
)
dist_between_tpgs_aff[kg, kg_zones := i.kg_zones, on = "continent"]

## LM ----
lm_btw_tpgs <- lm(dist ~ kg_zones,
                  data = dist_between_tpgs_aff)
summary(lm_btw_tpgs) 
par(mfrow=c(2,2))
plot(lm_btw_tpgs) # Okay?
estimate <- round(coef(lm_btw_tpgs)["kg_zones"], digits = 3)
p_val <-
  round(summary(lm_btw_tpgs)$coefficients[, "Pr(>|t|)"]["kg_zones"], digits = 3)

# GLM
# glm_btw_tpgs <- glm(dist ~ kg_zones,
#                     data = dist_between_tpgs_aff,
#                     family = quasibinomial)
# summary(glm_btw_tpgs)
# null_glm_btw_tpgs <- glm(dist ~ 1,
#                          data = dist_between_tpgs_aff,
#                          family = quasibinomial)
# 
# # We can use the anova test and specify
# # the chisquare test for comparing the deviances
# # technically, this is a Likelihood ratio test
# anova(glm_btw_tpgs, null_glm_btw_tpgs, test = "Chisq")
# 
# # If we divide the Residual Deviance by the
# # degrees of freedom we yield a dispersion parameter
# # of approximately 1 which is very good
# # (remember that >2 or <0.5 would indicate over- or underdispersion)
# glm_btw_tpgs$deviance / glm_btw_tpgs$df.resid
# # Still underdispersion with quasibinomial model
# 
# # QQ plot
# qr <- qres.binom(glm_tvar_env)
# # we compute Dunn-Smyth-Residuals that are most appropriate for GLMs
# qqnorm(qr)

# Plot
ggplot(dist_between_tpgs_aff,
       aes(x = kg_zones, y = dist)) +
  geom_point(aes(color = continent),
             position = position_jitter(width = 0.1, seed = 1234)) +
  geom_smooth(method = "lm",#glm
              # method.args = list(family = "quasibinomial"),
              formula = y ~ x,
              se = FALSE,
              color = "black") +
  scale_color_d3(name = "Region",
                 labels = c("AUS", "EUR", "NA", "NZ", "SAf")) +
  scale_x_continuous(breaks = c(9, 14, 17, 18, 20)) + 
  stat_summary(fun = "mean", 
               color = "grey",
               geom = "crossbar") +
  labs(x = "Number of Köppen Geiger zones", 
       y = "Distance between mean trait profiles of TPGs") +
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
  ) +
  annotate(
    "text",
    x = 12,
    y = 0.58,
    label = paste0("Estimate = ", estimate, "\np-value = ", p_val),
    color = "black",
    size = 4,
    fontface = "bold"
  ) # +
  # annotate(
  #   "curve",
  #   x = 22.05,
  #   y = 0.426,
  #   xend = 22,
  #   yend = 0.36,
  #   curvature = 0.2,
  #   arrow = arrow(length = unit(2, "mm"))
  # ) +
  # annotate(
  #   "text",
  #   x = 21.5,
  #   y = 0.35,
  #   label = "TPG_NA4 vs. TPG_NA5"
  # )
ggsave(
  filename = file.path(
    data_paper,
    "Graphs",
    "Variation_between_tpgs_vs_kg.png"
  ),
  width = 35,
  height = 20,
  units = "cm"
)

# Plot, arranged according to the number of Köppen-Geiger zones
# TODO: Add a line/arrow in the plot indicating the direction of environmental heterogeneity
# and the number of Köppen-Geiger climate zones
# dist_between_tpgs_aff[, continent := factor(continent,
#                                             levels = c("NZ", "SA", "AUS", "EU", "NOA"))]
# set.seed(1234)
# ggplot(dist_between_tpgs_aff, aes(x = continent, y = dist)) +
#   geom_violin() +
#   geom_point(position = position_jitter(width = 0.05, seed = 1234)) +
#   stat_summary(fun = "mean", color = "red") +
#   labs(x = "Continent or region",
#        y = "Distance between mean trait profiles of TPGs") +
#   scale_x_discrete(labels = c("NZ", "SA", "AUS", "EUR", "NA")) +
#   theme_bw() +
#   theme(
#     axis.title = element_text(size = 16),
#     axis.text.x = element_text(family = "Roboto Mono",
#                                size = 14),
#     axis.text.y = element_text(family = "Roboto Mono",
#                                size = 14),
#     legend.title = element_text(family = "Roboto Mono",
#                                 size = 16),
#     legend.text = element_text(family = "Roboto Mono",
#                                size = 14),
#     strip.text = element_text(family = "Roboto Mono",
#                               size = 14)
#   ) +
#   annotate(
#     "segment",
#     x = 0.5,
#     y = 0,
#     xend = 5.5,
#     yend = 0,
#     size = 5.5,
#     linejoin = "mitre",
#     arrow = arrow(type = "closed", length = unit(0.01, "npc"))
#   ) +
#   annotate(
#     "text",
#     x = 3,
#     y = 0,
#     label = "Environmental heterogeneity",
#     color = "white",
#     size = 4,
#     fontface = "bold"
#   ) +
#   annotate(
#     "text",
#     x = 1,
#     y = 0.08,
#     label = "KG-zones: 9",
#     color = "black",
#     size = 4,
#     fontface = "bold"
#   ) +
#   annotate(
#     "text",
#     x = 2,
#     y = 0.08,
#     label = "KG-zones: 14",
#     color = "black",
#     size = 4,
#     fontface = "bold"
#   ) +
#   annotate(
#     "text",
#     x = 3,
#     y = 0.08,
#     label = "KG-zones: 17",
#     color = "black",
#     size = 4,
#     fontface = "bold"
#   ) +
#   annotate(
#     "text",
#     x = 4,
#     y = 0.08,
#     label = "KG-zones: 18",
#     color = "black",
#     size = 4,
#     fontface = "bold"
#   ) +
#   annotate(
#     "text",
#     x = 5,
#     y = 0.08,
#     label = "KG-zones: 22",
#     color = "black",
#     size = 4,
#     fontface = "bold"
#   ) +
#   annotate(
#     "curve",
#     x = 5.05,
#     y = 0.426,
#     xend = 4.8,
#     yend = 0.36,
#     curvature = 0.2,
#     arrow = arrow(length = unit(2, "mm"))
#   ) +
#   annotate("text",
#            x = 4.8,
#            y = 0.35,
#            label = "TPG_NA4 vs. TPG_NA5")
# ggsave(
#   filename = file.path(
#     data_paper,
#     "Graphs",
#     "Variation_between_tpgs.png"
#   ),
#   width = 40,
#   height = 20,
#   units = "cm"
# )

