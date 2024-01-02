# __________________________________________________________________________________________________
# Comparing family trait profiles with the number of Köppen Geiger zones ----
# Calc. a distance between trait profiles of families and their distance to the 
# mean trait profile of their tpgs 
# __________________________________________________________________________________________________

## Data prep. ----
# Load trait data & groups, calc. mean trait profile
trait_CONT <- readRDS(file.path(data_cache, "trait_dat_grp_assig.rds"))
trait_CONT[, mean_affinity := mean(affinity),
           by = c("continent" , "group", "trait")]

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

# Calc distance to mean affinity
trait_CONT[, overlap_per_trait := fun_overlap_ind(p = affinity,
                             q = mean_affinity),
           by = c("continent",
                  "group",
                  "family",
                  "grouping_feature")]
trait_CONT[, dist_to_mean_tp := sqrt(mean(overlap_per_trait)), by = c("continent",
                                                                      "group",
                                                                      "family")] 
# Test
# trait_CONT[continent == "AUS" &
#              group == 1, fun_overlap_ind(p = affinity,
#                                          q = mean_affinity),
#            by = c("family", "grouping_feature")] %>% 
#   .[, sqrt(mean(V1)), by = "family"]

# Create dataset that only contains distances to mean trait profile
dist_within_tpgs <-
  unique(trait_CONT[, .(continent, family, group, dist_to_mean_tp)])
dist_within_tpgs[kg, kg_zones := i.kg_zones, on = "continent"]

## Regression against KG zones ----
dist_within_tpgs[, round(range(dist_to_mean_tp),digits = 3), by = "continent"]
dist_within_tpgs[, .(mean_dist = mean(dist_to_mean_tp),
                     median_dist = median(dist_to_mean_tp)), by = "continent"]

# Simple regression
# null_lm <- lm(dist_to_mean_tp ~ 1,
#               data = dist_within_tpgs)
lm_tvar_env <- lm(dist_to_mean_tp ~ kg_zones,
                  data = dist_within_tpgs)
summary(lm_tvar_env) # R^2 super low
par(mfrow = c(2, 2))
plot(lm_tvar_env) # Okay?
estimate <- round(coef(lm_tvar_env)["kg_zones"], digits = 3)
p_val <-
  round(summary(lm_tvar_env)$coefficients[, "Pr(>|t|)"]["kg_zones"], digits = 3)

# GLM
# glm_tvar_env <- glm(dist_to_mean_tp ~ kg_zones,
#                     data = dist_within_tpgs,
#                     family = quasibinomial)
# summary(glm_tvar_env)
# null_glm <- glm(dist_to_mean_tp ~ 1,
#               data = dist_within_tpgs, 
#               family = quasibinomial)
# 
# # We can use the anova test and specify
# # the chisquare test for comparing the deviances
# # technically, this is a Likelihood ratio test
# anova(glm_tvar_env, null_glm, test = "Chisq")
# 
# # If we divide the Residual Deviance by the
# # degrees of freedom we yield a dispersion parameter
# # of approximately 1 which is very good
# # (remember that >2 or <0.5 would indicate over- or underdispersion)
# glm_tvar_env$deviance / glm_tvar_env$df.resid
# # Still underdispersion with quasibinomial model
# 
# # QQ plot
# qr <- qres.binom(glm_tvar_env)
# # we compute Dunn-Smyth-Residuals that are most appropriate for GLMs
# qqnorm(qr)

## Plot ----
# dist_within_tpgs[continent == "EU", ] %>% 
#   .[min(dist_to_mean_tp) == dist_to_mean_tp, ]
ggplot(dist_within_tpgs,
       aes(x = kg_zones, y = dist_to_mean_tp)) +
  geom_point(aes(color = continent),
             position = position_jitter(width = 0.1, seed = 1234)) +
  geom_smooth(method = "lm",#glm
              #method.args = list(family = "quasibinomial"),
              formula = y ~ x,
              se = FALSE,
              color = "black") +
  stat_summary(fun = "mean", 
               color = "grey",
               geom = "crossbar") +
  scale_color_d3(name = "Continent or \n Region",
                 labels = c("AUS", "EUR", "NA", "NZ", "SAf")) +
  scale_x_continuous(breaks = c(9, 14, 17, 18, 20))+ 
  labs(x = "Number of Köppen Geiger zones", 
       y = "Distance to mean trait profile of corresponding TPG") +
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
    y = 0.42,
    label = paste0("Estimate = ", estimate, "\np-value = ", p_val),
    color = "black",
    size = 4,
    fontface = "bold"
  ) #+ 
  # annotate(
  #   "curve",
  #   x = 18.05,
  #   y = 0.17,
  #   xend = 18,
  #   yend = 0.12,
  #   curvature = 0.2,
  #   arrow = arrow(length = unit(2, "mm"))
  # ) + 
  # annotate(
  #   "text",
  #   x = 18,
  #   y = 0.11,
  #   label = "Corduliidae"
  # )
ggsave(
  filename = file.path(
    data_paper,
    "Graphs",
    "Variation_within_tpgs_vs_kg.png"
  ),
  width = 35,
  height = 20,
  units = "cm"
)
