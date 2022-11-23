# __________________________________________________________________________________________________
# Analysis null models ----
# TODO: - keep the iteration?
#       - non symmetric overlap
# Aim would be to reject the hypothesis 
# that the observed differences in niche space (perhaps accounting for central position and distribution) 
# are pure chance 
# - i.e. would arise from any distribution of taxa across continents. 
# The null model should therefore shuffle families across continents/regions, 
# taking the numbers into account.
# __________________________________________________________________________________________________
ovl_results <- readRDS(file.path(data_cache, "ovl_results.rds"))

## Convex hulls ----
### Non-symmetric overlap ----
ovl_results_non_symm <- lapply(ovl_results, function(y)
  as.data.table(colMeans(y$Overlap), keep.rownames = TRUE)) %>%
  rbindlist(.)
setnames(ovl_results_non_symm,
         c("V1", "V2"),
         c("continent", "overlap"))
ovl_results_non_symm[, type := "simulated"]

# load results from real dataset
real_ovl_non_symm <- readRDS(file.path(data_cache, "real_overlap_non_symm_pcoa.rds"))
real_ovl_non_symm[, type := "real"]
ovl_results_non_symm <- rbind(ovl_results_non_symm,
                              real_ovl_non_symm)


# Calculate p-value:
# Calculate number of test statistics as or more extreme than our initial ("real") test statistics 
# and divide by the total number of test-statistics calculated.
# In this case overlap values that are similar or smaller than the "real" overlap  
p_val_overlap <- list()
for(i in unique(ovl_results_non_symm$continent)) {
  p_val_overlap[[i]] <- ovl_results_non_symm[continent == i &
                                               type == "simulated" &
                                               overlap <= ovl_results_non_symm[continent == i &
                                                                                 type == "real", overlap], .N] /
    ovl_results_non_symm[continent == i &
                           type == "simulated", .N]
}
p_val_overlap <- cbind(p_val_overlap) %>% 
  as.data.table(., keep.rownames = TRUE)
setnames(p_val_overlap, "rn", "continent")
p_val_overlap[, `:=`(y = c(0.725, 0.825, 0.76, 0.805, 0.87),
                     x = rep(15, 5))]

# plot
wrap_names <- c(
  "AUS" = "AUS",
  "EU" = "EUR",
  "NOA" = "NA",
  "NZ" = "NZ",
  "SA" = "SA"
)
ggplot(ovl_results_non_symm, aes(y = overlap)) +
  geom_density(data = ~ .x[type == "simulated",]) +
  geom_hline(
    data = ~ .x[type == "real",],
    aes(yintercept = overlap),
    linetype = "dashed",
    color = "steelblue",
    show.legend = TRUE
  ) +
  geom_text(data = p_val_overlap,
            mapping = aes(x = x, y = y, label = paste("p-value: ",p_val_overlap))) +
  labs(x = "Density",
       y = "Mean overlap between convex hulls") +
  coord_flip() +
  facet_wrap( ~ as.factor(continent), 
              labeller = as_labeller(wrap_names)) +
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
    "Null_models_mean_ovl_cnxh.png"
  ),
  width = 35,
  height = 20,
  units = "cm"
)

# Functional symmetric overlap
# ovl_results <- readRDS(file.path(data_cache, "ovl_results.rds"))
# 
# ovl_results_symmetric <- lapply(ovl_results, function(y)
#   as.data.table(colMeans(y$Overlap_symmetric), keep.rownames = TRUE)) %>%
#   rbindlist(.)
# setnames(ovl_results_symmetric,
#          c("V1", "V2"),
#          c("continent", "symmetric_overlap"))
# ovl_results_symmetric[, type := "simulated"]
# 
# # load results from real dataset
# real_ovl_symm <-
#   readRDS(file.path(data_cache, "real_overlap_symm_pcoa.rds"))
# real_ovl_symm[, type := "real"]
# ovl_results_symmetric <- rbind(ovl_results_symmetric, real_ovl_symm)
# 
# # Distribution with symmetric overlap
# ggplot(ovl_results_symmetric, aes(y = symmetric_overlap)) +
#   geom_density(data = ~ .x[type == "simulated",]) +
#   geom_hline(
#     data = ~ .x[type == "real",],
#     aes(yintercept = symmetric_overlap),
#     linetype = "dashed",
#     color = "steelblue",
#     show.legend = TRUE
#   ) +
#   scale_linetype_manual(name="true value") +
#   labs(x = "Density", 
#        y = "Mean symmetric overlap between convex hulls") +
#   coord_flip() +
#   facet_wrap( ~ as.factor(continent)) +
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
#                               size = 14),
#     panel.grid = element_blank()
#   )

## Ellipses ----
ellipses_results <- readRDS(file = file.path(data_cache, "ellipses_results.rds"))
ellipses_results <- lapply(ellipses_results, function(y)
  as.data.frame(y)) %>%
  lapply(., function(z) as.data.table(z, keep.rownames = TRUE)) %>% 
  rbindlist(., idcol = "iteration")
setnames(ellipses_results, 
         "rn",
         "measure")
setnames(ellipses_results,
         names(ellipses_results),
         sub("X1\\.", "", names(ellipses_results)))

# Calculate % overlap for each comparison
# Mean for each continent
ellipses_results <-
  melt(
    ellipses_results,
    id.vars = c("iteration", "measure"),
    variable.name = "comparison"
  ) %>%
  dcast(., iteration + comparison ~ measure, value.var = "value")
ellipses_results[, prop_overlap := overlap/area.1]
ellipses_results[, first_continent := sub("([A-Z])(\\_)(.+)", "\\1", comparison)]
ellipses_results[, type := "simulated"]
ellipses_results_mean <-
  ellipses_results[, .(mean_overlap = mean(prop_overlap),
                       type),
                   by = c("iteration",
                          "first_continent")] 

# load results from real dataset
real_ovl_ellipses <- readRDS(file.path(data_cache, "real_ovl_ellipses.rds"))
real_ovl_ellipses[, comparison := sub("^1\\.", "", comparison)]
real_ovl_ellipses[, first_continent := sub("([A-Z])(\\_)(.+)", "\\1", comparison)]
real_ovl_ellipses[, type := "real"]
real_ovl_ellipses_mean <-
  real_ovl_ellipses[, .(mean_overlap = mean(prop_overlap), type),
                    by = "first_continent"]

ellipses_results_mean <-
  rbind(ellipses_results_mean, real_ovl_ellipses_mean,
        fill = TRUE)
ellipses_results_mean <- unique(ellipses_results_mean)

# Calc p-value
p_val_ellipses <- list()
for (i in unique(ellipses_results_mean$first_continent)) {
  p_val_ellipses[[i]] <- ellipses_results_mean[first_continent == i &
                                                 type == "simulated" &
                                                 mean_overlap <= ellipses_results_mean[first_continent == i &
                                                                                         type == "real",
                                                                                       mean_overlap], .N] /
    ellipses_results_mean[first_continent == i &
                            type == "simulated", .N]
}
p_val_ellipses <- cbind(p_val_ellipses) %>% 
  as.data.table(., keep.rownames = TRUE)
setnames(p_val_ellipses, "rn", "first_continent")
p_val_ellipses[, `:=`(y = c(0.8, 0.77, 0.78, 0.72, 0.895),
                     x = rep(8.5, 5))]

# plot
ggplot(ellipses_results_mean, aes(y = mean_overlap)) +
  geom_density(data = ~ .x[type == "simulated",]) +  
  geom_hline(
    data = ~ .x[type == "real",],
    aes(yintercept = mean_overlap),
    linetype = "dashed",
    color = "steelblue",
    show.legend = TRUE
  ) +
  geom_text(data = p_val_ellipses,
            mapping = aes(x = x, y = y, 
                          label = paste("p-value: ", p_val_ellipses))) +
  labs(x = "Density",
       y = "Mean overlap between ellipses") +
  coord_flip() +
  facet_wrap( ~ first_continent, 
              labeller = as_labeller(wrap_names)) +
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
    "Null_models_mean_ovl_ellipses.png"
  ),
  width = 35,
  height = 20,
  units = "cm"
)

# Comparison explained variances of dbrda ----
expl_var <- readRDS(file.path(data_cache, "expl_var.rds"))
# P value
p_val_explvar <- expl_var[type == "simulated" &
           var_constr >= expl_var[type == "real", var_constr], .N]/
  expl_var[type == "simulated", .N]
p_val_explvar <- as.data.table(p_val_explvar)
p_val_explvar[, `:=`(y = 0.16,
                     x = 40)]

# Plot
ggplot(expl_var, aes(y = var_constr)) +
  geom_density(data = ~ .x[type == "simulated", ]) +
  geom_hline(
    data = ~ .x[type == "real", ],
    aes(yintercept = var_constr),
    linetype = "dashed",
    color = "steelblue",
    show.legend = TRUE
  ) +
  geom_text(p_val_explvar, mapping = aes(
    y = y,
    x = x,
    label = paste("p-value: ", p_val_explvar)
  )) +
  labs(x = "Density",
       y = "Explained variances of continents and regions") +
  coord_flip() +
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
    "Null_models_expl_variance.png"
  ),
  width = 35,
  height = 20,
  units = "cm"
)

# RF results ----
# Load five most important traits for each null trait dataset
imp_traits_sim <- readRDS(file.path(data_cache, "sim_five_most_imp_traits.rds"))

# most occurring most important traits 
# In every dataset the same?
imp_traits_sim_n <- imp_traits_sim[, .N, by = c("continent", "trait")]

# bf spherical was never among the five most important traits in AUS, EU, NA and SA
# imp_traits_sim_n[trait == "bf_spherical", ]
imp_traits_sim_n <- rbind(imp_traits_sim_n, data.table(
  continent = c("AUS", "EU", "NOA", "SA"),
  trait = rep("bf_spherical", 4),
  N = c(0, 0, 0, 0)
))
imp_traits_sim_n[, grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", trait)]
five_most_abund_traits <- imp_traits_sim_n[order(continent, -N), .SD[1:5], by = "continent"]
imp_traits_sim_n[five_most_abund_traits,
                 five_most_abund := i.trait,
                 on = c("continent", "trait")]
# map colors for plotting
imp_traits_sim_n[, color := fcase(!is.na(five_most_abund),
                                  "mediumpurple1",
                                  is.na(five_most_abund),
                                  "gray70")]

# plot
wrap_names <- c(
  "AUS" = "AUS",
  "EU" = "EUR",
  "NOA" = "NA",
  "NZ" = "NZ",
  "SA" = "SA",
  "resp" = "Resp.",
  "size" = "Size",
  "feed" = "Feed. m.",
  "locom" = "Locom.",
  "bf" = "Body f.",
  "volt" = "Volt."
)

ggplot(imp_traits_sim_n,
       aes(x = trait,
           y = N,
           fill = color)) +
  geom_col() +
  geom_text(mapping = aes(
    x = as.factor(trait),
    y = N + 15,
    label = trait,
    hjust = 0.05
  ),
  size = 4.1) +
  facet_grid(grouping_feature ~ continent,
             labeller = as_labeller(wrap_names),
             scales = "free") +
  coord_flip() +
  scale_fill_identity(guide = "none") +
  lims(y = c(0, 900)) +
  labs(x = "",
       y = "Number of times a trait was among the most five important traits for TPG selection") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(family = "Roboto Mono",
                               size = 14),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(family = "Roboto Mono",
                                size = 16),
    legend.text = element_text(family = "Roboto Mono",
                               size = 14),
    strip.text = element_text(family = "Roboto Mono",
                              size = 14),
    legend.position = "none",
    panel.grid = element_blank()
  )
ggsave(
  filename = file.path(data_paper,
                       "Graphs",
                       "Null_models_rf_most_imp_traits.png"),
  width = 50,
  height = 30,
  units = "cm"
)