# __________________________________________________________________________________________________
# Distribution of each trait for all families in each region
# __________________________________________________________________________________________________

# Data 
traits_ww <- readRDS(file.path(
  data_cache,
  "preproc_traits.rds")
)

# __________________________________________________________________________________________________
# Trait distributions ----
# __________________________________________________________________________________________________

## For aggregated trait data ----
# across all aggreagted families per region and grouping feature
for(region in c("AUS", "EU", "NOA", "NZ")) {
  # preprocessing
  overview <- traits_ww[[region]]$Ana1
  overview$family <- rownames(overview)
  setDT(overview)
  
  overview_lf <- melt(overview,
                      id.vars = "family",
                      variable.name = "traits")
  overview_lf[, grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", traits)]
  
  overview_lf[, trait_label := as.character(sub("([a-z]{1,})(\\_)(.+)", "\\3", traits))]
  
  # plotting
  facet_names <- c(
    "bf" = "Body form",
    "feed" = "Feeding mode",
    "locom" = "Locomotion",
    "ovip" = "Oviposition",
    "resp" = "Respiration",
    "size" = "Size",
    "volt" = "Voltinism"
  )
  
  pl <- ggplot(overview_lf, aes(x = factor(trait_label),
                                y = value)) +
    geom_jitter(width = 0.05) +
    geom_violin(alpha = 0.2) +
    facet_wrap(grouping_feature ~ .,
               scales = "free",
               labeller = as_labeller(facet_names)) +
    labs(x = "Traits",
         y = "Affinity") +
    ggtitle(paste("Trait affinity distribution", region)) +
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
  ggsave(
    plot = pl,
    filename = file.path(
      data_paper,
      "Graphs",
      paste0("Trait_distribution_", region, ".png")
    ),
    width = 52,
    height = 30,
    units = "cm"
  )
}

## For not-aggregated harmonized trait data ----

### Australia ----
aus_preaggr <- readRDS(file = file.path(data_in, "Trait_AUS_harmonized_preaggre.rds"))

aus_preaggr <- melt(
  aus_preaggr,
  id.vars = c("species", "genus", "family", "order", "taxa", "unique_id"),
  variable.name = "traits"
) 
aus_preaggr[, grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", traits)]

# subset to relevant families
aus_preaggr <- aus_preaggr[family %in% rownames(traits_ww$AUS$Ana1),]
aus_preaggr <- aus_preaggr[!is.na(value), ] 

# plot trait distribution
aus_preaggr[family %in% c("Chironomidae", 
                          "Leptophlebiidae", 
                          "Conoesucidae", 
                          "Gripopterygidae",
                          "Baetidae",
                          "Leptoceridae"),] %>%
  .[grouping_feature != "dev", ] %>%
  ggplot(., aes(x = as.factor(traits), y = value)) +
  geom_jitter(width = 0.2) +
  geom_violin() +
  stat_summary(fun = "median", col = "red") +
  stat_summary(fun = "mean", col = "forestgreen") +
  facet_grid(family ~ grouping_feature ,
             scales = "free") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(
      family = "Roboto Mono",
      size = 14,
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
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

# Calculate sd per family
aus_preaggr[, s_dev := sd(value), by = .(family, traits)]
aus_preaggr[family %in% c("Chironomidae", 
                          "Leptophlebiidae", 
                          "Conoesucidae", 
                          "Gripopterygidae",
                          "Baetidae",
                          "Leptoceridae"),] %>%
  ggplot(., aes(x = as.factor(family), y = s_dev)) +
  geom_point() +
  facet_wrap(. ~ traits) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(
      family = "Roboto Mono",
      size = 14,
      angle = 90,
      hjust = 1,
      vjust = 0.3
    ),
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
unique(aus_preaggr[, .(family, order, traits, grouping_feature, value, s_dev)]) %>% 
  .[, sum(s_dev, na.rm = TRUE), by = family] %>% 
  .[V1 >= 10, ]

# Idea: Compare distributions for the same traits based on cdf!
# https://stats.stackexchange.com/questions/25764/clustering-distributions
test <- aus_preaggr[family == "Leptophlebiidae" & traits == "volt_uni", value]
P = ecdf(test) 
plot(P)

# __________________________________________________________________________________________________
# Mean-var relationship ----
# __________________________________________________________________________________________________
trait_data_bind <- readRDS(file.path(
  data_cache,
  "trait_data_ww_bind.rds"
))

# preproc
trait_data_bind[, id := paste0(continent, "_", family)]
df <- copy(trait_data_bind)

df_lf <- melt(
  df,
  id.vars = c("family",
              "order",
              "continent",
              "id"),
  variable.name = "trait",
  value.name = "affinity"
)

# Plot mean-var relationship 
df_lf[, .(mean_vals = mean(affinity),
          var_vls = var(affinity)),
      by = .(continent, trait)] %>%
  ggplot(., aes(x = mean_vals,
                y = var_vls,
                color = continent)) +
  facet_grid( ~ as.factor(continent)) +
  geom_point() +
  geom_smooth(method = "glm") +
  theme_bw()










