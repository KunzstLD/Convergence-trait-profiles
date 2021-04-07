# __________________________________________________________________________________________________
# Distribution of each trait for all families in each region
# __________________________________________________________________________________________________

# Data 
traits_ww <- readRDS(file.path(
  data_cache,
  "preproc_traits.rds")
)

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
  overview_lf$traits %>% unique
  
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
