# __________________________________________________________________________________________________
# Plotting of convex hulls and calculation of intersection
# Use pcoa scores from 02_1_PCOA.R
# __________________________________________________________________________________________________

# Load PCoA scores and PCoA object
pcoa_scores <- readRDS(file = file.path(data_cache, "pcoa_scores.rds"))
pcoa_scores[, family := sub("(.+)\\_", "", id)]
comb_pcoa <- readRDS(file.path(data_cache, "comb_pcoa.rds"))

# Calculate convex hull
hull <- pcoa_scores %>%
  group_by(continent) %>%
  slice(chull(A1, A2))
setDT(hull)
# saveRDS(hull, file.path(data_cache, "hull_pcoa.rds"))

# With convex hulls
# Base plot of PCOA scores (first two axes)
# needs comb_pcoa from 02_01_PCOA.R
pcoa_base_pl <- ggplot(pcoa_scores, aes(x = A1, y = A2)) +
  geom_point(alpha = 0.15) +
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

# Combine both plots?
pcoa_final_pl <- pcoa_base_pl + geom_polygon(
  data = hull,
  alpha = 0.05,
  aes(color = continent)
)
# pcoa_base_pl + geom_polygon(
#   data = hull,
#   alpha = 0.05,
#   aes(color = continent)
# ) +
#   facet_wrap(~continent)
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
pcoa_final_pl + geom_text_repel(
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
    ),]
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

# Calculate overlap ----
# library(geometry)
hull_split <- split(hull[, .(A1, A2)], f = hull$continent)
Overlap(hull_split)
Overlap(hull_split, symmetric = TRUE)

# Save for comparison with null models
real_overlap_symm <- colMeans(Overlap(hull_split, symmetric = TRUE)) %>% 
  as.data.table(., keep.rownames = TRUE)
setnames(real_overlap_symm,
         c("rn", "."), 
         c("continent", "symmetric_overlap"))
saveRDS(real_overlap_symm, file.path(data_cache, "real_overlap_symm_pcoa.rds"))

real_overlap_non_symm <- colMeans(Overlap(hull_split)) %>%
  as.data.table(., keep.rownames = TRUE)
setnames(real_overlap_non_symm,
         c("rn", "."),
         c("continent", "overlap"))
saveRDS(real_overlap_non_symm,
        file.path(data_cache, "real_overlap_non_symm_pcoa.rds"))
