# __________________________________________________________________________________________________
# SIBER Ellipses ----
# __________________________________________________________________________________________________
pcoa_scores <- readRDS(file = file.path(data_cache, "pcoa_scores.rds"))

# Create SIBER object
pcoa_siber <- createSiberObject(pcoa_scores[, .(
  iso1 = A1,
  iso2 = A2,
  group = continent,
  community = 1
)])

# Calculate group metrics
# - Convex hull total area
# - SEA - standard ellipse area
# - SEAc - sample size corrected standard ellipse area
group.ML <- groupMetricsML(pcoa_siber)

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3
ellipses.posterior <- siberMVN(pcoa_siber, parms, priors)

# The posterior estimates of the ellipses for each group can be used to
# multimodal data will be a problem
# histogram(SEA.B[, 5])
SEA.B <- siberEllipses(ellipses.posterior)
siberDensityPlot(SEA.B, xticklabels = colnames(group.ML),
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)

### Overlap between ellipses ----
# The first ellipse is referenced using a character string representation where 
# in "x.y", "x" is the community, and "y" is the group within that community.
# e.g, ellipse1 <- "1.AUS" , ellipse2 <- "1.EU"
# The overlap of the maximum likelihood fitted standard ellipses
# sea.overlap <- maxLikOverlap(ellipse1, ellipse2, pcoa_siber,
#                              p.interval = NULL, n = 100)
# # The overlap between the corresponding 95% prediction ellipses is given by:
# ellipse95.overlap <- maxLikOverlap(ellipse1, ellipse2, pcoa_siber,
#                                    p.interval = 0.95, n = 100)
real_ovl_ellipses <- calc_ellipses(scores = pcoa_scores)
real_ovl_ellipses <- lapply(real_ovl_ellipses, function(y)
  as.data.table(y, keep.rownames = TRUE)) %>%
  rbindlist(., idcol = "comparison")
setnames(real_ovl_ellipses, 
         "rn",
         "measure")
real_ovl_ellipses <- dcast(real_ovl_ellipses, comparison ~ measure, value.var = "y")
real_ovl_ellipses[, prop_overlap := overlap/area.1]
saveRDS(real_ovl_ellipses,
        file.path(data_cache, "real_ovl_ellipses.rds"))

# Overlap with posterior draws?
# Computationally intensive, would need to average the different draws?

## Plot ellipses ----
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

# Standard ellipses with normal distribution (how is this done exactly?)
# pcoa_base_pl + stat_ellipse(
#   aes(
#     color = continent
#   ),
#   alpha = 0.1,
#   level = 0.95,
#   type = "norm",
#   geom = "polygon"
# ) +
#   facet_wrap(~continent)
# ggsave(
#   filename = file.path(
#     data_paper,
#     "Graphs",
#     "PCOA_continent_ellipses_combined.png"
#   ),
#   width = 35,
#   height = 20,
#   units = "cm"
# )

# Posterior ellipses 
# Slighly modified from Andrew Jackson
# https://github.com/AndrewLJackson/SIBER/blob/master/vignettes/Plot-posterior-ellipses.Rmd
# how many of the posterior draws do you want?
n.posts <- 15
# decide how big an ellipse you want to draw
p.ell <- 0.95
# for a standard ellipse use
# p.ell <- pchisq(1,2)
# a list to store the results
all_ellipses <- list()
# loop over groups
for (i in 1:length(ellipses.posterior)) {
  # a dummy variable to build in the loop
  ell <- NULL
  post.id <- NULL
  
  for (j in 1:n.posts) {
    # covariance matrix
    Sigma  <- matrix(ellipses.posterior[[i]][j, 1:4], 2, 2)
    
    # mean
    mu     <- ellipses.posterior[[i]][j, 5:6]
    
    # ellipse points
    out <- ellipse::ellipse(Sigma, centre = mu , level = p.ell)
    
    
    ell <- rbind(ell, out)
    post.id <- c(post.id, rep(j, nrow(out)))
    
  }
  ell <- as.data.frame(ell)
  ell$rep <- post.id
  all_ellipses[[i]] <- ell
}
ellipse_df <- bind_rows(all_ellipses, .id = "id")

# now we need the group and community names
# extract them from the ellipses.posterior list
group_comm_names <- names(ellipses.posterior)[as.numeric(ellipse_df$id)]
# split them and conver to a matrix, NB byrow = T
split_group_comm <- matrix(unlist(strsplit(group_comm_names, "[.]")),
                           nrow(ellipse_df), 2, byrow = TRUE)
ellipse_df$community <- split_group_comm[,1]
ellipse_df$continent     <- split_group_comm[,2]
ellipse_df <- dplyr::rename(ellipse_df, iso1 = x, iso2 = y)

# rename columns of ellipse_df to match the aesthetics
pcoa_base_pl + facet_wrap(~as.factor(continent)) +
  geom_polygon(data = ellipse_df,
               mapping = aes(iso1, 
                             iso2,
                             group = rep,
                             color = as.factor(continent),
                             fill = NULL),
               fill = NA)
ggsave(
  filename = file.path(
    data_paper,
    "Graphs",
    "PCOA_continent_ellipses_posterior.png"
  ),
  width = 35,
  height = 20,
  units = "cm"
)