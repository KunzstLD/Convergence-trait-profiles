# Ellipses with SIBER
# Start in 02_01_PCOA.R after base plot
# Ellipses seem to wide, especially for NZ
pcoa_siber <- createSiberObject(pcoa_scores[, .(iso1 = A1, iso2 = A2, group = continent, community = 1)])
 
# # Plot with ellipses
pcoa_base_pl + stat_ellipse(
  aes(
    color = continent
  ),
  alpha = 0.1,
  level = 0.95,
  type = "norm",
  geom = "polygon"
)
# 
# # Calculate group metrics
# # SEA - standard ellipse area
# # SEAc - sample size corrected standard ellipse area
group.ML <- groupMetricsML(pcoa_siber)
# 
# 
# # fit the ellipses which uses an Inverse Wishart prior
# # on the covariance matrix Sigma, and a vague normal prior on the 
# # means. Fitting is via the JAGS method.
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
# 
# # The posterior estimates of the ellipses for each group can be used to
# # Not sure if I understand these 100 %, New Zealand seems to have a few outliers 
# # multimodal data will be a problem
# # histogram(SEA.B[, 5])
SEA.B <- siberEllipses(ellipses.posterior)
siberDensityPlot(SEA.B, xticklabels = colnames(group.ML),
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)
# 
# # Overlap between ellipses
# # The first ellipse is referenced using a character string representation where 
# # in "x.y", "x" is the community, and "y" is the group within that community.
# # So in this example: community 1, group 2
ellipse1 <- "1.AUS" 
# 
# # Ellipse two is similarly defined: community 1, group3
ellipse2 <- "1.EU"
# 
# # The overlap of the maximum likelihood fitted standard ellipses are 
# # estimated using
sea.overlap <- maxLikOverlap(ellipse1, ellipse2, pcoa_siber,
                             p.interval = NULL, n = 100)
 
# # the overlap betweeen the corresponding 95% prediction ellipses is given by:
ellipse95.overlap <- maxLikOverlap(ellipse1, ellipse2, pcoa_siber,
                                   p.interval = 0.95, n = 100)

# # so in this case, the overlap as a proportion of the non-overlapping area of 
# # the two ellipses, would be
prop.95.over <- ellipse95.overlap[3] / (ellipse95.overlap[2] +
                                          ellipse95.overlap[1] -
                                          ellipse95.overlap[3])
