
# Mean trait profile per TPG
trait_CONT[continent == "AUS" & group %in% c(7,8), mean(affinity), by = c("group", "trait")] %>% 
  .[order(group, -V1), ]
tpgs_mean <-
  trait_CONT[, .(mean_profile = mean(affinity)), by = c("continent", "group", "trait")] 

tpgs_mean[continent == "EU" & group == 1, ]

# cor(tpgs_mean[continent == "AUS" & group == 1, mean_profile],
#     tpgs_mean[continent == "EU" & group == 1, mean_profile],
#     method = "spearman") # == 1

cor_results_AUS_EU <- matrix(nrow = 8, ncol = 9)
for(i_AUS in 1:8) {
  for (i_EU in 1:9) {
    cor_results_AUS_EU[i_AUS, i_EU] <-
      cor(tpgs_mean[continent == "AUS" & group == i_AUS, mean_profile],
          tpgs_mean[continent == "EU" & group == i_EU, mean_profile],
          method = "spearman")
  }
}

cor_results_AUS_NOA <- matrix(nrow = 8, ncol = 8)
for (i_AUS in 1:8) {
  for (i_NOA in 1:8) {
    cor_results_AUS_NOA[i_AUS, i_NOA] <-
      cor(tpgs_mean[continent == "AUS" &
                      group == i_AUS, mean_profile],
          tpgs_mean[continent == "NOA" &
                      group == i_NOA, mean_profile],
          method = "spearman")
  }
}

cor_results_AUS_NZ <- matrix(nrow = 8, ncol = 7)
for (i_AUS in 1:8) {
  for (i_NZ in 1:7) {
    cor_results_AUS_NZ[i_AUS, i_NZ] <-
      cor(tpgs_mean[continent == "AUS" &
                      group == i_AUS, mean_profile],
          tpgs_mean[continent == "NZ" &
                      group == i_NZ, mean_profile],
          method = "spearman")
  }
}

cor_results_AUS_SA <- matrix(nrow = 8, ncol = 10)
for (i_AUS in 1:8) {
  for (i_SA in 1:10) {
    cor_results_AUS_SA[i_AUS, i_SA] <-
      cor(tpgs_mean[continent == "AUS" &
                      group == i_AUS, mean_profile],
          tpgs_mean[continent == "SA" &
                      group == i_SA, mean_profile],
          method = "spearman")
  }
}

cor_results_ls <- list("AUS_EU" = cor_results_AUS_EU,
                       "AUS_NOA" = cor_results_AUS_NOA,
                       "AUS_NZ" = cor_results_AUS_NZ, 
                       "AUS_SA" = cor_results_AUS_SA)
pl_list <- list()
for(i in names(cor_results_ls)) {
  cor_results <- as.data.table(cor_results_ls[[i]])
  setnames(
    x = cor_results,
    old = paste0("V", 1:ncol(cor_results)),
    new = paste0(sub("AUS\\_","", i), "_TPG", 1:ncol(cor_results))
  )
  cor_results[, TPG_AUS := paste0("AUS_TPG", 1:nrow(cor_results))]
  
  pl_list[[i]] <- melt(cor_results, id.vars = "TPG_AUS") %>%
    .[order(TPG_AUS), ] %>%
    ggplot(., aes(x = TPG_AUS,
                  y = variable,
                  fill = value)) +
    geom_tile() +
    scale_fill_gradient2(
      limits = c(-1, 1),
      breaks = seq(-1, 1, by = 0.25),
      low = "lightblue",
      mid = "white",
      high = "darkblue",
      name = "Spearman correlation \n of mean trait profiles"
    ) +
    labs(x = "", y = "") +
    theme_bw() +
    theme(
      axis.title = element_text(size = 22),
      axis.text.x = element_text(family = "Roboto Mono",
                                 size = 14),
      axis.text.y = element_text(family = "Roboto Mono",
                                 size = 14),
      legend.title = element_text(family = "Roboto Mono",
                                  size = 14),
      legend.text = element_text(family = "Roboto Mono",
                                 size = 12),
      strip.text = element_text(family = "Roboto Mono",
                                size = 16),
      panel.grid = element_blank()
    )
}
lapply(names(pl_list), function(y)
  ggsave(
    filename = paste0(y, ".png"),
    plot = pl_list[[y]],
    path = "/home/kunzst/Schreibtisch"
  )
)
pl_list$AUS_EU
pl_list$AUS_NOA
pl_list$AUS_NZ
pl_list$AUS_SA


# From Martin Wilkes
get.cor <- function(x){
  samplei <- x[c(sample(1:nrow(tr), nrow(tr), replace=F)),]
  samplei <- as.data.frame(scale(samplei, center=T, scale=F))
  cor(samplei, method="spearman")
} #Function samples ntaxa (number of taxa in trait database) from resampled trait distributions and gives pairwise trait correlation matrix

cors.samples <- list()
for(i in 1:n.cors){
  cors.samples[[i]] <- get.cor(resampled.traits)
} #Null correlation matrix for comparison with observed correlation matrix
cors.samples



tpgs_mean_wf <- dcast(tpgs_mean, continent + group ~ trait, value.var = "mean_profile")
tpgs_mean_wf[, id := paste0(continent, "_", group)]
continent_names <- factor(tpgs_mean_wf$continent)
tpgs_mean_wf_cp <- copy(tpgs_mean_wf)
tpgs_mean_wf_cp[, c("group", "continent") := NULL]
setDF(tpgs_mean_wf_cp)

# Vegan way
# TODO: Do this with ade4 respecting the fuzzy coding system! (i.e. also normalize again)
# Problem that TPGs are also compared with themselves in distance matrix
# How do we get variance explained by continent
tpgs_dist <- vegdist(scale(tpgs_mean_wf_cp), method = "euclidean")

# PCA
tpgs_pca <- rda(tpgs_mean_wf_cp)
summary(tpgs_pca)
plot(tpgs_pca)

# RDA
tpgs_rda_continents <-
  rda(tpgs_mean_wf_cp ~ as.factor(continent),
      data = tpgs_mean_wf)
summary(tpgs_rda_continents)

# PCoA
tpgs_pcoa <- dbrda(tpgs_dist~1)
summary(tpgs_pcoa)

# Distance-based RDA
tpgs_dbrda <- dbrda(tpgs_mean_wf_cp ~ as.factor(continent), tpgs_mean_wf, distance = "euclidean")
summary(tpgs_dbrda)

 ## Basic Analysis
vare.cap <- capscale(varespec ~ N + P + K + Condition(Al), varechem,
                     dist="bray")
summary(vare.cap)
plot(vare.cap)
anova(vare.cap)

dbRDA=capscale(species001 ~ MAT+MWMT+MCMT+TD+lnMAP+lnMSP+lnAHM+lnSHM, environment,
               dist="bray")
plot(dbRDA)
anova(dbRDA)








