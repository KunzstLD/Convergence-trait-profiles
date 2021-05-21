# ______________________________________________________________________________
#### Hierarchical clustering & Visualization ####
# TODO: test the influence of ward.D, ward.D2 and
# other, see also Everitt and
# https://cran.r-project.org/web/packages/dendextend/vignettes/Cluster_Analysis.html#animals---attributes-of-animals
# TODO: Change to genus or family, use rather all data/aggregated data?
# TODO: Implement chord distance for fuzzy codes
# ______________________________________________________________________________

# Load hc_output
hc_output <- readRDS(file.path(
  data_cache,
  "hc_output_temperate.rds"
))

# Read orig. data for taxonomic and trait information
trait_subsets <- readRDS(file.path(data_cache, "trait_subsets.rds"))
trait_agg <- readRDS(file.path(data_cache, "trait_agg_genus.rds"))

# ---- Dendrogram & distance matrices -----------------------------------------

# extract distance matrices
ls_dist_mat <- lapply(hc_output, function(y) y[["distance_matrix"]])

# construct dendrogram and extract labels
dend_output <- list()
for (i in names(ls_dist_mat)) {
  dist_mat <- ls_dist_mat[[i]]

  # dendrogram
  hc_taxa <- hclust(as.dist(dist_mat), method = "ward.D2")

  # get labels of dendrogram
  # merge the according order information to the labels
  label <- hc_taxa %>%
    as.dendrogram() %>%
    labels()

  # save output
  dend_output[[i]] <- list(
    "dend_label" = label,
    "hc_element" = hc_taxa
  )
}

# Dendrogram with species as labels
# TODO: Needs adjustment for species names
png(
  file = file.path(
    data_grp,
    paste0("Dendrogram_", sub("\\.rds", "", "temperate_EU"), ".png")
  ),
  width = 1100,
  height = 1300,
  res = 100
)
dend_output$temperate$hc_element %>%
  as.dendrogram() %>%
  color_branches(k = hc_output$temperate$optimal_nog) %>%
  hang.dendrogram(hang_height = -0.5) %>%
  set("labels_cex", 0.75) %>%
  dendextend::ladderize() %>%
  plot(horiz = TRUE)
dev.off()

#---- Save clustered groups for further analysis -------------------------

# get groups
dend_taxa <- as.dendrogram(dend_output$temperate$hc_element)

# grouping of taxa
# PTODO: could merge cutree grouping
data_cluster <-
  data.table(
    genus = names(cutree(dend_taxa,
      k = hc_output$temperate$optimal_nog
    )),
    group = cutree(dend_taxa, k = hc_output$temperate$optimal_nog)
  )

# merge back taxonomic and trait information
cols <- grep("genus|family|order|feed.*|resp.*|volt.*|locom.*|ovip.*",
  names(trait_subsets),
  value = TRUE
)

# add taxonomic information and traits
data_cluster[trait_subsets,
  `:=`(
    family = i.family,
    order = i.order
  ),
  on = "genus"
]

saveRDS(
  object = data_cluster,
  file = file.path(
    data_cache,
    "data_cluster_temperate.rds"
  )
)

# ________________________________________________________________
#### Display traits as dendrogram & trait profiles as heatmap ####
# ________________________________________________________________

# clustering of traits -> transpose trait data
# trait_transpose <- as.data.frame(t(data[, -grep("order|family", names(data))]))
# # change columns to ordered
# cols <- names(trait_transpose)
# trait_transpose[, cols] <- lapply(trait_transpose[, cols], as.ordered)
#
# # transform "binary" columns to numeric
# col_two_levels <- trait_transpose %>%
#   lapply(levels) %>%
#   lapply(., function(y)
#     length(y)[length(y) == 2]) %>%
#   unlist %>%
#   names
# if(length(col_two_levels) == 1) {
#   trait_transpose[[col_two_levels]] <- sapply(trait_transpose[[col_two_levels]], as.numeric)
# }
# if (length(col_two_levels) > 1) {
#   trait_transpose[, col_two_levels] <- lapply(trait_transpose[, col_two_levels], as.numeric)
# }
#
# # dist mat
# gower_dist_trait <- gowdis(trait_transpose,
#                            ord = "podani")
# # hc object
# hc_trait <- hclust(gower_dist_trait, method = "ward.D2")
#
# # plot & safe
# png(
#   file = file.path(data_out, "Graphs", paste0(
#     "Dendrogram_traits_podani_", name_dataset,".png"
#   )),
#   width = 1100,
#   height = 1300,
#   res = 100
# )
# hc_trait %>% as.dendrogram() %>%
#   #  color_branches(k = optimal_nog) %>%
#   hang.dendrogram(hang_height = 0.01) %>%
#   set("labels_cex", 0.75) %>%
#   dendextend::ladderize() %>%
#   plot(horiz = TRUE, main = paste("Dendrogram traits", name_dataset))
# dev.off()

#### heatmap
# does not work properly yet
# some_col_func <- function(n) {
#   colorspace::diverge_hcl(n,
#     h = c(246, 40),
#     c = 96, l = c(65, 90)
#   )
# }
# # par(mar = c(3,3,3,3))
# # library(gplots)
# cols <- names(data[,-1])
# data[, cols] <- lapply(data[,cols], as.numeric)
#
# gplots::heatmap.2(as.matrix(data[,-1]),
#                   main = "Macroinvertebrate traits",
#                   srtCol = NULL,
#                   Rowv = dend_taxa_metric,
#                   Colv = dend_trait,
#                   trace="row", hline = NA, tracecol = "darkgrey",
#                   margins =c(6,3),
#                   key.xlab = "no / yes",
#                   denscol = "grey",
#                   density.info = "density",
#                   col = some_col_func
# )