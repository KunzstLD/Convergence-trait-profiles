# __________________________________________________________________________________________________
# PcOA
# TODO: check scripts for function to find suitable nr. of axis
# https://github.com/tanogc/overarching_functional_space
# TODO: Between group analysis better suited than discriminant analysis? 
# What about within groups analysis?
# __________________________________________________________________________________________________

#### Test example
traits_ww <- readRDS(file.path(
  data_cache,
  "preproc_traits.rds")
)
test <- dist.ktab(traits_ww$AUS, type = "F") 

test_pco <- ade4::dudi.pco(test, scannf = FALSE, nf = 10)
print(test_pco)
scatter(test_pco)

# Explained variance by each axis
round(test_pco$eig[1:10]/sum(test_pco$eig),2)

# Cumulative variance explained by all FS axes 
cumsum(test_pco$eig)[10]/sum(test_pco$eig) 

# __________________________________________________________________________________________________
#### PCoA & Discriminant analysis with trait data from all 4 regions ####
# __________________________________________________________________________________________________
trait_data_bind <- readRDS(file.path(data_cache,
                                     "trait_data_ww_bind.rds"))

# preproc
trait_data_bind[, id := paste0(continent, "_", family)]
continent_names <- factor(trait_data_bind$continent)
trait_data_bind[, c("family", "order", "continent") := NULL]
setDF(trait_data_bind)

# Add row.names
taxa_reg_names <- trait_data_bind$id
row.names(trait_data_bind) <- taxa_ref_names 
trait_data_bind$id <- NULL

# data with only two levels/integer/factor variables have to be numeric
col_int_fac <-
  names(Filter(function(y) {
    is.integer(y) | is.factor(y)
  }, trait_data_bind))
trait_data_bind[, col_int_fac] <- apply(trait_data_bind[, col_int_fac], 2, as.double)

# Convert to ktab object
vec <- sub("\\_.*", "\\1", names(trait_data_bind))
blocks <- rle(vec)$lengths
trait_data_bind <- prep.fuzzy(trait_data_bind, blocks)
trait_data_bind <- ktab.list.df(list(trait_data_bind))
comb_dist <- dist.ktab(trait_data_bind, type = "F") 

# PcOA
comb_pcoa <- dudi.pco(comb_dist, scannf = FALSE)
scatter(comb_pcoa)

# Explained variance by each axis
round(comb_pcoa$eig[1:10]/sum(comb_pcoa$eig),2)

# Cumulative variance explained by first ten axes 
cumsum(comb_pcoa$eig)[10]/sum(comb_pcoa$eig) 

# Discriminant analysis: between class comparisons ----
# TODO: Check what to do when assumptions not met
# TODO: What do the group scores mean?
# Prior to discriminant analysis the assumption of 
# homogeneous within-group dispersion should be tested
# Differences are not too big, but assumption of homogeneous within-group dispersipon
# is actually not met!
bd_comb <- betadisper(comb_dist,
                      group = continent_names)
bd_comb
boxplot(bd_comb)
permutest(bd_comb)
 
comb_discrim <- discrimin(comb_pcoa, continent_names, scannf = FALSE)
plot(test_discrim, row.pellipses.col = adegpar()$ppalette$quali(4))
# [c(4,6)]
# Explanations: 
# - Row_scores and classes:
#     "Shows the projections of the sites onto the plane defined by the axes of the
#     Discriminant Analysis"
#     test_discrim$li
# - Class scores: 
#     Group scores (four points -> four regions)
#     test_discrim$gc
# - Unconstrained axis: 
#     "The following graph on the left (Unconstrained axes) shows the projection of the axes
#      of the initial Principal Component Analysis ($cp data frame) to understand the
#     relationships between the initial PCA and the Discriminant Analysis"
# - Eigenvalues: Usual eigenvalues barchart
# - Loadings: 
#     "Coefficients of the variables"
#     test_discrim$fa
# - Columns: 
#     "Covariance between traits and axes of the analysis"
#     test_discrim$va

# Test for stat. significance of differences between seasons
# ~ 4.5 % of the explained inertia comes from the difference between continents
(rt_discrim <- randtest(test_discrim))



