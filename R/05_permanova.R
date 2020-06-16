# ________________________________________________________
# PERMANOVA 
# ________________________________________________________

# from 04_FCA script
trait_comb <- rbindlist(file_list,
                        idcol = "region",
                        use.names = TRUE)


# clean region col
trait_comb[, region := sub("([A-z]{5})(\\_)([A-z]{2,})(\\_)([A-z]{3})(.+)", "\\3", region)]

# PERMANOVA data prep

# define region vec
region <- trait_comb$region

trait_comb[, family_p_region := paste0(family, "_",region)]
trait_comb[, c("family", "order", "region") := NULL]
setDF(trait_comb)
rownames(trait_comb) <- trait_comb$family_p_region

# del family_p_region col
trait_comb$family_p_region <- NULL

# distance matrix
# TODO: says chord distance in help, investigate  
trait_comb_dist <-  vegdist(decostand(trait_comb, "norm"), "euclidean")/sqrt(2)

# alternative: euclidean distance
# trait_comb_dist <- vegdist(trait_comb, method = "euclidean")


# PERMANVOA
# Taxa in trait databases differ stat. significantly 
# between regions
# 4.6 % of variance explained by region
pmv <- adonis(trait_comb_dist ~ region,
       permutations = 999#,
       #method = "euclidean"
       )
# test assumptions
bd <- betadisper(trait_comb_dist, region)
bd
# relatively similar distances

# However, distances are stat. significantly different (i.e. dispersion is different)
permutest(bd)

