# ________________________________________________________
# Fuzzy correspondence analysis 
# ________________________________________________________

# bind trait data
trait_comb <- rbindlist(file_list,
                        idcol = "region",
                        use.names = TRUE)


# clean region col
trait_comb[, region := sub("([A-z]{5})(\\_)([A-z]{2,})(\\_)([A-z]{3})(.+)", "\\3", region)]

# create dataset with taxonomic information
taxa_comb <- trait_comb[, .(family, order, region)]

# create column with family name
# and region
# delete taxonomic information
# convert to df
# create rownames
trait_comb[, family_p_region := paste0(family, "_",region)]
trait_comb[, c("family", "order", "region") := NULL]
setDF(trait_comb)
rownames(trait_comb) <- trait_comb$family_p_region

# del family_p_region col
trait_comb$family_p_region <- NULL


####  fuzzy correspondence analysis ####
vec <- sub("\\_.*", "", names(trait_comb))
blocks <- rle(vec)$lengths

trait_comb <- prep.fuzzy.var(trait_comb, blocks)

# TODO: throws error that I quite don't understand
# chek out what's going wrong!
trait_fca <- dudi.fca(trait_comb[1:277, ])
trait_comb[277, ]

# PERMANOVA