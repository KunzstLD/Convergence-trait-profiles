# __________________________________________________________________________________________________
#### Analysis of TPGs ####
# Identify traits that characterize the groups obtained through cluster analysis
# Clustering + Ordination  (principal coordinate analysis) together? 
# Consider correlations among traits (see Wilkes Paper!)
# __________________________________________________________________________________________________

# Read in trait data with group assignment
trait_CONT <- readRDS(file.path(data_cache, "trait_dat_grp_assig.rds"))

# Rm dev traits for now
trait_CONT <- trait_CONT[!trait %in% c("dev_hemimetabol", "dev_holometabol"),] 


# --- TPG's ----------------------------------------------------------------------------------------

# Unique families per continent and group
trait_CONT[, nr_families := uniqueN(family), by = .(continent, group)]

# Calculate proportion of taxa per continent, TPG and trait. Traits that are expressed at least with 
# an affinity of 0.5 per TPG should presumably distinguish the TPG's  
trait_CONT[affinity > 0.5, prop_taxa_express_05 := .N/nr_families, 
           by = .(continent,group, trait)]
disting_traits <-
  unique(trait_CONT[prop_taxa_express_05 >= 0.5, ] %>%
           .[order(-prop_taxa_express_05), .(continent, 
                                             group,
                                             trait,
                                             grouping_feature,
                                             prop_taxa_express_05,
                                             nr_families)], by = c("continent", "group", "trait"))
disting_traits <- disting_traits[order(continent, group, -prop_taxa_express_05), ] 

# Calculate nr of traits per tpg for upcoming comparisons
disting_traits[, nr_traits_group := .N, by = .(continent, group)]

# Are there similarities across continents?
# A) between the clustered groups
output_tpgs <- list()

for(cont in c("AUS", "EU", "NOA", "NZ")) {
  groups <- unique(disting_traits[continent == cont, group])
  similar_tpgs <- list()
  
  for (i in groups) {
    tpg <- disting_traits[continent == cont & group == i, trait]
    
    # length of subsetted tpg (equal to the nr_traits_group for this subset!)
    l_tpg <- length(tpg)
    
    similar_tpgs[[i]] <-
      disting_traits[trait %in% tpg,] %>%
      .[order(continent, group),] %>%
      .[, .(nr_comp = .N,
            trait,
            prop_taxa_express_05,
            nr_families,
            nr_traits_group),
        by = .(continent, group)] %>%
      .[nr_comp == l_tpg,]
  }
  output_tpgs[[cont]] <- similar_tpgs
}

# 1) select only tpgs where nr_traits_group equals nr_comp, i.e. the whole tpgs are identical 
# across regions, not only parts of it

# 1.1) Distinguish low and high affinities?
# 1.2) Pilliere found differences within orders (e.g. Ephemeroptera) -> something worth to check out


# B) In terms of grouping features contributing to the global distance
# What does the value exactly mean?
# From the docs of kdist.cor(): provides the correlations between the (squared) distances obtained
# for each trait and the global (squared) distances obtained by mixing all the traits
# (= contributions of traits to the global distances);
global_dist <- readRDS(file.path(data_cache, "global_dist.rds"))
