#############

# Filter for families from NOA that cannot be used
# Look if these are related to the used taxa -> e.g. by shared order (?)
trait_non_agg[order %in% aq_insects, uniqueN(family), by = dataset]

trait_non_agg[order %in% aq_insects &
                dataset == "Traits_US_LauraT_pp_harmonized.rds",] 
trait_NOA <- trait_non_agg[order %in% aq_insects & dataset == "Traits_US_LauraT_pp_harmonized.rds",]
trait_SA <- trait_non_agg[order %in% aq_insects & dataset == "Trait_SA_pp_harmonised.rds",]
trait_NZ <- trait_non_agg[order %in% aq_insects & dataset == "Trait_NZ_pp_harmonized.rds",]
trait_CONT <- readRDS(file.path(data_cache, "trait_dat_grp_assig.rds"))


# Families not used in NOA
trait_NOA[, uniqueN(family), by = order]

trait_NOA[!family %in% trait_data_bind[continent == "NOA", family], ] %>% 
  .[, .(order, .N), by = family] %>% View()

# Not used families -> on family-level
trait_NOA[!family %in% trait_data_bind[continent == "NOA", family], ] %>% 
  .[is.na(species) & is.na(genus), family] %>% unique()

fam_noa <- trait_NOA[!family %in% trait_data_bind[continent == "NOA", family], .(order, .N), by = family] %>% 
  unique(.) %>% 
  .[order(-N, order), family] #%>% 
#  .[N >= 20, family]

# Families not used in SA
trait_SA[, uniqueN(family), by = order]

trait_SA[!family %in% trait_data_bind[continent == "SA", family], ] %>% 
  .[, .(order, .N), by = family]

fam_SA <- trait_SA[!family %in% trait_data_bind[continent == "SA", family], .(order, .N), by = family] %>% 
  unique(.) %>% 
  .[order(-N, order), family]

# Families with most taxa missing in NOA
# Chironomidae, Elmidae, Hydrophilidae, Tipulidae, Corduliidae, Caenidae, Gerridae, Ameletidae 
trait_NZ[family %in% c("Chironomidae", 
                       "Elmidae", 
                       "Hydrophilidae",
                       "Tipulidae",
                       "Corduliidae",
                       "Caenidae",
                       "Gerridae",
                       "Ameletidae"), family] %>% unique(.)

# From all 89 families that could not be used from NOA, 23 are displayed in the NZ dataset
trait_NZ[family %in% fam_noa, family] %>% unique(.)

# NOA: How is it in the global dataset? -> 42 covered through other continents
trait_CONT[family %in% fam_noa, family] %>% unique(.)
trait_CONT[family %in% fam_noa, .(family, order, continent, group, trait, affinity)] %>%
  dcast(., ...~trait, value.var = "affinity") %>% 
  .[, .(family, .N), by = .(order, continent)]

# ...without Australia: 39
trait_CONT[continent != "AUS" & family %in% fam_noa, family] %>% unique(.)

# SA: 61 could not be used, thereof 43 covered by other continents
trait_CONT[family %in% fam_SA, family] %>% unique(.)
# .. without Australia: 40
trait_CONT[continent != "AUS" & family %in% fam_SA, family] %>% unique(.)

#############