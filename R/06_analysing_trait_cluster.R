
#### trait cluster ####

# merge table with traits to cluster data
trait_cluster <- mapply(function(x, y)
  merge(x,
        y,
        by.x = "taxa",
        by.y = "family"),
  data_cluster,
  file_list,
  SIMPLIFY = FALSE)

# bind into one data.table
trait_cluster <- rbindlist(trait_cluster,
                           use.names = TRUE,
                           idcol = "region")

# remove order column and rename
trait_cluster[, order.y := NULL]  
setnames(trait_cluster, "order.x", "order") 

# nr. of clusters per region
trait_cluster[, max(groups), by = region]
trait_cluster[, .N, by = region]

# TODO: find clusters with similar combinations of traits
trait_cluster %>% 
  melt(., id.vars = c("region",
                      "taxa",
                      "order",
                      "groups")) %>% 
  .[value == 1, ] %>% 
  .[region %in% c("AUS", "NZ"), ] %>% 
  .[order(region, groups), ] %>% 
  .[, unique(variable), by = .(region, groups)] %>% 
  .[groups == 6, ]

trait_cluster %>% 
  melt(., id.vars = c("region",
                      "taxa",
                      "order",
                      "groups")) %>%
  .[value >= 0.5, ] %>% 
  .[, .N, by = .(region, groups, variable)] %>% 
  .[order(region, groups), ] %>% 
  .[N >= 10, ] %>% 
  ggplot(., aes(x = as.factor(variable), y = N))+
  geom_point()+
  facet_wrap(groups~region) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 75, hjust = 1))

# possible profile: volt_uni, resp_gil, ovip_aqu
trait_cluster %>% 
  melt(., id.vars = c("region",
                      "taxa",
                      "order",
                      "groups")) %>%
  .[value == 1, ] %>% 
  .[variable %in% c("volt_uni", "resp_gil", "ovip_aqu"), ] %>% 
  .[order(region, groups), ] %>% View

#
trait_cluster[locom_crawl == 1 & resp_gil == 1 & ovip_aqu == 1, ] %>% 
  .[order(region, groups), ]


# plot?
trait_cluster[order(region, groups), ] %>%
  .[, order.y := NULL] %>% 
  setnames(., "order.x", "order") %>% 
  melt(., id.vars = c("region",
                      "taxa",
                      "order",
                      "groups")) %>% 
  .[groups %in% c(1,2) & region == "AUS", ] %>% 
  ggplot(., aes(x = as.factor(variable), y = value))+
  geom_point(aes(color = groups))+
  facet_wrap(~region)
  


#### Most important traits ####
most_important_variables






