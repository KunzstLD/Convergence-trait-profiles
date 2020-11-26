# ________________________________________________________
# Analysis HC with approrpiate distance 
# for fuzzy coded traits 

# TODO: Analyse per Order?
# TODO: Try using chord distance with blocks
# (was it really the chord distance I was trying to use?)
# TODO: Optimal number of groups to test for?
# ________________________________________________________

# load data
preproc_data <- readRDS(file.path(data_cache,
                                  "preproc_data.rds"))

preproc_data_omit <- readRDS(file.path(data_cache,
                                       "preproc_data_omit.rds"))

# for interactive use
if (interactive()) {
  omitted <-
    readline("Use data where NAs are OMITTED? \n Otherwise the FULL data with NAs substituted with 0
             are used")
  
  if (omitted == "T" | omitted == "TRUE" | omitted == "True") {
    data_type <- "omit"
  } else{
    data_type <- "full"
  }
  rm(omitted)
} 
# else choose data_type manually
# data_type = "omit"
data <- switch(data_type, "full" = preproc_data, "omit" = preproc_data_omit)


# lists to store values
hc_output <- list()
 
set.seed(1234)
for (i in names(data)) {

  # dist. matrix
  # binary variables treated as interval scaled
  dist_mat <- as.matrix(daisy(x = data[[i]]))

  # Optimal number of groups
  # Using the gap statistic
  gap <- clusGap(
    x = dist_mat,
    FUN = mycluster_hc,
    K.max = 15,
    B = 100
  )
  # plot(gap)
  # determines location of maximum

  # store optimal nr of groups, gap stat. and dist. matrix
  optimal_nog <- maxSE(gap$Tab[, "gap"],
                       gap$Tab[, "SE.sim"],
                       method = "Tibs2001SEmax")

  hc_output[[i]] <- list(
    "distance_matrix" = dist_mat,
    "gap_statistic" = gap,
    "optimal_nog" = optimal_nog
  )
}
saveRDS(object = hc_output,
        file = file.path(data_cache,
                         paste0("hc_output_", data_type, ".rds")))