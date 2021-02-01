# Run Set_up.R script to load prerequisites for data preparation and analysis
source("./Comparison_between_climatic_zones/R/Set_up.R")

# ______________________________________________________________________________
# 1) Checks:
source("./Comparison_between_climatic_zones/R/01_checks.R")
# -> taxa on the same taxonomic level?
# -> Which orders are compared?
# ______________________________________________________________________________

# 2) Data preparation
source("./Comparison_between_climatic_zones/R/02_data_prep_analysis.R")

# 3) Hierarchical cluster analysis (HCA)
source("./Comparison_between_climatic_zones/R/03_HC.R")


# 4) Visualisation HCA