# Run Set_up.R script to load prerequisites for data preparation and analysis
source("./Comparison_between_climatic_zones/R/Set_up.R")

# 1) Checks: 
source("./Comparison_between_climatic_zones/R/02_checks.R")
# -> taxa on the same taxonomic level?
# -> Which orders are compared?

# 2) Data preparation
source("./Comparison_between_climatic_zones/R/01_data_prep_analysis.R")