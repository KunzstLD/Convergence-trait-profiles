
library(sp)
library(rgdal)

# Loading data
genus_occ <- fread(file.path(data_in, "NOA", "Genus_Occurrences.csv"))

# Identifying the different CRS
unique(genus_occ$Horizontal_datum)

# Seperating the different CRS
wgs <- genus_occ[genus_occ$Horizontal_datum == "WGS84", ]
nad83 <- genus_occ[genus_occ$Horizontal_datum == "NAD83", ]
nad27 <- genus_occ[genus_occ$Horizontal_datum == "NAD27", ]


# Transforming data into a spatial object
coords.wgs <- SpatialPoints(wgs[, c("Longitude", "Latitude")])
genus.wgs.spatial <- SpatialPointsDataFrame(coords.wgs, wgs)

proj4string(genus.wgs.spatial) <- CRS("+proj=longlat +elips=WGS84")

coords.nad83 <- SpatialPoints(nad83[, c("Longitude", "Latitude")])
genus.nad83.spatial <- SpatialPointsDataFrame(coords.nad83, nad83)
proj4string(genus.nad83.spatial) <- CRS("+proj=longlat +elips=NAD83")

coords.nad27 <- SpatialPoints(nad27[, c("Longitude", "Latitude")])
genus.nad27.spatial <- SpatialPointsDataFrame(coords.nad27, nad27)
proj4string(genus.nad27.spatial) <- CRS("+proj=longlat +elips=NAD27")

# Reprojecting to WGS84
nad83.to.wgs <- spTransform(genus.nad83.spatial, CRS("+proj=longlat +elips=WGS84"))
nad27.to.wgs <- spTransform(genus.nad27.spatial, CRS("+proj=longlat +ellps=WGS84"))
?spTransform


nad27.to.wgs@data[nad27.to.wgs$Unique_ID %in% c(3117, 3126, 3218, 3887), ]
genus_occ[genus_occ$Unique_ID %in% c(3117, 3126, 3218, 3887), ]


genus_occ[Horizontal_datum == "NAD83" & Latitude > 45, ]

nad83.to.wgs@data[nad83.to.wgs@data$Unique_ID %in% c(225), ]


ids <- nad27.to.wgs@data[, "Unique_ID"]
setDT(genus_occ)
all(genus_occ[Unique_ID %in% ids, "Longitude"] == nad27.to.wgs@data[, "Longitude"])


#### Use new script from Luic! ####
# read Conus_data_KoppenGeiger
conus_kg <- fread(file.path(data_in, "NOA", "Conus_data_KoppenGeiger.csv"))


# Q1) What is with the taxa that have no letter_code/full_description -> Ocean?
conus_kg[is.na(Letter_code), ]

conus_kg$Full_description %>% unique

# Plausability check
conus_kg[Full_description == "Arid_desert_cold", ]

# Create major climate regions
conus_kg[!is.na(Full_description), major_climate_regions := sub("(\\w)\\_.*","\\1",Full_description)]
conus_kg[!duplicated(Genus), ]