# Loading data
genus_occ <- fread(file.path(data_in, "Genus_Occurrences.csv"))

# Identifying the different CRS
unique(genus_occ$Horizontal_datum)

# Seperating the different CRS
wgs <- genus_occ[genus_occ$Horizontal_datum == "WGS84",]
nad83 <- genus_occ[genus_occ$Horizontal_datum == "NAD83",]
nad27 <- genus_occ[genus_occ$Horizontal_datum == "NAD27",]

# Creating a spatial object using sf
wgs.spatial.sf <- st_as_sf(wgs, coords = c( x = "Longitude",
                                            y = "Latitude"),
                           crs = 4326 )

nad83.spatial.sf <- st_as_sf(nad83, coords = c( x = "Longitude",
                                                y = "Latitude"),
                             crs = 4269)

nad27.spatial.sf <- st_as_sf(nad27, coords = c( x = "Longitude",
                                                y = "Latitude"),
                             crs = 4267)

# Transforming into WGS84 using sf
nad83.to.wgs.sf <- st_transform(nad83.spatial.sf, crs= crs(wgs.spatial.sf))
nad27.to.wgs.sf <- st_transform(nad27.spatial.sf, crs= crs(wgs.spatial.sf))

# Unite spatial data frames
all.wgs.sf <- rbind(wgs.spatial.sf,nad27.to.wgs.sf,nad83.to.wgs.sf)
all.wgs.sf$Horizontal_datum <- "WGS84"

# Loading Koppen-Geiger-Data
KG_regions <- readGDAL(file.path(data_in, "Beck_KG_V1_present_0p083.tif"))
# plot(KG_regions)

# Checking KG data to have the same CRS as Conus data (WGS84)
st_crs(KG_regions) # Yes, it also uses WGS84 as CRS

# Turning KG_regions from a SpatialGridDataFrame to a Raster object 
# to utilize raster functions
KG_regions.raster <- raster(KG_regions)

# Extracting the Koppen Geiger Index data for each Conus Observation
all.wgs.sf$KoppenGeigerIndex <- raster::extract(KG_regions.raster, all.wgs.sf)

# Loading Koppen Geiger legend data 
KG_legend <- read.csv(file.path(data_in, "KG_legend.txt"))

# Adding Koppen Geiger descriptions to Conus observations
all.wgs.sf <- base::merge(all.wgs.sf,
                          KG_legend, 
                          by.x = "KoppenGeigerIndex",
                          by.y = "Index")

# Manually adding missing values
# Looked up coordinates "-117.835, 33.5814" in QGIS using the 
# Koppen Geiger Raster Map. 
# Results:  6,BSh,Arid_steppe_hot  
all.wgs.sf$Letter_code[c(1:3,5,6,8:10)] <- "BSh"
all.wgs.sf$KoppenGeigerIndex[c(1:3,5,6,8:10)] <- 6
all.wgs.sf$Full_description[c(1:3,5,6,8:10)] <- "Arid_steppe_hot"

# For coordinates: "-94.80336, 29.592028"
# Results: 14,Cfa,Temperate_nodryseason_hotsummer
all.wgs.sf$Letter_code[c(4,7)] <- "Cfa"
all.wgs.sf$KoppenGeigerIndex[c(4,7)] <- 14
all.wgs.sf$Full_description[c(4,7)] <- "Temperate_nodryseason_hotsummer"

# Exporting file
setDT(all.wgs.sf)
saveRDS(all.wgs.sf, file.path(data_in, "Conus_data_KoppenGeiger.rds"))
fwrite(all.wgs.sf, file.path(data_in, "Conus_data_KoppenGeiger.csv"))
