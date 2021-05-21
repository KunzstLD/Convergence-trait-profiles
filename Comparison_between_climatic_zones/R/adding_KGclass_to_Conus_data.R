
## Loading packages
library(sf)
library(maptools)
library(raster)
library(sp)
library(tiff)
library(ijtiff)
library(rgdal)
# Loading data
genus_occ <- read.csv("Conus/Genus_Occurrences.csv")

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
KG_regions <- readGDAL("koppen_geiger/Beck_KG_V1_present_0p083.tif")

# Checking KG data to have the same CRS as Conus data (WGS84)
st_crs(KG_regions) # Yes, it also uses WGS84 as CRS



# # Creating a bounding box around Conus data
# bbox.USA <- st_bbox(all.wgs.sf)
# 
# # Checking if the bounding box fits the desired area of the Koppen 
# # Geiger data
# plot(KG_regions,
#      xlim = bbox.USA[c(1, 3)],
#      ylim = bbox.USA[c(2, 4)]) # Looks good
# 
# 
# # Croping the Koppen Geiger Data to the required bounding box
# KG_USA <- crop(KG_regions.raster, bbox.USA)



# Turning KG_regions from a SpatialGridDataFrame to a Raster object 
# to utilize raster functions
KG_regions.raster <- raster(KG_regions)

# Extracting the Koppen Geiger Index data for each Conus Observation
all.wgs.sf$KoppenGeigerIndex <- raster::extract(KG_regions.raster, all.wgs.sf)

# Loading Koppen Geiger legend data 
KG_legend <- read.csv("koppen_geiger/KG_legend.txt")

# Adding Koppen Geiger descriptions to Conus observations
all.wgs.sf <- base::merge(all.wgs.sf,KG_legend, 
                          by.x = "KoppenGeigerIndex",
                          by.y = "Index",
                          all = T)

# Exporting file
write.csv(all.wgs.sf,"C:\\Users\\euble\\Documents\\Convergence-trait-profiles\\Conus_data_KoppenGeiger.csv")
