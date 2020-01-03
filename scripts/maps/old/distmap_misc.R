#library(sf)
#library(rworldmap)

## Install ggmap:
if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("dkahle/ggmap", ref = "tidyup", force=TRUE)

## Google map API:
ggmap::register_google(key = "AIzaSyDbUxCT8W_YxLhpPsDT25mJQK7aOPPy1Cg")

## Shapefile info:
st_geometry_type(shapefile.gri)
st_crs(shapefile.gri)
st_bbox(shapefile.gri)
shapefile.gri

## Underlying map:
map <- getMap(resolution = 'low') %>% fortify() # rworldmap function to get a map

## Set long/lat limits:
longmin <- 42
longmax <- 50
latmin <- -26
latmax <- -19
myLocation <- c(42, -26, 50, -19) # Zoomed out
myLocation <- c(46.55039, -24.97052, 46.56164, -24.95795) # Zoomed in
my.map <- get_map(location = myLocation, source = "stamen",
                  maptype = "terrain", crop = FALSE)
my.map <- get_map(location = myLocation, source = "osm")

## Shapefiles:
shapefile_gri <- 'metadata/range_data/gri/data_0.shp'
shapefile_mur <- 'metadata/range_data/mur/data_0.shp'
shapefile.gri <- st_read(shapefile_mur)
shapefile.mur <- st_read(shapefile_gri)

p <- ggplot() +
  geom_map(aes(long, lat, map_id = id),
           colour = 'gray50', # colour defines colour for land and country borders
           data = map, map = map, fill = 'white') +
  geom_sf(data = shapefile.gri, size = 0, color = "black", fill = col.gri, alpha = 0.8) +
  geom_sf(data = shapefile.mur, size = 0, color = "black", fill = col.mur, alpha = 0.8) +
  coord_sf() +
  #theme_invis +
  geom_point(data = locs, aes(x = long, y = lat),
             shape = 21, colour = 'black', fill = 'black', size = 3) +
  scale_x_continuous(limits = c(longmin, longmax), expand = c(0, 0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0, 0))
p

shapedata_gri <- readOGR(dsn = shapefile_dir_gri, layer = "data_0")
shapedata_mur <- readOGR(dsn = shapefile_dir_mur, layer = "data_0")
#proj4string(shapedata_gri) # Check current projection
#shpData <- spTransform(shpData, CRS("+proj=longlat +datum=WGS84"