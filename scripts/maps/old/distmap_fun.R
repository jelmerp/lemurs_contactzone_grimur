## Function to load/install packages:
getPackage <- function(pkg) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
  return(TRUE)
}

## Load packages:
packages <- c('ggplot2', 'plyr', 'maps', 'rworldmap', 'grid', 'sp', 'rgdal',
              'rgeos', 'ggmap', 'maptools', 'extrafont')
sapply(packages, getPackage)
# install.packages("gpclib", type="source")
# gpclibPermit()


## Function to read a shapefile:
read.shapefile <- function(layer, folder = './shapefiles/', to.mercator = FALSE, read.method = 'rgdal') {

  if(read.method == 'rgdal') { # use rgdal function readOGR to read in shapefile
    dist <- readOGR(dsn = folder, layer = layer)
  }
  if(read.method == 'maptools') { # use maptools function readShapeSpatial to read in shapefile
    dist <- readShapeSpatial(fn = paste0(folder, layer))
  }

  #proj4string(dist) # check projection
  if(to.mercator == TRUE) dist <- spTransform(dist, CRS("+proj=merc"))

  return(dist)
}

## Function to prepare shapefile for plotting:
prep.shapefile <- function(shapefile, crop_longlat = FALSE) {
  # shapefile <- dist.mur

  ## See https://github.com/hadley/ggplot2/wiki/plotting-polygon-shapefiles for following steps

  shapefile@data$id <- rownames(shapefile@data)
  shapefile.points <- fortify(shapefile, region = 'id')
  shapefile.df <- join(shapefile.points, shapefile@data, by = 'id')
  #shapefile.df <- subset(shapefile.df, SEASONAL != 3) # get rid of wintering grounds

  if(crop_longlat == TRUE) {
    shapefile.df <- shapefile.df %>%
      dplyr::filter(long > longmin & long < longmax & lat > latmin & lat < latmax)
  }

  return(shapefile.df)
}