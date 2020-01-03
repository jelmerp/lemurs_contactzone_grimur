################################################################################
##### SET-UP #####
################################################################################
setwd('/home/jelmer/Dropbox/sc_lemurs/')
#source('scripts/maps/distmap_fun.R')

library(sf)
library(tidyverse)
library(rworldmap)

## Files:
shapefile_mur <- 'metadata/range_data/mur/data_0.shp'
shapefile_gri <- 'metadata/range_data/gri/data_0.shp'
outfile_map <- 'metadata/maps/map1.pdf'

## Prep sampling locations:
locs <- read.delim('metadata/coords/sites.txt', header = TRUE, as.is = TRUE)

## Set long/lat limits:
longmin <- 42
longmax <- 50
latmin <- -19
latmax <- -26

## Set distribution colours:
col.mur <- '#FF00FF'
col.gri <- '#FFB90F'
loc.size <- 11

## Define a suitable map theme for ggplot:
theme_invis <- theme_bw()
theme_invis$plot.background = element_blank()
theme_invis$panel.background = element_blank()
theme_invis$panel.grid.minor = element_blank()
theme_invis$panel.grid.major = element_blank()
theme_invis$panel.border = element_blank()
theme_invis$axis.line = element_blank()
theme_invis$axis.ticks = element_blank()
theme_invis$axis.text.x = element_blank()
theme_invis$axis.text.y = element_blank()
theme_invis$axis.title.x = element_blank()
theme_invis$axis.title.y = element_blank()
theme_invis$panel.margin = unit(0, 'lines')
theme_invis$plot.margin = unit(c(0, 0, -1, -1), 'lines')
theme_invis$legend.position = 'none'
## Font: # see http://blog.revolutionanalytics.com/2012/09/how-to-use-your-favorite-fonts-in-r-charts.html
#loadfonts() # For ps files instead of pdf: loadfonts(device="postscript")

################################################################################
##### NEW METHODS #####
################################################################################
shapefile.gri <- st_read(shapefile_mur)
shapefile.mur <- st_read(shapefile_gri)
#st_geometry_type(shapefile.gri)
#st_crs(shapefile.gri)
#st_bbox(shapefile.gri)
#shapefile.gri

## Map:
map_raw <- getMap(resolution = 'low') # rworldmap function to get a map
map <- fortify(map_raw) # convert map data into df for ggplot2

p <- ggplot() +
  geom_map(aes(long, lat, map_id = id),
                      colour = 'gray50', # colour defines colour for land and country borders
                      data = map, map = map, fill = 'white') +
  geom_sf(data = shapefile.gri, size = 0, color = "black", fill = col.gri, alpha = 0.8) +
  geom_sf(data = shapefile.mur, size = 0, color = "black", fill = col.mur, alpha = 0.8) +
  coord_sf() +
  theme_invis +
  scale_x_continuous(limits = c(longmin, longmax), expand = c(0, 0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0, 0))
p




################################################################################
##### READ FILES #####
################################################################################

## Load carrion/hooded/orientalis crow geographical distribution:
dist.mur <- read.shapefile(layer = 'data_0', folder = dir_shapefile_mur) %>%
  prep.shapefile(., crop_longlat = FALSE)
dist.gri <- read.shapefile(layer = 'data_0', folder = dir_shapefile_gri) %>%
  prep.shapefile(., crop_longlat = FALSE)

## Prep sampling locations:
locs <- read.delim('metadata/coords/sites.txt', header = TRUE, as.is = TRUE)

# samples <- read.delim('sample_locations/individuals.txt', header = TRUE)
# locs$coloc <- factor(paste(locs$country_code, locs$location_code, sep = '_')) # 'coloc' is the combination of country and location codes
# samples$coloc <- factor(paste(samples$country, samples$population2, sep = '_')) # create the same column in the samples file
# samples <- subset(samples, subset = species == 'C' | species == 'H' | species == 'O' | species == 'P' | species == 'Y',
#                   select = c(coloc, sample_code)) # get the samples of interest
# counts <- data.frame(table(samples$coloc)) # get counts of nr of samples per location (coloc)
# locs <- subset(locs, select = c(lat, long, coloc, location_colour))
# locs$freq <- counts$Freq[match(locs$coloc, counts$Var1)] # align the counts per location with the locations in the locs file, where lat/long are found
# locs <- subset(locs, freq > 0)


################################################################################
##### PLOT MAP ####
################################################################################
## Get world map and "fortify" it:
map_raw <- getMap(resolution = 'low') # rworldmap function to get a map
map <- fortify(map_raw) # convert map data into df for ggplot2

## Plot the base map:
p <- ggplot() # p <- ggplot(map, aes(long, lat))
p <- p + geom_map(aes(long, lat, map_id = id),
                  colour = 'gray50', # colour defines colour for land and country borders
                  data = map, map = map, fill = 'white')
p <- p + theme_invis # add invisible theme suitable for map
p <- p + scale_x_continuous(limits = c(longmin, longmax), expand = c(0, 0)) # long
p <- p + scale_y_continuous(limits = c(latmin, latmax), expand = c(0, 0)) # lat

## Add distributions:
p <- p + geom_polygon(aes(long, lat, group = group), data = dist.mur, fill = col.mur, alpha = 0.8)
p <- p + geom_polygon(aes(long, lat, group = group), data = dist.gri, fill = col.gri, alpha = 0.8)

## Add sampling points:
p <- p + geom_point(data = locs, aes(x = long, y = lat),
                    shape = 21, colour = 'black', fill = 'black', size = 1)
#p <- p + geom_point(aes(x = long, y = lat, fill = location_colour),
#                    data = locs, shape = 21, colour = 'black', size = 7)
#p <- p + scale_fill_identity(breaks = locs$location_colour)

## Add text annotation - sampling locations:
p <- p + annotate(geom = 'text', x = 23.3, y = 45, fontface = 'plain', label = 'cnx2', size = loc.size) # Bulgary
p <- p + geom_segment(aes(x = 89.2, xend = 90.5, y = 56.2, yend = 58.8), colour = 'white') # arrow to Russian hybrids

## Add text annotation - sampling locations:
p <- p + annotate(geom = 'text', x = -9.5, y = 62, fontface = 'bold.italic', label = 'C. c. corone', size = 12)


################################################################################
##### PLOT MAP ####
################################################################################
ggsave(outfile_map, width = 20,
       height = (latmax - latmin) / (longmax - longmin) * 1.5 * 20) # 1.8

system(paste0('xdg-open ', outfile_map))



################################################################################
##### OTHER PLOTTING OPTIONS ####
################################################################################
#p <- p + theme(panel.background = element_rect(fill = 'lightblue1')) # colour for waterbodies
#p <- p + coord_equal(ratio = (latmax - latmin) / (longmax - longmin)) # edit aspect ratio to a more realistic one
#p <- p + scale_size(range = c(3, 10)) # bounds for point size variation

# ## Mercator projections - can't get them to work
#map_raw <- spTransform(map_raw, CRS("+proj=merc"))
#map_raw <- map(database = "world", projection = "mercator", xlim = c(longmin, longmax), ylim = c(latmin, latmax))

# ## Select only certain countries - problem is very approximate outlines:
# countries <- read.delim('sample_locations/country_codes.txt', header = TRUE)
# map_data <- joinCountryData2Map(countries, joinCode = "ISO3", nameJoinColumn = "code")
# map_raw <- map_data[map_data$NAME %in% countries$country,]
# map <- fortify(map_raw) # convert map data into df for ggplot2

# ## Plot conically:
# pc <- p + coord_map('conic', lat0 = mean(latmin, latmax))
# p <- p + scale_x_continuous(limits = c(0, 160), expand = c(0, 0)) # long
# p <- p + scale_y_continuous(limits = c(25, 70), expand = c(0, 0)) # lat
# ggsave(plot = pc, file = 'crowmap_conic.png', width = 20, height = (latmax - latmin) / (longmax - longmin) * 2 * 20)
# system('xdg-open crowmap_conic.png')

# ## Plot using mercator projection: # Doesnt work, creates lines through map
# pm <- p + coord_map('mercator')
# ggsave(plot = pm, file = 'crowmap_mercator1.png', width = 20, height = 20)
# ggsave(plot = pm, file = 'crowmap_mercator1.png', width = 20, height = (latmax - latmin) / (longmax - longmin) * 1.8 * 20)
# system('xdg-open crowmap_mercator1.png')

# ## Plotting in ggmap with google maps:
# map <- get_googlemap(center = c(lon = 79.8, lat = 45.0), zoom = 3)
# attr(map, 'bb')[1:4] <- c('ll.lat' = 5, 'll.lon' = -10, 'ur.lat' = 70, 'ur.lon' = 130)
# p <- ggmap(map, extent = 'device')
#
# map <- get_googlemap(center = 'Moscow', maptype = 'terrain', zoom = 2, color = 'bw', style = "style=feature:all|element:labels|visibility:off")
# attr(map, 'bb')[1:4] <- c('ll.lat' = 5, 'll.lon' = -10, 'ur.lat' = 70, 'ur.lon' = 130)
# p <- ggmap(map,  extent = 'device')
# p <- p + geom_point(aes(x = long, y = lat, fill = location_colour), data = locs, shape = 21, colour = 'black', size = 7)
# p <- p + geom_polygon(aes(long, lat, group = group), data = dist.cx)
