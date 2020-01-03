################################################################################
##### SET-UP #####
################################################################################
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')
library(ggmap)
library(rgdal) # readOGR
library(tidyverse)
library(ggsn) #devtools::install_github('oswaldosantos/ggsn')
library(ggrepel)

## Files:
shapefile_dir_gri <- '../metadata/range_data/mur/'
shapefile_dir_mur <- '../metadata/range_data/gri/'
infile_locs <- '../metadata/coords/sites.txt'

#outfile_map_pdf <- 'maps/distmap.pdf'
#outfile_map_png <- 'maps/distmap.png'
outfile_map <- 'maps/distmap.eps'
outfile_map_distOnly <- 'maps/distmap_distOnly.eps'

## Set distribution colours:
col.mur <- '#FF00FF'
col.gri <- '#D69B12' #'#FFB90F'
loc.size <- 11


################################################################################
##### READ SHAPEFILE, LOAD BASE MAP, PREP SAMPLE LOCS #####
################################################################################
## Sampling locations:
locs <- read.delim(infile_locs, as.is = TRUE) %>%
  filter(site != 'Tsimelahy',
         site != 'Hazofotsy',
         site != 'Ambatoabo')

## Shapefiles:
shapedata_gri <- readOGR(dsn = shapefile_dir_gri, layer = "data_0")
shapedata_mur <- readOGR(dsn = shapefile_dir_mur, layer = "data_0")

## Google map:
ggmap::register_google(key = "AIzaSyB0hxxyLujH5apnsRW21X5MlqO8_LnbjGc")
my.map <- get_googlemap(center = c(lon = 46, lat = -23),
                     zoom = 7, scale = 2, # 3=continent, 10=city, 21=building
                     maptype = 'satellite', # roadmap, terrain, hybrid, satellite
                     color = 'color')

## Lat and lon:
lon.min <- 43.1
lon.max <- 48.7
lat.min <- -25.7
lat.max <- -19.7


################################################################################
##### CREATE PLOT #####
################################################################################
p.dist <- ggmap(my.map, extent = 'device') +
  geom_polygon(aes(x = long, y = lat, group = id),
               data = shapedata_gri,
               fill = col.gri, alpha = 0.5, size = 0) +
  geom_polygon(aes(x = long, y = lat, group = id),
               data = shapedata_mur,
               fill = col.mur, alpha = 0.3, size = 0) +
  geom_point(aes(x = long, y = lat,
                 fill = hzproj_type, shape = hzproj_type2),
             data = locs, colour = 'black', size = 12, stroke = 6) +
  scale_x_continuous(limits = c(lon.min, lon.max), expand = c(0, 0)) +
  scale_y_continuous(limits = c(lat.min, lat.max), expand = c(0, 0)) +
  scale_shape_manual(values = c(8, 21), guide = FALSE) +
  scale_fill_manual(values = c('black', '#FFB90F', '#FF00FF'),
                    labels = c('hybridzone', 'griseorufus', 'murinus'),
                    guide = FALSE) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'lines'),
        panel.border = element_rect(colour = "grey20", fill = NA, size = 1)) +
  ggsn::scalebar(x.min = lon.min,
                 x.max = lon.max - 0.4,
                 y.min = lat.min + 0.3,
                 y.max = lat.max,
                 dist = 50,
                 dist_unit = 'km',
                 location = 'bottomright',
                 transform = TRUE,
                 model = 'WGS84',
                 st.size = 9,
                 border.size = 1.5,
                 box.fill = c("gray60", "white"),
                 st.color = "white")
p.dist


################################################################################
##### CREATE PLOT #####
################################################################################
p.dist2 <- ggmap(my.map, extent = 'device') +
  geom_polygon(aes(x = long, y = lat, group = id),
               data = shapedata_gri,
               fill = col.gri, alpha = 0.6, size = 0) +
  geom_polygon(aes(x = long, y = lat, group = id),
               data = shapedata_mur,
               fill = col.mur, alpha = 0.3, size = 0) +
  scale_x_continuous(limits = c(lon.min, lon.max), expand = c(0, 0)) +
  scale_y_continuous(limits = c(lat.min, lat.max), expand = c(0, 0)) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'lines'),
        panel.border = element_rect(colour = "grey20", fill = NA, size = 1)) +
  ggsn::scalebar(x.min = lon.min, x.max = lon.max - 0.4,
                 y.min = lat.min + 0.3, y.max = lat.max,
                 dist = 50, dist_unit = 'km', location = 'bottomright',
                 transform = TRUE, model = 'WGS84',
                 st.size = 9, border.size = 1.5,
                 box.fill = c("gray60", "white"), st.color = "white")
p.dist2


################################################################################
##### SAVE MAP ####
################################################################################
ggsave(outfile_map, width = 20, height = 15,
       device = cairo_ps, fallback_resolution = 150)
system(paste0('xdg-open ', outfile_map))

ggsave(outfile_map_distOnly, width = 20, height = 15,
       device = cairo_ps, fallback_resolution = 150)
system(paste0('xdg-open ', outfile_map_distOnly))

#outfile_map_eps_highres <- paste0(outfile_map_eps, '.highres.eps')
#ggsave(outfile_map_eps_highres, width = 20, height = 15,
#       device = cairo_ps, fallback_resolution = 600)
#ggsave(outfile_map_pdf, width = 20, height = 15)
#ggsave(outfile_map_png, width = 20, height = 15)

################################################################################
################################################################################
#p.dist + ggsn::north2(p, x = 0.25, y = 0.9, symbol = 4, scale = 0.05)

#geom_label_repel(data = locs, colour = 'black', size = 5, #point.padding = 5,
#           aes(x = long, y = lat, label = site)) +

## Add text annotation - sampling locations:
#p <- p + annotate(geom = 'text', x = 23.3, y = 45, fontface = 'plain',
#                  label = 'cnx2', size = loc.size) # Bulgary
#p <- p + geom_segment(aes(x = 89.2, xend = 90.5, y = 56.2, yend = 58.8),
#                      colour = 'white') # arrow