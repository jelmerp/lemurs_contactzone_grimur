################################################################################
##### SET-UP #####
################################################################################
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')
library(ggmap)
library(rgdal) # readOGR
library(tidyverse)
library(ggsn) #devtools::install_github('oswaldosantos/ggsn')
library(ggrepel)
library(ggpubr)
library(cowplot)

## Files:
shapefile_dir_gri <- '../metadata/range_data/mur/'
shapefile_dir_mur <- '../metadata/range_data/gri/'
infile_locs <- '../metadata/coords/sites.txt'
infile_coords <- 'metadata/gps/Andohahela_coords_all.txt'

outfile_map <- 'finalfigures/Fig1_map.eps'

## Set distribution colours:
col.mur <- '#FF00FF'
col.gri <- '#D69B12' #'#FFB90F'


################################################################################
##### OVERVIEW MAP #####
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
dist_map <- get_googlemap(center = c(lon = 46, lat = -23),
                     zoom = 7, scale = 2, # 3=continent, 10=city, 21=building
                     maptype = 'satellite', # roadmap, terrain, hybrid, satellite
                     color = 'color')

## Lat and lon:
lon.min <- 43.1
lon.max <- 48.7
lat.min <- -25.7
lat.max <- -19.7

## Create plot:
p_dist <- ggmap(dist_map, extent = 'device') +
  geom_polygon(aes(x = long, y = lat, group = id),
               data = shapedata_gri,
               fill = col.gri, alpha = 0.5, size = 0) +
  geom_polygon(aes(x = long, y = lat, group = id),
               data = shapedata_mur,
               fill = col.mur, alpha = 0.3, size = 0) +
  geom_point(data = locs,
             aes(x = long, y = lat, fill = hzproj_type, shape = hzproj_type2),
             colour = 'black',
             size = 12,
             stroke = 6) +
  scale_x_continuous(limits = c(lon.min, lon.max), expand = c(0, 0)) +
  scale_y_continuous(limits = c(lat.min, lat.max), expand = c(0, 0)) +
  scale_shape_manual(values = c(8, 21), guide = FALSE) +
  scale_fill_manual(values = c('black', '#FFB90F', '#FF00FF'),
                    labels = c('hybridzone', 'griseorufus', 'murinus'),
                    guide = FALSE) +
  theme_void() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'lines'),
        panel.border = element_rect(colour = "grey20", fill = NA, size = 1)) +
  ggsn::scalebar(x.min = lon.min + 0.15,
                 x.max = lon.max,
                 y.min = lat.min + 0.25,
                 y.max = lat.max,
                 dist = 50,
                 dist_unit = 'km',
                 location = 'bottomleft',
                 transform = TRUE,
                 model = 'WGS84',
                 st.size = 9,
                 border.size = 1.5,
                 box.fill = c("gray60", "white"),
                 st.color = "white")
#p_dist


################################################################################
##### CONTACT ZONE MAP #####
################################################################################
## Prep sampling locations:
coords <- read.delim(infile_coords, as.is = TRUE)
sp.cols <- c('#D69B12', '#FF00FF')
loclabel.size <- 14

## Coord center and bounds:
lat.center <- -24.88038
lat.min <- -24.97801
lat.max <- -24.78632
lon.center <- 46.61519
lon.min <- 46.50373
lon.max <- 46.7207

cz_map <- get_googlemap(center = c(lon = lon.center, lat = lat.center),
                        zoom = 10, scale = 2,
                        maptype = 'terrain',
                        color = 'color')

p_cz <- ggmap(cz_map, extent = 'device') +
  geom_point(data = coords,
             aes(x = lon, y = lat, fill = species.cor),
             colour = 'black',
             size = 9,
             stroke = 1,
             shape = 21) +
  scale_x_continuous(limits = c(lon.min, lon.max), expand = c(0, 0)) +
  scale_y_continuous(limits = c(lat.min, lat.max), expand = c(0, 0)) +
  scale_fill_manual(values = sp.cols) +
  theme_void() +
  theme(legend.position = 'none',
        plot.margin = unit(c(2, 0.4, 0.8, 0.4), 'lines'),
        panel.border = element_rect(colour = "grey1", fill = NA, size = 1)) +
  ggsn::scalebar(x.min = lon.min,
                 x.max = lon.max - 0.015,
                 y.min = lat.min + 0.012,
                 y.max = lat.max,
                 dist = 2,
                 dist_unit = 'km',
                 location = 'bottomright',
                 transform = TRUE,
                 model = 'WGS84',
                 st.size = 10,
                 border.size = 1.5,
                 box.fill = c("gray60", "white"),
                 st.color = "black") +
  annotate(geom = 'text', x = 46.615, y = -24.795,
           fontface = 'bold', size = 16,
           label = 'contact zone area') +
  annotate(geom = 'text', x = 46.567, y = -24.85,
           fontface = 'plain', size = loclabel.size,
           label = 'Hazofotsy\n(parapatric)') +
  annotate(geom = 'text', x = 46.67, y = -24.84,
           fontface = 'plain', size = loclabel.size,
           label = 'Ambatoabo\n(parapatric)') +
  annotate(geom = 'text', x = 46.553, y = -24.938,
           fontface = 'plain', size = loclabel.size,
           label = 'Mangatsiaka\n(sympatric)') +
  annotate(geom = 'text', x = 46.66, y = -24.935,
           fontface = 'plain', size = loclabel.size,
           label = 'Tsimelahy\n(sympatric)')


################################################################################
##### SAVE MAP ####
################################################################################
p <- ggarrange(p_dist, p_cz,
               ncol = 1, nrow = 2, heights = c(1, 0.6)) +
  draw_line(x = c(0.273, 0.56), y = c(0.358, 0.45), colour = 'grey1', size = 1.5) +
  draw_line(x = c(0.728, 0.598), y = c(0.358, 0.45), colour = 'grey1', size = 1.5)

ggsave(outfile_map, width = 19, height = 24,
       device = cairo_ps, fallback_resolution = 150)
system(paste0('xdg-open ', outfile_map))



################################################################################
################################################################################
#outfile_map_eps_highres <- paste0(outfile_map_eps, '.highres.eps')
#ggsave(outfile_map_eps_highres, width = 20, height = 15,
#       device = cairo_ps, fallback_resolution = 600)

#p.dist + ggsn::north2(p, x = 0.25, y = 0.9, symbol = 4, scale = 0.05)

#p <- p + geom_segment(aes(x = 89.2, xend = 90.5, y = 56.2, yend = 58.8),
#                      colour = 'white') # arrow

#geom_label_repel(data = locs, colour = 'black', size = 5, #point.padding = 5,
#           aes(x = long, y = lat, label = site)) +