################################################################################
##### SET-UP #####
################################################################################
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')
library(ggmap)
#library(rgdal) # readOGR() #install.packages('rgdal')
library(tidyverse)
library(ggsn) #devtools::install_github('oswaldosantos/ggsn')
library(ggpubr)

## Files:
infile_coords <- 'metadata/gps/Andohahela_coords_all.txt'
plotdir <- '/home/jelmer/Dropbox/sc_lemurs/hybridzone/maps/'
outfile_map_eps <- paste0(plotdir, '/CZmap.eps')

## Prep sampling locations:
coords <- read.delim(infile_coords, as.is = TRUE)

## Labels:
sp.cols <- c('#D69B12', '#FF00FF')

## Coord center and bounds:
lat.center <- -24.88038
lat.min <- -24.98906
lat.max <- -24.78632
lon.center <- 46.61519
lon.min <- 46.51575
lon.max <- 46.70765

## Register API key for Google maps:
ggmap::register_google(key = "AIzaSyB0hxxyLujH5apnsRW21X5MlqO8_LnbjGc")


################################################################################
##### MAP #####
################################################################################
cz_map <- get_googlemap(center = c(lon = lon.center, lat = lat.center),
                        zoom = 10, scale = 2, # zoom:3=continent, 10=city, 21=building
                        maptype = 'terrain', # roadmap, terrain, hybrid, satellite
                        color = 'color')

cz <- ggmap(cz_map, extent = 'device') +
  geom_point(data = coords,
             aes(x = lon, y = lat, fill = species.cor),
             colour = 'black',
             size = 3,
             stroke = 1,
             shape = 21) +
  scale_x_continuous(limits = c(lon.min, lon.max), expand = c(0, 0)) +
  scale_y_continuous(limits = c(lat.min, lat.max), expand = c(0, 0)) +
  scale_fill_manual(values = sp.cols) +
  theme(legend.position = 'none',
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'lines'),
        panel.border = element_rect(colour = "grey20", fill = NA, size = 1)) +
  ggsn::scalebar(x.min = lon.min + 0.01,
                 x.max = lon.max,
                 y.min = lat.min,
                 y.max = lat.max - 0.01,
                 dist = 2,
                 dist_unit = 'km',
                 location = 'topleft',
                 transform = TRUE,
                 model = 'WGS84',
                 st.size = 4,
                 border.size = 1,
                 box.fill = c("gray40", "gray20"),
                 st.color = "black")
cz

## Save map:
#p <- ggarrange(p.Mtk, p.Tsm, ncol = 2, nrow = 1, widths = c(1, 1.673))
#ggexport(p, filename = outfile_map_msats_eps)
#system(paste0('xdg-open ', outfile_map_msats_eps))

