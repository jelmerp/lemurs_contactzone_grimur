##########################################################################
##### SET-UP #####
##########################################################################
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')
library(ggmap)
library(rgdal) # readOGR()
library(tidyverse)
library(ggsn) #devtools::install_github('oswaldosantos/ggsn')
library(ggpubr)

## Files:
infile_coords_Mtk <- 'metadata/gps/Mangatsiaka_coords_dec.txt'
infile_coords_Tsm <- 'metadata/gps/Tsimelahy_coords_dec.txt'
plotdir <- '/home/jelmer/Dropbox/sc_lemurs/presentations/posters/2019_03_GordonConf/plots/'
outfile_map_pdf <- paste0(plotdir, '/HZmap1.pdf')
outfile_map_png <- paste0(plotdir, '/HZmap1.png')
outfile_map_eps <- paste0(plotdir, '/HZmap1.eps')

## Prep sampling locations:
coords_Mtk <- read.delim(infile_coords_Mtk, header = TRUE, as.is = TRUE)
coords_Tsm <- read.delim(infile_coords_Tsm, header = TRUE, as.is = TRUE)
#coords_Mtk %>% filter(species.short == 'hybrid')

## Labels:
sp.labs <- c(expression(italic("griseorufus")), expression(italic("microcebus")))
hyb.labs <- c('hybrid', 'non-hybrid')
sp.cols <- c('#D69B12', '#FF00FF')
msat.name <- 'previous\nmicrosatellite\nassignment:'
radseq.name <- '\nRADseq\nassignment:'

## Mangatsiaka lat and lon:
lat.center.Mtk <- -24.96275
lon.center.Mtk <- 46.55692
lon.min.Mtk <- 46.55267
lon.max.Mtk <- 46.56099
lat.min.Mtk <- -24.96914 #-24.96686
lat.max.Mtk <- -24.96000
lon.diff.Mtk <- 46.56099 - 46.55267
lon.min.Mtk.scalebar <- lon.min.Mtk + (0.05 * (lon.max.Mtk - lon.min.Mtk))
lat.max.Mtk.scalebar <- lat.max.Mtk - (0.04 * (lat.max.Mtk - lat.min.Mtk))
lon.max.Mtk.scalebar <- lon.max.Mtk
lat.min.Mtk.scalebar <- lat.min.Mtk

## Tsimelahy lat and lon:
lat.center.Tsm <- -24.95425
lon.center.Tsm <- 46.61117
lon.min.Tsm <- 46.61004
lon.max.Tsm <- 46.62222
lat.min.Tsm <- -24.95884
lat.max.Tsm <- -24.94711
lon.min.Tsm.scalebar <- lon.min.Tsm + (0.05 * (lon.max.Tsm - lon.min.Tsm))
lat.max.Tsm.scalebar <- lat.max.Tsm - (0.04 * (lat.max.Tsm - lat.min.Tsm))
lon.max.Tsm.scalebar <- lon.max.Tsm
lat.min.Tsm.scalebar <- lat.min.Tsm

##########################################################################
##### MTK PLOT #####
##########################################################################
## Google map:
ggmap::register_google(key = "AIzaSyB0hxxyLujH5apnsRW21X5MlqO8_LnbjGc")
map.Mtk <- get_googlemap(center = c(lon = lon.center.Mtk, lat = lat.center.Mtk),
                        zoom = 15, scale = 2, # 3=continent, 10=city, 21=building
                        maptype = 'satellite', # roadmap, terrain, hybrid, satellite
                        color = 'color')

## Plot:
p.Mtk <- ggmap(map.Mtk, extent = 'device') +
  geom_point(data = coords_Mtk, colour = 'black', size = 4, stroke = 1,
             aes(x = lon, y = lat, shape = msat.hybrid, fill = species.cor)) +
  scale_x_continuous(limits = c(lon.min.Mtk, lon.max.Mtk), expand = c(0, 0)) +
  scale_y_continuous(limits = c(lat.min.Mtk, lat.max.Mtk), expand = c(0, 0)) +
  scale_fill_manual(values = sp.cols) +
  scale_shape_manual(values = c(25, 21)) +
  labs(title = 'contact zone - site 1\n(Mangatsiaka)') +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = 'plain'),
        legend.position = 'none',
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'lines'),
        panel.border = element_rect(colour = "grey20", fill = NA, size = 2)) +
  ggsn::scalebar(x.min = lon.min.Mtk.scalebar, x.max = lon.max.Mtk.scalebar,
                 y.min = lat.min.Mtk.scalebar, y.max = lat.max.Mtk.scalebar,
                 dist = 100, dist_unit = 'm', location = 'topleft',
                 transform = TRUE, model = 'WGS84', st.size = 2, border.size = 0.8,
                 box.fill = c("gray60", "white"), st.color = "white")
#p.Mtk

##########################################################################
##### TSM PLOT #####
##########################################################################
## Google map:
map.Tsm <- get_googlemap(center = c(lon = lon.center.Tsm, lat = lat.center.Tsm),
                        zoom = 15, scale = 2, # 3=continent, 10=city, 21=building
                        maptype = 'satellite', # roadmap, terrain, hybrid, satellite
                        color = 'color')

## Plot:
p.Tsm <- ggmap(map.Tsm, extent = 'device') +
  geom_point(data = coords_Tsm, colour = 'black', size = 4, stroke = 1,
             aes(x = lon, y = lat, shape = msat.hybrid, fill = species.cor)) +
  scale_x_continuous(limits = c(lon.min.Tsm, lon.max.Tsm), expand = c(0, 0)) +
  scale_y_continuous(limits = c(lat.min.Tsm, lat.max.Tsm), expand = c(0, 0)) +
  scale_shape_manual(values = c(25, 21), labels = hyb.labs, name = msat.name) +
  scale_fill_manual(values = sp.cols, labels = sp.labs, name = radseq.name) +
  guides(fill = guide_legend(override.aes = list(shape = 21),
                             label.hjust = -1)) +
  labs(title = 'contact zone - site 2\n(Tsimelahy)') +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = 'plain'),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = 'bold'),
        legend.key.size = unit(0.6, "cm"),
        legend.position = 'right',
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-5, 0, -5, -5),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'lines'),
        panel.border = element_rect(colour = "grey20", fill = NA, size = 2)) +
  scalebar(x.min = lon.min.Tsm.scalebar, x.max = lon.max.Tsm.scalebar,
           y.min = lat.min.Tsm.scalebar, y.max = lat.max.Tsm.scalebar,
           dist = 100, dist_unit = 'm', location = 'topleft',
           transform = TRUE, model = 'WGS84', st.size = 2, border.size = 0.8,
           box.fill = c("gray60", "white"), st.color = "white")
#p.Tsm

##########################################################################
##### SAVE MAP ####
##########################################################################
p <- ggarrange(p.Mtk, p.Tsm, ncol = 2, nrow = 1, widths = c(1, 1.673))
ggexport(p, filename = outfile_map_eps)
system(paste0('xdg-open ', outfile_map_eps))

#ggexport(p, filename = outfile_map_pdf)
#ggexport(p, filename = outfile_map_png)