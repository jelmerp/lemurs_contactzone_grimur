################################################################################
##### SET-UP - GENERAL #####
################################################################################
## Libraries and scripts:
library(tidyverse)
library(ggpubr)
library(cowplot)
library(here)

## Input files:
infile_inds <- here('metadata/radseq_metadata_link/lookup_IDshort.txt')
infile_cols <- '/home/jelmer/Dropbox/sc_lemurs/metadata/colors/colors.species.txt'

## Output files:
figdir <- here('analyses/admixture/figs_final/')
figfile_eps <- file.path(figdir, '/Fig_structure.eps')

## Read metadata:
cols.sp.df <- read.delim(infile_cols, header = T, as.is = T) %>%
  select(-species)
inds.df <- read.delim(infile_inds, header = T, as.is = T) %>%
  merge(., cols.sp.df, by = 'species.short') %>%
  rename(labcol = color)

col.mur <- '#FF00FF'
col.gri <- '#D69B12'
col.hyb <- 'grey50'


################################################################################
##### SET-UP - ADMIXTURE #####
################################################################################
source('/home/jelmer/Dropbox/scripts/genomics/admixture/admixture_plot_fun.R')
indir.admix <- 'analyses/admixture/output/'
setID <- 'r03.all.mac3.FS6'

## Sympatric vs parapatric sites:
inds.df$site.type[inds.df$site.type == 'parapatric' &
                    inds.df$species.short == 'mmur'] <- 'zparapatric'
#table(inds.df$site.type)

## Shorter species abbrev:
inds.df$species.short2 <- gsub('^m', '', inds.df$species.short)
inds.df$species.short2 <- gsub('hyb', 'hybrid?', inds.df$species.short2)


################################################################################
##### ADMIXTURE PLOTS #####
################################################################################
## CV-plot:
(p_cv <- k.plot(setID) +
  labs(y = 'CV error') +
  theme(axis.title.x = element_text(margin = unit(c(0.1, 0, 0, 0), 'cm')),
        axis.title.y = element_text(margin = unit(c(0, 0.1, 0, 0), 'cm')),
        plot.title = element_blank(),
        plot.margin = margin(1, 0.2, 0.5, 0.2, 'cm')))

## Q-plot:
my.grouplab.bgcol <- c(col.gri, col.gri, col.hyb, col.mur, col.mur)
site.type.labs <- c(parapatric = 'para', sympatric = 'sym', zparapatric = 'para')
grouplab.labeller <- labeller(site.type = site.type.labs)

p_k2 <- Qdf(setID,
            indir.Qdf = indir.admix,
            K = 2,
            sort.by = 'species.short2') %>%
  ggax.v(.,
         group.column = c('species.short2', 'site.type'),
         barcols = c(col.gri, col.mur),
         grouplab.bgcol = my.grouplab.bgcol,
         grouplab.labeller = grouplab.labeller,
         ylab = 'RADseq')

(p_mt <- Qdf(setID,
            indir.Qdf = indir.admix,
            K = 2,
            sort.by = 'species.short2') %>%
    ggax.v(.,
           group.column = c('species.short2', 'site.type'),
           barcols = c(col.gri, col.mur),
           ylab = 'mtDNA') +
    theme(strip.background = element_blank(),
          strip.text = element_blank()))

p_k <- ggarrange(p_k2, p_mt, ncol = 1, nrow = 2, heights = c(1, 0.25))
p_admix <- ggarrange(p_cv, p_k, nrow = 2, heights = c(0.4, 0.6))


################################################################################
##### SET-UP - PCA #####
################################################################################
source('/home/jelmer/Dropbox/scripts/genomics/PCA/PCA_R_fun.R')

fileID <- 'r03.all.mac1.FS6'
ID.type <- 'ID.short'
keep.set <- 'all'
pca.df.dir <- 'analyses/PCA/dfs/'
subset.ID = 'all'

infile_pca.df <- paste0(pca.df.dir, '/', fileID, '_', subset.ID, '.txt')
infile_eigenvalues <- paste0(pca.df.dir, '/', fileID, '_',
                             subset.ID, '_eigenvalues.txt')

sp.levels <- c('mgri', 'mmur', 'mhyb')
pca.df <- read.table(infile_pca.df, header = TRUE) %>%
  mutate(species.short = factor(species.short, levels = sp.levels)) %>%
  arrange(species.short)

eig.df <- read.table(infile_eigenvalues, header = TRUE)


################################################################################
##### PCA PLOT #####
################################################################################
col.by.labs <- c(expression(italic("griseorufus")),
                 expression(italic("murinus")), 'hybrid?')
shape.by.labs <- c('parapatric', 'sympatric')

(p_pca <- pcplot(pca.df,
              eigenvalues = eig.df$eig,
              col.by = 'species.short',
              col.by.name = 'microsats:',
              col.by.labs = col.by.labs,
              cols = c(col.gri, col.mur, col.hyb),
              shape.by = 'site.type',
              shape.by.name = 'site type:',
              shape.by.labs = shape.by.labs,
              shapes = c(1, 0),
              legpos = 'top',
              plot.save = FALSE) +
  theme(legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = margin(0.6, 0.2, 0, 1, 'cm')))


################################################################################
##### MAP - SET-UP #####
################################################################################
library(ggmap)
library(ggsn) #devtools::install_github('oswaldosantos/ggsn')

## Files:
infile_coords_Mtk <- 'metadata/gps/Mangatsiaka_coords_dec.txt'
infile_coords_Tsm <- 'metadata/gps/Tsimelahy_coords_dec.txt'

## Prep sampling locations:
coords_Mtk <- read.delim(infile_coords_Mtk, header = TRUE, as.is = TRUE)
coords_Tsm <- read.delim(infile_coords_Tsm, header = TRUE, as.is = TRUE)

## Labels:
sp.labs <- c(expression(italic("griseorufus")), expression(italic("   murinus")))
hyb.labs <- c('hybrid', 'non-hybrid')
sp.cols <- c('#D69B12', '#FF00FF')
msat.name <- '\nmicrosat\nassignment:'
radseq.name <- '\nRADseq\nassignment:'

## Mangatsiaka lat and lon:
lat.center.Mtk <- -24.96275
lon.center.Mtk <- 46.55692
lon.min.Mtk <- 46.55267
lon.max.Mtk <- 46.56099
lat.min.Mtk <- -24.96914
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

ggmap::register_google(key = "AIzaSyB0hxxyLujH5apnsRW21X5MlqO8_LnbjGc")


################################################################################
##### PLOT MAPS #####
################################################################################
## Mtk:
map_Mtk <- get_googlemap(center = c(lon = lon.center.Mtk, lat = lat.center.Mtk),
                         zoom = 15, # 3=continent, 10=city, 21=building
                         scale = 2,
                         maptype = 'satellite', # roadmap, terrain, hybrid, satellite
                         color = 'color')

p_Mtk <- ggmap(map_Mtk, extent = 'device') +
  geom_point(data = coords_Mtk,
             aes(x = lon, y = lat, shape = msat.hybrid, fill = species.cor),
             colour = 'black',
             size = 4,
             stroke = 1) +
  scale_x_continuous(limits = c(lon.min.Mtk, lon.max.Mtk), expand = c(0, 0)) +
  scale_y_continuous(limits = c(lat.min.Mtk, lat.max.Mtk), expand = c(0, 0)) +
  scale_fill_manual(values = sp.cols) +
  scale_shape_manual(values = c(25, 21)) +
  labs(title = 'sympatric site 1: Mangatsiaka') +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = 'plain'),
        legend.position = 'none',
        panel.border = element_rect(colour = "grey20", fill = NA, size = 2),
        plot.margin = unit(c(1, 0.2, 0.2, 0.2), 'cm')) +
  ggsn::scalebar(x.min = lon.min.Mtk.scalebar,
                 x.max = lon.max.Mtk.scalebar,
                 y.min = lat.min.Mtk.scalebar,
                 y.max = lat.max.Mtk.scalebar,
                 dist = 150,
                 dist_unit = 'm',
                 location = 'topleft',
                 transform = TRUE,
                 model = 'WGS84',
                 st.size = 5,
                 border.size = 0.8,
                 box.fill = c("gray60", "white"),
                 st.color = "white")

## Tsm:
map_Tsm <- get_googlemap(center = c(lon = lon.center.Tsm, lat = lat.center.Tsm),
                         zoom = 15,
                         scale = 2,
                         maptype = 'satellite',
                         color = 'color')

p_Tsm <- ggmap(map_Tsm, extent = 'device') +
  geom_point(data = coords_Tsm,
             aes(x = lon, y = lat, shape = msat.hybrid, fill = species.cor),
             colour = 'black',
             size = 4,
             stroke = 1) +
  scale_x_continuous(limits = c(lon.min.Tsm, lon.max.Tsm),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(lat.min.Tsm, lat.max.Tsm),
                     expand = c(0, 0)) +
  scale_shape_manual(values = c(25, 21),
                     labels = hyb.labs,
                     name = msat.name) +
  scale_fill_manual(values = sp.cols,
                    labels = sp.labs,
                    name = radseq.name) +
  guides(fill = guide_legend(override.aes = list(shape = 21),
                             label.hjust = -1)) +
  labs(title = 'sympatric site 2: Tsimelahy') +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = 'plain'),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = 'bold'),
        legend.key.size = unit(0.6, "cm"),
        legend.position = 'right',
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-5, 0, -5, -5),
        panel.border = element_rect(colour = "grey20", fill = NA, size = 2),
        plot.margin = unit(c(1, 0, 0.2, 0.2), 'cm')) +
  scalebar(x.min = lon.min.Tsm.scalebar,
           x.max = lon.max.Tsm.scalebar,
           y.min = lat.min.Tsm.scalebar,
           y.max = lat.max.Tsm.scalebar,
           dist = 150,
           dist_unit = 'm',
           location = 'topleft',
           transform = TRUE,
           model = 'WGS84',
           st.size = 5,
           border.size = 0.8,
           box.fill = c("gray60", "white"),
           st.color = "white")

p_map <- ggarrange(p_Mtk, p_Tsm,
                   ncol = 2, nrow = 1, widths = c(1, 1.5), heights = c(1, 1))


################################################################################
##### COMBINE PLOTS #####
################################################################################
p_top <- ggarrange(p_admix, p_pca, ncol = 2, widths = c(0.5, 0.5))
p <- ggarrange(p_top, p_map, ncol = 1, nrow = 2, heights = c(0.5, 0.5))
p <- p +
  draw_plot_label(label = c('A', 'B', 'C'),
                  size = 24,
                  x = c(0, 0.52, 0), y = c(1, 1, 0.48))

## Save plot:
ggsave(p, filename = figfile_eps, width = 11, height = 11)
system(paste0('xdg-open ', figfile_eps))
