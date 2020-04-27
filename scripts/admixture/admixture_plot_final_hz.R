################################################################################
##### SET-UP #####
################################################################################
## Libraries and scripts:
library(here)
source('/home/jelmer/Dropbox/scripts/genomics/admixture/admixture_plot_fun.R')

## SetID:
setID <- 'r03.all.mac3.FS6'

## Input files:
infile_inds <- here('metadata/radseq_metadata_link/lookup_IDshort.txt')
infile_cols <- here('metadata/lemurs_metadata_link/lookup/popcols.txt')
indir <- here('analyses/admixture/output/')

## Output files:
figdir <- here('analyses/admixture/figs_final/')

## Read metadata:
cols.sp.df <- read.delim(infile_cols, as.is = TRUE) %>%
  select(-species)
inds.df <- read.delim(infile_inds, as.is = TRUE) %>%
  merge(., cols.sp.df, by = 'species.short') %>%
  rename(labcol = color)

## Sympatric vs parapatric sites:
inds.df$site.type <- gsub('Mangatsiaka', 'sympatric', inds.df$site)
inds.df$site.type <- gsub('Tsimelahy', 'sympatric', inds.df$site.type)
inds.df$site.type <- gsub('Hazofotsy', 'parapatric', inds.df$site.type)
inds.df$site.type <- gsub('Ambatoabo', 'parapatric', inds.df$site.type)
inds.df$site.type[inds.df$site.type == 'parapatric' &
                    inds.df$species.short == 'mmur'] <- 'zparapatric'
inds.df$site.type[grep('patric', inds.df$site.type, invert = TRUE)] <- NA
#table(inds.df$site.type)

## Shorter species abbrev:
inds.df$species.short2 <- gsub('^m', '', inds.df$species.short)
inds.df$species.short2 <- gsub('hyb', 'hyb?', inds.df$species.short2)
#table(inds.df$species.short2)

## Colours:
col.mur <- '#FF00FF'
col.gri <- '#D69B12'


################################################################################
##### K2 PLOT #####
################################################################################
## CV-plot:
k.plot(setID)

## Plot prep:
my.grouplab.bgcol <- c(col.gri, col.gri, 'grey50', col.mur, col.mur)
site.type.labs <- c(parapatric = 'parapatric', sympatric = 'sympatric',
                    zparapatric = 'parapatric')
grouplab.labeller <- labeller(site.type = site.type.labs)

## K=2:
k2 <- Qdf(setID, K = 2, sort.by = 'species.short2') %>%
  ggax.v(.,
         group.column = c('species.short2', 'site.type'),
         barcols = c(col.gri, col.mur),
         grouplab.bgcol = my.grouplab.bgcol,
         grouplab.labeller = grouplab.labeller,
         ylab = 'RADseq')

mt <- Qdf(setID, K = 2, sort.by = 'species.short2') %>%
  ggax.v(.,
         group.column = c('species.short2', 'site.type'),
         barcols = c(col.gri, col.mur),
         ylab = 'mtDNA') +
  theme(strip.background = element_blank(),
        strip.text = element_blank())
mt

p <- ggarrange(k2, mt, ncol = 1, nrow = 2, heights = c(1, 0.2))
figfile <- paste0(figdir, '/', setID, '_K2_combined.eps')
ggsave(p, filename = figfile, width = 9, height = 5)
system(paste0('xdg-open ', figfile))


################################################################################
##### K3 PLOT #####
################################################################################
k2 <- Qdf(setID, K = 2, sort.by = 'species.short2') %>%
  ggax.v(.,
         group.column = c('species.short2', 'site.type'),
         barcols = c(col.gri, col.mur),
         grouplab.bgcol = my.grouplab.bgcol,
         grouplab.labeller = grouplab.labeller,
         ylab = 'RADseq: K=2')

(k3 <- Qdf(setID, K = 3, sort.by = 'species.short2') %>%
  ggax.v(.,
         group.column = c('species.short2', 'site.type'),
         barcols = c(col.mur, col.gri, 'grey40'),
         ylab = 'RADseq: K=3') +
  theme(strip.background = element_blank(),
        strip.text = element_blank()))

(mt <- Qdf(setID, K = 2, sort.by = 'species.short2') %>%
  ggax.v(.,
         group.column = c('species.short2', 'site.type'),
         barcols = c(col.gri, col.mur),
         ylab = 'mtDNA') +
  theme(strip.background = element_blank(),
        strip.text = element_blank()))

p <- ggarrange(k2, k3, mt, ncol = 1, nrow = 3, heights = c(1, 0.8, 0.3))
figfile <- paste0(figdir, '/', setID, '_K3_combined.eps')
ggsave(p, filename = figfile, width = 8, height = 6)
system(paste0('xdg-open ', figfile))


################################################################################
##### OLD PLOTTING #####
################################################################################
#plots <- plots + draw_plot_label(label = c('A', 'B', 'C'), size = 24,
#                                 x = c(0, 0, 0.4), y = c(1, 0.55, 0.9))

# plot1k <- function(ind.set, my.k = 2, indlab.column = 'ID.short', ...) {
#   admix.ggplot(ind.set = ind.set, K = my.k,
#                ID.type = 'ID.short', indlab.colum = indlab.column,
#                col.by = 'species.short', col.labs = TRUE,
#                plotwidth = 6, plotheight = 8, ...)
# }
# notany <- sapply(ind.sets, plot1k)
# notany <- mapply(plot1k, ind.set, my.k = 1:10,
#                  indlab.column = 'supersite2', indlab.size = 8)
#
# plot1k(ind.set, my.k = 2, indlab.column = 'supersite2', indlab.size = 8)
