################################################################################
##### SET-UP #####
################################################################################
setwd('/home/jelmer/Dropbox/sc_lemurs/proj/hybridzone/')
source('/home/jelmer/Dropbox/scripts/genomics/admixture/admixture_plot_fun.R')

## SetID:
setID <- 'r03.all.mac3.FS6'

## Files:
infile_inds <- '/home/jelmer/Dropbox/sc_lemurs/radseq/metadata/lookup_IDshort.txt'
infile_cols <- '/home/jelmer/Dropbox/sc_lemurs/metadata/colors/colors.species.txt'
indir <- 'analyses/admixture/output/'
figdir <- 'analyses/admixture/figures/'

## Get metadata:
#cols.sp.df <- read.delim(infile_cols, as.is = TRUE)
inds.df <- read.delim(infile_inds, as.is = TRUE)

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
##### PLOT #####
################################################################################
## CV-plot:
kplot <- k.plot(setID)

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
         file.save = TRUE,
         figfile = paste0(figdir, setID, '_K2.eps'))

## K=3:
k3 <- Qdf(setID, K = 3, sort.by = 'species.short2') %>%
  ggax.v(.,
         group.column = c('species.short2', 'site.type'),
         grouplab.bgcol = my.grouplab.bgcol,
         grouplab.labeller = grouplab.labeller,
         file.save = FALSE,
         figfile = paste0(figdir, setID, '_K3.png'))


################################################################################
##### COMBINE PLOTS #####
################################################################################
## Combine plots:
p <- ggarrange(kplot, k2, k3,
               nrow = 3, heights = c(1.7, 1, 1))
#plots <- plots + draw_plot_label(label = c('A', 'B', 'C'), size = 24,
#                                 x = c(0, 0, 0.4), y = c(1, 0.55, 0.9))

## Save plot:
figfile <- paste0('analyses/admixture/figures/', fileID, '.combined.eps')
ggexport(p, filename = figfile, width = 8, height = 8)
system(paste0('xdg-open ', figfile))
