################################################################################
#### SET-UP #####
################################################################################
## Libraries and scripts:
library(data.table)
library(png)
library(grid)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(TeachingDemos)
library(plyr)
library(reshape2)
library(tidyverse)
library(here)

script_plotfun_general <- here('scripts/generalscripts_lemurs_link/gphocs/gphocs_plot_fun.R')
script_plotfun_hz <- here('scripts/gphocs/gphocs_plot_fun_hz.R')
source(script_plotfun_general)
source(script_plotfun_hz)

## Settings:
setID <- 'hz.mur3gri2c'
gentime <- 3.75
mutrate.gen <- 1.64e-8
mutrate.year <- mutrate.gen / gentime
m.scale <- 1000
t.scale <- 0.0001

## Input files:
infile_logs <- paste0(here(), '/analyses/gphocs/output/', setID, '/mergedLogs.txt')
infile_cols <- here('analyses/gphocs/popinfo/ghocs_cols.txt')
infile_inds.df <- '/home/jelmer/Dropbox/sc_lemurs/radseq/metadata/lookup_IDshort.txt'

## Output files:
plotdir <- '/home/jelmer/Dropbox/sc_lemurs/presentations/posters/2019_03_GordonConf/plots/'
figfile_png <- file.path(plotdir, 'gphocs_demo.png')
figfile_eps <- file.path(plotdir, 'gphocs_demo.eps')

## Read log file:
Log <- as.data.frame(fread(infile_logs, stringsAsFactors = TRUE))
Log <- Log %>% dplyr::filter(runID == 'noMig' | runID == 'g2m2anc')

## Read metadata - pop colors, and get intermed color for root:
popcols.df <- read.delim(infile_popcols, as.is = TRUE)
col.gri.a <- popcols.df$popcol[popcols.df$popname.short == 'a.gri']
col.mur.a <- popcols.df$popcol[popcols.df$popname.short == 'a.mur']
midval <- get.midpoint(col.gri.a, col.mur.a)
popcols.df$popcol[popcols.df$popname.short == 'a.root'] <- midval

inds.df <- read.delim(infile_inds.df, as.is = TRUE)

## Pops:
pops.org <- c('mmur_w', 'mmur_hz', 'mmur_gan', 'mgri_hz', 'mgri_sw',
              'anc_mmur_se', 'anc_mmur', 'anc_mgri', 'anc_root')
pops <- c('mmur.w', 'mmur.hz',  'mmur.gan', 'mgri.hz', 'mgri.sw',
          'anc.mmur.se','anc.mmur',  'anc.mgri', 'anc.root')
cbind(pops.org, pops)
kidpops <- c('murW', 'murHZ', 'murGan', 'griHZ', 'griSW', 'a.murSE',
             'a.gri', 'a.mur')
parentpops <- c('a.mur', 'a.murSE', 'a.murSE', 'a.gri', 'a.gri', 'a.mur',
                'a.root', 'a.root')
cbind(kidpops, parentpops)
allpops <- unique(c(kidpops, parentpops))
currentpops <- c('murW', 'murHZ', 'murGan', 'griHZ', 'griSW')
ancpops <- c('a.root', 'a.gri', 'a.mur', 'a.murSE')


################################################################################
#### CREATE PLOT #####
################################################################################
## Prepare df underlying plot:
runID.focal <- 'g2m2anc'
ttp <- ttPrepMur3(subset(Log, runID == runID.focal),
                  x.even = FALSE, pop.spacing = 75, extra.time.root = 25)
ttp$pop <- factor(ttp$pop, levels = levels(Log$pop))
ttp <- arrange(ttp, pop)

## Plot:
p <- dplot(tt = ttp, ann.pops = FALSE, x.min = 0, yticks.by = 100,
           popnames.adj.horz = rep(0, nrow(ttp)), popnames.adj.vert = 15,
           popnames.col = popnames.col, popnames.size = 5, x.extra = 5,
           saveplot = FALSE, plot.title = '') +
  geom_segment(aes(y = ttp$y.max[ttp$pop == 'a.mur'], yend = ttp$y.max[ttp$pop == 'a.mur'],
                   x = ttp$x.max[ttp$pop == 'a.mur'], xend = ttp$x.min[ttp$pop == 'a.gri']),
               colour = 'grey50') +
  theme(axis.text.y = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = unit(c(15, 0, 0, 0), 'mm')),
        axis.title.y = element_text(size = 32, margin = unit(c(0, 5, 0, 0), 'mm')),
        panel.border = element_rect(colour = 'grey20', size = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

ggsave(filename = outfile_pcplot.eps, width = 9, height = 10)
system(paste("xdg-open", outfile_pcplot.eps))


################################################################################
## Migration summary:
(m.sum <- filter(Log, var == '2Nm', runID == runID.focal) %>%
    group_by(migfrom, migto, var) %>%
    summarise(Nm_min = round(hpd.min(val), 4),
              Nm_mean = round(mean(val), 4),
              Nm_max = round(hpd.max(val), 4)))
