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

## Output files:
plotdir <- here('analyses/gphocs/figs_final/')
figfile_png <- file.path(plotdir, 'gphocs.png')
figfile_eps <- file.path(plotdir, 'gphocs.eps')

## Read log file:
Log <- as.data.frame(fread(infile_logs, stringsAsFactors = TRUE))

## Read metadata - pop colors, and get intermed color for root:
popcols.df <- read.delim(infile_popcols, as.is = TRUE)
col.gri.a <- popcols.df$popcol[popcols.df$popname.short == 'a.gri']
col.mur.a <- popcols.df$popcol[popcols.df$popname.short == 'a.mur']
midval <- get.midpoint(col.gri.a, col.mur.a)
popcols.df$popcol[popcols.df$popname.short == 'a.root'] <- midval

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
#### DEMOGRAPHY PLOT #####
################################################################################
## Prepare df underlying plot:
runID.focal <- 'g2m2anc'
ttp <- ttPrepMur3(subset(Log, runID == runID.focal),
                  x.even = FALSE, pop.spacing = 75, extra.time.root = 25)
ttp$pop <- factor(ttp$pop, levels = levels(Log$pop))
(ttp <- arrange(ttp, pop))

## Plot:
(d <- dplot(tt = ttp,
            ann.pops = FALSE,
            x.min = 0,
            yticks.by = 100,
            x.extra = 5,
            saveplot = FALSE) +
  geom_segment(aes(y = ttp$y.max[ttp$pop == 'a.mur'],
                   yend = ttp$y.max[ttp$pop == 'a.mur'],
                   x = ttp$x.max[ttp$pop == 'a.mur'],
                   xend = ttp$x.min[ttp$pop == 'a.gri']),
               colour = 'grey50') +
  theme(axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 22, margin = unit(c(10, 0, 0, 0), 'mm')),
        axis.title.y = element_text(size = 22, margin = unit(c(0, 2, 0, 0), 'mm')),
        panel.border = element_rect(colour = 'grey20', size = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = margin(0.5, 0.3, 0.1, 0.1, 'cm')))


################################################################################
#### MIGRATION PLOT #####
################################################################################
poplevels <- c('murW', 'murHZ', 'murGan', 'a.mur', 'a.murSE',
               'griHZ', 'griSW','a.gri', 'a.root')
miglevels <- c("griHZ_2_murHZ", "murHZ_2_griHZ",
               "a.gri_2_a.murSE", "a.murSE_2_a.gri",
               "a.gri_2_a.mur", "a.mur_2_a.gri")
Log$pop <- factor(Log$pop, levels = poplevels)
Log$migpattern <- factor(Log$migpattern, miglevels)

runID.focal <- 'g2m2anc'
mlog <- subset(Log, runID == runID.focal)
pop.labs.m <- c("gri-C > mur-C", "mur-C > gri-C",
                "C > B", "B > C",
                "C > A", "A > C")

(Nm <- vplot(subset(mlog, var == '2Nm'),
            xvar = 'migpattern',
            fillvar = 'cn',
            colvar = 'cn',
            yvar = 'val',
            pop.labs = pop.labs.m,
            ylab = '2Nm',
            rotate.x.ann = TRUE,
            yticks.by = 0.02,
            linecols = 'red',
            saveplot = FALSE) +
    geom_vline(xintercept = 2.5, color = "grey30", size = 1) +
    geom_vline(xintercept = 4.5, color = "grey30", size = 1) +
    theme(axis.text.x = element_text(size = 16),
          axis.title.x = element_blank(),
          plot.margin = margin(0.5, 0.2, 0.1, 1, 'cm')))


################################################################################
#### TAU PLOT #####
################################################################################
taulog <- subset(Log,
                 var == 'tau' &
                   runID %in% c('noMig', runID.focal) &
                   ! pop %in% c('a.root', 'murW'))
#pop.labs.tau <- c('gri-W — gri-C', 'mur-C — mur-E', 'mur-W — B')
pop.labs.tau <- c('A', 'B', 'C')

(tau <- vplot(taulog,
              xvar = 'pop',
              fillvar = 'cn',
              colvar = 'runID',
              yvar = 'cval',
              col.labs = c('G-PhoCS: mig', 'G-PhoCS: iso'),
              pop.labs = pop.labs.tau,
              rotate.x.ann = FALSE,
              yticks.by = 50,
              linecols = NULL,
              legpos = 'none',
              legcolname = "",
              rm.leg.col = FALSE,
              xlab = "",
              ylab = cvar('tau'),
              saveplot = FALSE) +
    labs(y = 'div. time (ka ago)') +
    theme(axis.text.x = element_text(size = 18),
          axis.title.x = element_blank(),
          plot.margin = margin(0.4, 0.2, 0.1, 1, 'cm')))


################################################################################
#### THETA PLOT #####
################################################################################
thlog <- subset(Log, var == 'theta' & runID %in% c('noMig', runID.focal))
pop.labs.theta <- c('mur-W', 'mur-C', 'mur-E', 'A', 'B',
                    'gri-C', 'gri-W','C', 'root')

(theta <- vplot(thlog,
                xvar = 'pop',
                fillvar = 'cn',
                colvar = 'runID',
                yvar = 'cval',
                col.labs = c('migration', 'isolation'),
                pop.labs = pop.labs.theta,
                rotate.x.ann = TRUE,
                yticks.by = 50,
                linecols = NULL,
                legpos = 'bottom',
                legcolname = "model",
                rm.leg.col = FALSE,
                xlab = "",
                ylab = cvar('theta'),
                saveplot = FALSE) +
    geom_vline(xintercept = 5.5, color = "grey30", size = 1) +
    geom_vline(xintercept = 8.5, color = "grey30", size = 1) +
    theme(axis.text.x = element_text(size = 16),
          axis.title.x = element_blank(),
          legend.box.margin = margin(0, 0, 0, 0),
          plot.margin = margin(0.4, 0.2, 0.1, 1, 'cm')))


################################################################################
#### COMBINE PLOTS #####
################################################################################
rightside <- ggarrange(Nm, tau, theta,
                       ncol = 1, nrow = 3,
                       heights = c(1.3, 1, 1.2)) +
  draw_plot_label(label = c('B', 'C', 'D'),
                  size = 24,
                  x = c(0, 0, 0), y = c(1, 0.67, 0.38))
p <- ggarrange(d, rightside,
               ncol = 2, widths = c(0.5, 0.5)) +
  draw_plot_label(label = 'A', size = 24, x = 0, y = 1)

ggsave(p, filename = figfile_eps, width = 12, height = 10)
system(paste('xdg-open', figfile_eps))

#ggexport(rightside, filename = figfile_png, width = 500, height = 800)
#ggexport(p, filename = figfile_png, width = 800, height = 800)
#system(paste('xdg-open', figfile_png))
