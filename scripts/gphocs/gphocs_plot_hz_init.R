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
Log$pop <- factor(Log$pop,levels = poplevels)
Log$migpattern <- factor(Log$migpattern, levels = migpatterns)

## Read metadata - pop colors, and get intermed color for root:
popcols.df <- read.delim(infile_popcols, as.is = TRUE)

## Prep pops:
poplevels <- c('griSW', 'griHZ', 'a.gri', 'murHZ', 'murGan',
               'murW', 'a.murSE', 'a.mur', 'a.root')
migpatterns <- c("griHZ_2_murHZ", "murHZ_2_griHZ", "a.gri_2_a.murSE",
                 "a.murSE_2_a.gri", "a.gri_2_a.mur", "a.mur_2_a.gri")

if(setID == 'hz.mur2gri2c') {
  pops.org <- c('mmur_w', 'mmur_hz', 'mgri_hz', 'mgri_sw', 'anc_mgri', 'anc_mmur', 'anc_root')
  pops <- c('mmur.w', 'mmur.hz', 'mgri.hz', 'mgri.sw', 'anc.mgri', 'anc.mmur', 'anc.root')
  kidpops <- c('murW', 'murHZ', 'griHZ', 'mgriSW', 'a.gri', 'a.mur')
  parentpops <- c('a.mur', 'a.mur', 'a.gri', 'a.gri', 'a.root', 'a.root')
  allpops <- unique(c(kidpops, parentpops))
  currentpops <- kidpops
}
if(setID == 'hz.mur3gri2c') {
  pops.org <- c('mmur_w', 'mmur_hz', 'mmur_gan', 'mgri_hz', 'mgri_sw', 'anc_mmur_se', 'anc_mmur', 'anc_mgri', 'anc_root')
  pops <- c('mmur.w', 'mmur.hz',  'mmur.gan', 'mgri.hz', 'mgri.sw',  'anc.mmur.se','anc.mmur',  'anc.mgri', 'anc.root')
  cbind(pops.org, pops)
  kidpops <- c('murW', 'murHZ', 'murGan', 'griHZ', 'mgriSW', 'a.murSE', 'a.gri', 'a.mur')
  parentpops <- c('a.mur', 'a.murSE', 'a.murSE', 'a.gri', 'a.gri', 'a.mur', 'a.root', 'a.root')
  cbind(kidpops, parentpops)
  allpops <- unique(c(kidpops, parentpops))
  currentpops <- kidpops
}


################################################################################
#### SUMMARIZE #####
################################################################################
(m.sum <- filter(Log, var == '2Nm', runID == 'g2m2anc') %>%
    group_by(migfrom, migto, var) %>%
    summarise(Nm.pt = round(mean(val), 3),
              Nm.min = round(hpd.min(val), 3),
              Nm.max = round(hpd.max(val), 3)))

(m.sum <- filter(Log, var == 'm.prop', runID == 'g2m2anc') %>%
    group_by(migfrom, migto, var) %>%
    summarise(Nm = round(mean(val) * 100, 3)))

(theta.sum <- filter(Log, var == 'theta') %>%
    group_by(pop, var) %>%
    summarise(value = round(mean(cval) / 1000),
              min = round(hpd.min(cval) / 1000),
              max = round(hpd.max(cval) / 1000)) %>%
    arrange(pop)) %>%
  print(n = 100)

(tau.sum <- filter(Log, var == 'tau', runID == 'g2m2anc') %>%
    group_by(pop, var) %>%
    summarise(min = round(hpd.min(cval) / 1000, 1),
              mean = round(mean(cval) / 1000, 1),
              max = round(hpd.max(cval) / 1000, 1)) %>%
    arrange(pop) %>%
    as.data.frame())


################################################################################
#### MIGRATION PLOT #####
################################################################################
runID.focal <- 'g2m2anc'
mlog <- subset(Log, runID == runID.focal)
pop.labs.m <- c("gri-C > mur-C", "mur-C > gri-C",
                "anc.gri > anc.mur-SE", "anc.mur-SE > anc.gri",
                "anc.gri > anc.mur", "anc.mur > anc.gri")

Nm <- vplot(subset(mlog, var == '2Nm'),
            xvar = 'migpattern',
            fillvar = 'cn',
            colvar = 'cn',
            yvar = 'val',
            y.max = 'max.hpd',
            pop.labs = pop.labs.m,
            ylab = '2Nm',
            rotate.x.ann = TRUE,
            yticks.by = 0.01,
            linecols = 'red') +
  theme(plot.margin = margin(1, 1, 0.2, 0.2, 'cm'))
Nm

mprop <- vplot(subset(Log, var == 'm.prop' & runID == runID.focal),
               xvar = 'migpattern',
               fillvar = 'cn',
               colvar = 'cn',
               yvar = 'val',
               y.max = 'max.hpd',
               pop.labs = pop.labs.m,
               ylab = 'migrant percentage',
               rotate.x.ann = TRUE,
               yticks.by = 0.002,
               linecols = 'red') +
  theme(plot.margin = margin(1, 0.2, 0.2, 1, 'cm'))
mprop

p <- ggarrange(Nm, mprop, ncol = 2, nrow = 1) +
  draw_plot_label(label = c("A", "B"),
                  size = 24, x = c(0, 0.5), y = c(1, 1))
figfile <- paste0(plotdir, '/final/', setID, '.', runID.focal, '_allMig.png')
ggexport(p, filename = figfile, width = 800, height = 450)
system(paste('xdg-open', figfile))


################################################################################
#### TAU + THETA PLOT #####
################################################################################
## Tau:
taulog <- subset(Log, var == 'tau' &
                   runID %in% c('noMig', runID.focal) &
                   ! pop %in% c('a.root', 'murW'))
pop.labs.tau <- c('anc.gri', 'anc.mur-SE', 'anc.mur')

tau <- vplot(taulog,
             xvar = 'pop',
             fillvar = 'cn',
             colvar = 'runID',
             yvar = 'cval',
             col.labs = c('G-PhoCS: mig', 'G-PhoCS: iso'),
             pop.labs = pop.labs.tau,
             rotate.x.ann = TRUE,
             yticks.by = 50,
             linecols = NULL,
             legpos = 'top',
             legcolname = "",
             rm.leg.col = FALSE,
             xlab = "",
             ylab = cvar('tau')) +
  theme(plot.margin = margin(0.8, 1.2, 0, 0.2, 'cm'))
tau

## Theta:
thlog <- subset(Log, var == 'theta' & runID %in% c('noMig', runID.focal))
pop.labs.theta <- c('gri-W', 'gri-C', 'anc.gri', 'mur-C', 'mur-E',
                    'mur-W', 'anc.mur-SE', 'anc.mur', 'root')

theta <- vplot(thlog,
              xvar = 'pop',
              fillvar = 'cn',
              colvar = 'runID',
              yvar = 'cval',
              col.labs = c('G-PhoCS: mig', 'G-PhoCS: iso'),
              pop.labs = pop.labs.theta,
              rotate.x.ann = TRUE,
              yticks.by = 50,
              linecols = NULL,
              legpos = 'top',
              legcolname = "",
              rm.leg.col = FALSE,
              xlab = "",
              ylab = cvar('theta')) +
  theme(plot.margin = margin(0.8, 1.2, 0, 0.2, 'cm'))
theta

p <- ggarrange(tau, theta, ncol = 1, nrow = 2) +
  draw_plot_label(label = c("A", "B"),
                  size = 24, x = c(0, 0), y = c(1, 0.5))
figfile <- paste0(plotdir, '/final/', setID, '.', runID.focal, '_tau.theta.png')
ggexport(p, filename = figfile, width = 500, height = 800)
system(paste('xdg-open', figfile))
