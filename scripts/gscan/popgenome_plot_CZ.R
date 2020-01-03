################################################################################
##### SET-UP #####
################################################################################
setwd('/home/jelmer/Dropbox/sc_lemurs/')
source('scripts/gscan/manplots_fun.R')
library(tidyverse)
library(ggpubr)
library(cowplot)

## FileID and pops:
fileID <- 'r03.all.mac1.FS6'; winID <- 'win10_step10'
do_sitestats <- TRUE

pop1 <- 'gri'
pop2 <- 'mur'
outgroup <- NA

scaffolds_to_select <- c(1:6, 33)

## Process:
fileID_win <- paste0(fileID, '_', winID)
pops <- c(pop1, pop2, pop3, pop4, pop5) %>% .[!is.na(.)] %>% sort(.)
npop <- length(pops)

## Files:
indir_pg <- 'hybridzone/analyses/gscan/popgenome/output/'
infile_winstats <- paste0(indir_pg, '/', fileID_win, '_winstats.txt')
infile_sitestats <- paste0(indir_pg, '/', fileID, '_sitestats.txt')

infile_scafs <- 'other/seqdata_misc/reference/mmur/scaffolds_withLength2.txt'
infile_cols <- 'metadata/colors/popcols.txt'

outdir_plots <- 'hybridzone/analyses/gscan/popgenome/figures/'
figfilebase <- paste0(outdir_plots, fileID_win, '_')

## Metadata - colours:
cat("##### Reading cols.file ...\n")
cols_df <- read.delim(infile_cols, header = TRUE, as.is = TRUE)
popcols <- cols_df$col[match(pops, cols_df$pop_short)]
names(popcols) <- pops
popcols

## Metadata - scaffolds:
xchrom <- 'NC_033692.1'
scaf_all <- read.table(infile_scafs, header = TRUE,
                       colClasses = c('character', 'integer', 'integer'))


################################################################################
##### READ & PREP DATA #####
################################################################################
## Read & prep windowstats:
winstats_wide <- read.delim(infile_winstats, as.is = TRUE)

scaf_win <- filter(scaf_all, scaffold %in% unique(winstats_wide$scaffold))
scaf_win$scafstart_run <- c(1, cumsum(as.numeric(scaf_win$length)) + 1)[1:nrow(scaf_win)]
winstats_wide$scaf_index <- scaf_win$scaf_index[match(winstats_wide$scaffold, scaf_win$scaffold)]
winstats_wide <- winstats_wide %>%
  arrange(scaf_index) %>%
  mutate(win_index = 1:n(),
         pi2_gri = (pi_gri  + nucdiv_gri) / 2,
         pi2_mur = (pi_mur + nucdiv_mur) / 2,
         start_run = scaf_win$scafstart_run[match(scaffold, scaf_win$scaffold)] + start,
         end_run = scaf_win$scafstart_run[match(scaffold, scaf_win$scaffold)] + end)

winstats <- winstats_wide %>%
  gather("stat", "val", -scaffold, -start, -end,
         -start_run, -end_run, -scaf_index, -win_index) %>%
  separate(stat, into = c('stat', 'pop'), sep = '_') %>%
  mutate(scaf_index = match(scaffold, scaf_all$scaffold))
cat('#### Win stats:\n')
head(winstats, n = 2)
table(winstats$stat)
table(winstats$pop)

## Data subset with selection of scaffolds:
scafsel <- filter(scaf_all, scaf_index %in% scaffolds_to_select) %>% pull(scaffold)
winstats_scafsel <- filter(winstats, scaffold %in% scafsel)

## Read & prep sitestats:
if(do_sitestats == TRUE) {
  sitestats_wide <- read.delim(infile_sitestats, as.is = TRUE)
  scaf_site <- filter(scaf_all, scaffold %in% unique(sitestats_wide$scaffold))
  scaf_site$scafstart_run <- c(1, cumsum(as.numeric(scaf_site$length)) + 1)[1:nrow(scaf_site)]
  sitestats_wide$scaf_index <- scaf_site$scaf_index[match(sitestats_wide$scaffold, scaf_site$scaffold)]
  sitestats_wide <- sitestats_wide %>%
    arrange(scaf_index) %>%
    rename(fst_gri.mur = fst) %>% ## FIX IN RUN-SCRIPT
    mutate(
      win_index = 1:n(),
      fstByMaf_gri.mur = fst_gri.mur / MAF_all, ## FIX FST BY MAF!!
      start_run = scaf_site$scafstart_run[match(scaffold, scaf_site$scaffold)] + start,
      end_run = scaf_site$scafstart_run[match(scaffold, scaf_site$scaffold)] + end
      )

  sitestats <- sitestats_wide %>%
    gather("stat", "val", -scaffold, -start, -end,
           -start_run, -end_run, -win_index, -scaf_index) %>%
    separate(stat, into = c('stat', 'pop'), sep = '_') %>%
    mutate(scaf_index = match(scaffold, scaf_all$scaffold))
  cat('#### Site stats:\n')
  print(head(sitestats, n = 2))
  print(table(sitestats$stat))
  print(table(sitestats$pop))

  ## Data subset with selection of scaffolds:
  sitestats_scafsel <- filter(sitestats, scaffold %in% scafsel)

  ## Data with MAF filter:
  highmaf_sites <- filter(sitestats, stat == 'MAF', pop == pop1, val > 0.05) %>% pull(start)
  length(highmaf_sites)
  sitestats_MAFsel <- filter(sitestats, start %in% highmaf_sites)
  sitestats_scafsel_MAFsel <- filter(sitestats_scafsel, start %in% highmaf_sites)
}


################################################################################
##### MANPLOTS #####
################################################################################
#poplabs <- c(expression(italic("griseorufus")), expression(italic("murinus")))

highmaf1 <- filter(sitestats, stat == 'MAF', pop == pop1, val > 0.2) %>% pull(start_run)
highmaf2 <- filter(sitestats, stat == 'MAF', pop == pop2, val > 0.2) %>% pull(start_run)
highmaf3 <- filter(sitestats, stat == 'fstByMaf', val < 1) %>% pull(start_run)
highmaf <- intersect(intersect(highmaf1, highmaf2), highmaf3)


fstw_sel <- ggman(winstats_scafsel, xvar = 'win_index', yvar = 'fst',
                  cols.lines = 'red', cols.points = 'grey50', drawpoints = TRUE,
                  line.size = 1, point.size = 2,
                  smoothpar = 0.15, my.ymax = 1,
                  scaf_df = scaf_win, scaffolds = scafsel)

fsts_sel <- ggman(sitestats_scafsel, xvar = 'win_index', yvar = 'fst',
                  drawpoints = TRUE, cols.points = 'grey50', cols.lines = 'red',
                  smoothpar = 0.1, my.ymin = 0.0, my.ymax = 1,
                  ylab = expression(paste('per-site ', F[ST])),
                  scaf_df = scaf_site, scaffolds = scafsel)

fm_sel <- ggman(sitestats_scafsel_MAFsel, xvar = 'win_index', yvar = 'fstByMaf',
                cols.lines = 'red', cols.points = 'grey50',
                drawpoints2 = TRUE, points2 = highmaf, cols.points2 = 'red',
                smoothpar = 0.1, my.ymax = 2.8,
                ylab = expression(paste(F[ST], ' by MAF')),
                scaf_df = scaf_site, scaffolds = scafsel)

dxy_sel <- ggman(winstats_scafsel, xvar = 'win_index', yvar = 'dxy',
                 cols.lines = 'red', cols.points = 'grey50',
                 drawpoints = TRUE, my.ymax = 4,
                 point.size = 2, smoothpar = 0.15,
                 scaf_df = scaf_site, scaffolds = scafsel)

pi_sel <- ggman(winstats_scafsel, xvar = 'win_index', yvar = 'pi',
                colvar.lines = 'pop', cols.lines = popcols, drawpoints = FALSE,
                my.ymin = 0, my.ymax = 2, smoothpar = 0.15,
                legplot = FALSE,
                scaf_df = scaf_site, scaffolds = scafsel)

taj_sel <- ggman(winstats_scafsel, xvar = 'win_index', yvar = 'tajD',
                 colvar.lines = 'pop', cols.lines = popcols, drawpoints = FALSE,
                 my.ymin = -1, my.ymax = 0.6, smoothpar = 0.3,
                 hline = 0, legpos = 'bottom',
                 scaf_df = scaf_site, scaffolds = scafsel)


pman <- ggarrange(fstw_sel, fsts_sel, fm_sel, dxy_sel,
                  pi_sel, taj_sel, ncol = 1)
figfile <- paste0(figfilebase, 'man_gri-mur_scafsel2.eps')
ggsave(figfile, width = 10, height = 8)
system(paste0('xdg-open ', figfile))
