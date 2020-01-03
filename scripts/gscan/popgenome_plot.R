## CREATE KARYOPLOT WITH SNP LOCATIONS
## GET PER-SPECIES NR SEG-SITES PER WINDOW? -- should have now
## CORRELATE PI FOR SAME WINDOW ACROSS POPS WITH WIDE FORMAT
## RND requires larger windows?
## CHANGE DENSITY PLOTS TO JOYPLOTS

################################################################################
##### SET-UP #####
################################################################################
setwd('/home/jelmer/Dropbox/sc_lemurs/')
source('scripts/gscan/manplots_fun.R')
library(tidyverse)
library(ggpubr)
library(cowplot)

## FileID and pops:
#fileID <- 'r03.all.mac1.FS6'; winID <- 'win10_step10'
fileID <- 'hzproj1.mac1.FS6'; winID <- 'win25_step5'
do_sitestats <- TRUE

pop1 <- 'griC'
pop2 <- 'griW'
pop3 <- 'murC'
pop4 <- 'murE'
pop5 <- 'murW'
outgroup <- 'ruf'

comps_within <- 'xx' ## TO DO
comps_between <- 'xx'

#scaffolds_to_select <- c(1:6, 33)
scaffolds_to_select <- 1:30

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
#poplabs <- c(expression(italic("griseorufus")), expression(italic("murinus")))

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
         #pi2_gri = (pi_gri  + nucdiv_gri) / 2,
         #pi2_mur = (pi_mur + nucdiv_mur) / 2,
         pi2_griC = (pi_griC + nucdiv_griC) / 2, # CREATE IN RUN SCRIPT!!!
         pi2_griW = (pi_griW  + nucdiv_griW) / 2,
         pi2_murC = (pi_murC + nucdiv_murC) / 2,
         pi2_murE = (pi_murE + nucdiv_murE) / 2,
         pi2_murW = (pi_murW + nucdiv_murW) / 2,
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
    #rename(fst_gri.mur = fst) %>% ## FIX IN RUN-SCRIPT
    mutate(win_index = 1:n(),
           #fstByMaf_gri.mur = fst / MAF_all, ## FIX FST BY MAF!!
           start_run = scaf_site$scafstart_run[match(scaffold, scaf_site$scaffold)] + start,
           end_run = scaf_site$scafstart_run[match(scaffold, scaf_site$scaffold)] + end)
  #head(sitestats_wide)

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
#### PER CHROM MEANS ####
################################################################################
winstats$chromtype <- ifelse(grepl(xchrom, winstats$scaffold), 'sex', 'auto')
cat('#### Win stats in auto vs sex chrom:\n')
winstats %>%
  group_by(chromtype, stat, pop) %>%
  summarize(mean = mean(val, na.rm = TRUE)) %>%
  arrange(stat) %>%
  print(n = Inf)

if(do_sitestats == TRUE) {
  sitestats$chromtype <- ifelse(grepl(xchrom, sitestats$scaffold), 'sex', 'auto')
  cat('#### Site stats in auto vs sex chrom:\n')
  sitestats %>%
    filter(stat != 'sipat') %>%
    group_by(chromtype, stat, pop) %>%
    summarize(mean = mean(val, na.rm = TRUE)) %>%
    arrange(stat)

  # cat('#### Site patterns in auto vs sex chrom:\n')
  # sitestats %>%
  #   filter(stat == 'sipat', pop == 'gri') %>%
  #   group_by(chromtype, val) %>%
  #   tally(val)
}


################################################################################
#### BETWEEN-POP STATS ####
################################################################################
## FST - by window:
fstw_all <- ggman(winstats, yvar = 'fst', xvar = 'win_index', #start_run
                  colvar.lines = 'pop', drawpoints = FALSE,
                  #cols.lines = 'red', cols.points = 'grey50',
                  smoothpar = 0.1, my.ymax = 1,
                  scaf_df = scaf_win)

fstw_sel <- ggman(winstats_scafsel, xvar = 'start_run', yvar = 'fst',
                  colvar.lines = 'pop', drawpoints = FALSE,
                  #cols.lines = 'grey20', cols.points = 'grey50', line.size = 0.5,point.size = 2,
                  smoothpar = 0.1, my.ymax = 1,
                  scaffolds = scafsel, scaf_df = scaf_win)

fst1 <- winstats %>% filter(stat == 'fst' & !grepl('griW', pop))
fst1sel <- winstats_scafsel %>% filter(stat == 'fst' & !grepl('griW', pop))
fst2 <- winstats %>% filter(stat == 'fst' & !grepl('griC', pop))
fst2sel <- winstats_scafsel %>% filter(stat == 'fst' & !grepl('griC', pop))
fst3 <- winstats %>% filter(stat == 'fst' & !grepl('murW', pop))
fst3sel <- winstats_scafsel %>% filter(stat == 'fst' & !grepl('murW', pop))
fst4sel <- winstats_scafsel %>% filter(stat == 'fst' & !grepl('murE', pop))

fstw1 <- ggman(fst1, xvar = 'start_run', yvar = 'fst',
               colvar.lines = 'pop', drawpoints = FALSE,
               smoothpar = 0.1, my.ymax = 0.6, scaf_df = scaf_win)
fstw1sel <- ggman(fst1sel, xvar = 'start_run', yvar = 'fst',
                  colvar.lines = 'pop', cols.lines = 'popcols',
                  drawpoints = FALSE, smoothpar = 0.1, my.ymax = 0.55,
                  scaf_df = scaf_win)
fstw2 <- ggman(fst2, xvar = 'start_run', yvar = 'fst',
               colvar.lines = 'pop', drawpoints = FALSE,
               smoothpar = 0.1, my.ymax = 0.55, scaf_df = scaf_win)
fstw2sel <- ggman(fst2sel, xvar = 'start_run', yvar = 'fst',
               colvar.lines = 'pop', cols.lines = 'popcols',
               drawpoints = FALSE,
               smoothpar = 0.1, my.ymax = 0.55, scaf_df = scaf_win)
fstw3 <- ggman(fst3, xvar = 'start_run', yvar = 'fst',
               colvar.lines = 'pop', drawpoints = FALSE,
               smoothpar = 0.1, my.ymax = 0.55, scaf_df = scaf_win)
fstw3sel <- ggman(fst3sel, xvar = 'start_run', yvar = 'fst',
               colvar.lines = 'pop', cols.lines = 'popcols',
               drawpoints = FALSE,
               smoothpar = 0.1, my.ymax = 0.6, scaf_df = scaf_win)
fstw4sel <- ggman(fst4sel, xvar = 'start_run', yvar = 'fst',
                  colvar.lines = 'pop', cols.lines = 'popcols',
                  drawpoints = FALSE,
                  smoothpar = 0.1, my.ymax = 0.6, scaf_df = scaf_win)

## FST - by site:
fsts_all <- ggman(sitestats, yvar = 'fst',
                  colvar.lines = 'pop', drawpoints = FALSE,
                  #cols.lines = 'red', cols.points = 'grey50', my.ymax = 1,
                  smoothpar = 0.1, my.ymin = 0.1, my.ymax = 0.5,
                  ylab = expression(paste('per-site ', F[ST])),
                  scaf_df = scaf_site)
fsts_sel <- ggman(sitestats_scafsel, yvar = 'fst',
                  colvar.lines = 'pop', drawpoints = FALSE,
                  #cols.lines = 'grey20', cols.points = 'grey50', line.size = 0.5, my.ymax = 1,
                  smoothpar = 0.1, my.ymin = 0.1, my.ymax = 0.5,
                  ylab = expression(paste('per-site ', F[ST])),
                  scaf_df = scaf_site)

## FST/MAF - by site:
highmaf1 <- filter(sitestats, stat == 'MAF', pop == pop1, val > 0.2) %>% pull(start_run)
highmaf2 <- filter(sitestats, stat == 'MAF', pop == pop2, val > 0.2) %>% pull(start_run)
highmaf3 <- filter(sitestats, stat == 'fstByMaf', val < 1) %>% pull(start_run)
highmaf <- intersect(intersect(highmaf1, highmaf2), highmaf3)

fm_all <- ggman(sitestats_MAFsel, yvar = 'fstByMaf',
                cols.lines = 'red', cols.points = 'grey50',
                drawpoints2 = TRUE, points2 = highmaf, cols.points2 = 'red',
                smoothpar = 0.1, my.ymax = 2.8,
                ylab = expression(paste(F[ST], ' by MAF')),
                scaf_df = scaf_site)

fm_sel <- ggman(sitestats_scafsel, yvar = 'fstByMaf',
                cols.lines = 'grey20', cols.points = 'grey50',
                drawpoints2 = TRUE, points2 = highmaf, cols.points2 = 'red',
                smoothpar = 0.1, my.ymax = 2.8, line.size = 0.5,
                ylab = expression(paste(F[ST], ' by MAF')),
                scaf_df = scaf_site)

fm_sel <- ggman(sitestats_scafsel_MAFsel, yvar = 'fstByMaf',
                cols.lines = 'grey20', cols.points = 'grey50',
                drawpoints2 = TRUE, points2 = highmaf, cols.points2 = 'red',
                smoothpar = 0.1, my.ymax = 2.8, line.size = 0.5,
                ylab = expression(paste(F[ST], ' by MAF')),
                scaf_df = scaf_site)

## DXY - by window:
dxy_all <- ggman(winstats, yvar = 'dxy',
                 #cols.lines = 'red', cols.points = 'grey50', my.ymax = 5,
                 colvar.lines = 'pop', drawpoints = FALSE, my.ymax = 2,
                 smoothpar = 0.1, scaf_df = scaf_win)

dxy_sel <- ggman(winstats_scafsel, yvar = 'dxy',
                 #cols.lines = 'grey20', cols.points = 'grey50', line.size = 0.5,
                 colvar.lines = 'pop', drawpoints = FALSE, my.ymax = 2,
                 point.size = 2, smoothpar = 0.1, scaf_df = scaf_win)

dxy1 <- winstats %>% filter(stat == 'dxy' & !grepl('griW', pop))
dxy1sel <- winstats_scafsel %>% filter(stat == 'dxy' & !grepl('griW', pop))
dxy2 <- winstats %>% filter(stat == 'dxy' & !grepl('griC', pop))
dxy2sel <- winstats_scafsel %>% filter(stat == 'dxy' & !grepl('griC', pop))
dxy3 <- winstats %>% filter(stat == 'dxy' & !grepl('murW', pop))
dxy3sel <- winstats_scafsel %>% filter(stat == 'dxy' & !grepl('griC', pop))
dxy4sel <- winstats_scafsel %>% filter(stat == 'dxy' & !grepl('murE', pop))

dxyw1 <- ggman(dxy1, xvar = 'start_run', yvar = 'dxy',
               colvar.lines = 'pop', drawpoints = FALSE,
               smoothpar = 0.1, scaf_df = scaf_win,
               my.ymin = 1.3, my.ymax = 3.7)
dxyw1sel <- ggman(dxy1sel, xvar = 'start_run', yvar = 'dxy',
               colvar.lines = 'pop', drawpoints = FALSE,
               smoothpar = 0.1, scaf_df = scaf_win,
               my.ymin = 1.3, my.ymax = 3.7)
dxyw2 <- ggman(dxy2, xvar = 'start_run', yvar = 'dxy',
               colvar.lines = 'pop', drawpoints = FALSE,
               smoothpar = 0.1, scaf_df = scaf_win,
               my.ymin = 1.3, my.ymax = 3.7)
dxyw3 <- ggman(dxy3, xvar = 'start_run', yvar = 'dxy',
               colvar.lines = 'pop', drawpoints = FALSE,
               smoothpar = 0.1, scaf_df = scaf_win,
               my.ymin = 1.3, my.ymax = 3.7)
dxyw4sel <- ggman(dxy4sel, xvar = 'start_run', yvar = 'dxy',
                  colvar.lines = 'pop', drawpoints = FALSE,
                  smoothpar = 0.1, scaf_df = scaf_win,
                  my.ymin = 1.3, my.ymax = 3.7)


################################################################################
#### WITHIN-POP STATS ####
################################################################################
## MAF - by site:
maf_all <- ggman(filter(sitestats, pop != 'all'), yvar = 'MAF',
                 colvar.lines = 'pop', cols.lines = popcols, smoothpar = 0.1,
                 colvar.points = 'pop', cols.points = popcols,
                 my.ymin = 0, my.ymax = 0.2, legplot = FALSE,
                 scaf_df = scaf_site)
maf_sel <- ggman(filter(sitestats_scafsel, pop != 'all'), yvar = 'MAF',
                 colvar.lines = 'pop', cols.lines = popcols, smoothpar = 0.2,
                 colvar.points = 'pop', cols.points = popcols,
                 my.ymin = 0, my.ymax = 0.25, legplot = FALSE,
                 scaf_df = scaf_site)

## Pi:
pi_all <- ggman(winstats, yvar = 'pi',
                colvar.lines = 'pop', cols.lines = popcols,
                #colvar.points = 'pop', cols.points = popcols, my.ymin = 0, my.ymax = 2,
                drawpoints = FALSE, my.ymin = 0.8, my.ymax = 2,
                smoothpar = 0.15, legplot = FALSE, scaf_df = scaf_win)

pi_sel <- ggman(winstats_scafsel, yvar = 'pi',
                colvar.lines = 'pop', cols.lines = popcols,
                #colvar.points = 'pop', cols.points = popcols, my.ymin = 0, my.ymax = 3,
                drawpoints = FALSE, my.ymin = 0.8, my.ymax = 2,
                smoothpar = 0.15, legplot = FALSE, scaf_df = scaf_win)

## Nucdiv:
nucdiv_all <- ggman(winstats, yvar = 'nucdiv',
                    colvar.lines = 'pop', cols.lines = popcols,
                    drawpoints = FALSE, my.ymin = 1, my.ymax = 2.2,
                    smoothpar = 0.15, legplot = FALSE, scaf_df = scaf_win)
nucdiv_sel <- ggman(winstats_scafsel, yvar = 'nucdiv',
                    colvar.lines = 'pop', cols.lines = popcols,
                    drawpoints = FALSE, my.ymin = 1, my.ymax = 2.2,
                    smoothpar = 0.15, legplot = FALSE, scaf_df = scaf_win)

## pi2
pi2_sel <- ggman(winstats_scafsel, yvar = 'pi2',
                 colvar.lines = 'pop', cols.lines = popcols,
                 drawpoints = FALSE, my.ymin = 0.9, my.ymax = 1.9,
                 smoothpar = 0.15, legplot = FALSE, scaf_df = scaf_win)

## Tajima's D:
taj_all <- ggman(winstats, yvar = 'tajD', #xvar = 'start_run',
                 colvar.lines = 'pop', cols.lines = popcols,
                 #colvar.points = 'pop', cols.points = popcols, my.ymin = -2.1, my.ymax = 2.1,
                 drawpoints = FALSE, my.ymin = -0.6, my.ymax = 1,
                 smoothpar = 0.15, hline = 0, legpos = 'bottom', scaf_df = scaf_win)
taj_sel <- ggman(winstats_scafsel, yvar = 'tajD',
                 colvar.lines = 'pop', cols.lines = popcols,
                 #colvar.points = 'pop', cols.points = popcols, my.ymin = -2.2, my.ymax = 2.2,
                 drawpoints = FALSE, my.ymin = -0.5, my.ymax = 1,
                 smoothpar = 0.2, hline = 0, legpos = 'bottom', scaf_df = scaf_win)

## Fu & Li's F:
fu_all <- ggman(winstats, yvar = 'fuLiF',
                colvar.lines = 'pop', cols.lines = popcols,
                #colvar.points = 'pop', cols.points = popcols, my.ymin = -2.1, my.ymax = 2.1,
                drawpoints = FALSE, my.ymin = -1, my.ymax = 0.8,
                smoothpar = 0.2, hline = 0, legpos = 'bottom', scaf_df = scaf_win)
fu_sel <- ggman(winstats_scafsel, yvar = 'fuLiF',
                colvar.lines = 'pop', cols.lines = popcols,
                #colvar.points = 'pop', cols.points = popcols, my.ymin = -2.2, my.ymax = 2.2,
                drawpoints = FALSE, my.ymin = -1, my.ymax = 0.8,
                smoothpar = 0.2, hline = 0, legpos = 'bottom', scaf_df = scaf_win)


################################################################################
#### INTROGRESSION STATS ####
################################################################################
## fd:
fdstats <- winstats %>% filter(stat == 'fd')
table(fdstats$pop)

fd_all <- ggman(winstats, yvar = 'fd', drawpoints = FALSE,
                colvar.lines = 'pop', #cols.lines = popcols,
                my.ymin = -0.12, my.ymax = 0.12,
                legpos = 'top', legtextsize = 10,
                smoothpar = 0.2, hline = 0, scaf_df = scaf_win)
fd_sel <- ggman(winstats_scafsel, yvar = 'fd', drawpoints = FALSE,
                colvar.lines = 'pop', #cols.lines = popcols,
                my.ymin = -0.12, my.ymax = 0.12,
                legpos = 'top', legtextsize = 10,
                smoothpar = 0.2, hline = 0, scaf_df = scaf_win)

## df:
df_all <- ggman(winstats, yvar = 'df',
                colvar.lines = 'pop', drawpoints = FALSE,
                legpos = 'top', legtextsize = 10,
                my.ymin = -0.1, my.ymax = 0.2,
                smoothpar = 0.2, hline = 0, scaf_df = scaf_win)
df_sel <- ggman(winstats_scafsel, yvar = 'df',
                colvar.lines = 'pop', drawpoints = FALSE,
                legpos = 'top', legtextsize = 10,
                my.ymin = -0.1, my.ymax = 0.2,
                smoothpar = 0.3, hline = 0, scaf_df = scaf_win)

## RNDmin:
RND_all <- ggman(winstats, yvar = 'RNDmin',
                 colvar.lines = 'pop', drawpoints = FALSE,
                 my.ymin = 0, my.ymax = 0.15,
                 smoothpar = 0.2, hline = 0, scaf_df = scaf_win)
RND_sel <- ggman(winstats_scafsel, yvar = 'RNDmin',
                 colvar.lines = 'pop', drawpoints = FALSE,
                 my.ymin = 0, my.ymax = 0.15,
                 smoothpar = 0.2, hline = 0,
                 scaf_df = scaf_win)


################################################################################
#### ARRANGE MANHATTAN PLOTS FOR FINAL FIGURE ####
################################################################################
# pman <- ggarrange(fstw_sel, fsts_sel, fm_sel, dxy_sel,
#                   maf_sel, pi_sel, taj_sel, ncol = 1)

pman <- ggarrange(fstw4sel, pi2_sel, taj_sel, ncol = 1)
figfile <- paste0(figfilebase, 'man_scafsel.eps')
ggsave(figfile, width = 8, height = 12)
system(paste0('xdg-open ', figfile))

pman <- ggarrange(fstw1sel, fstw2sel, fstw3sel, ncol = 1)
figfile <- paste0(figfilebase, 'man_scafsel2_2.eps')
ggsave(figfile, width = 10, height = 8)
system(paste0('xdg-open ', figfile))


################################################################################
#### PLOTS - DENSITY ####
################################################################################
pdens_taj <- densplot(winstats,
                      mystat = 'tajD',
                      fill_var = 'pop',
                      fill_name = 'species',
                      fill_vals = popcols,
                      #fill_labs = poplabs,
                      fill_legend = FALSE,
                      xlims = c(-3, 3), ylims = c(0, 0.5))

pdens_pi <- densplot(winstats,
                     mystat = 'pi',
                     fill_var = 'pop',
                     fill_name = 'species',
                     fill_vals = popcols,
                     #fill_labs = poplabs,
                     fill_legend = FALSE,
                     xlims = c(0, 3), ylims = c(0, 0.9))

pdens_nuc <- densplot(winstats,
                      mystat = 'nucdiv',
                      fill_var = 'pop',
                      fill_name = 'species',
                      fill_vals = popcols,
                      #fill_labs = poplabs,
                      fill_legend = FALSE,
                      xlab = 'nucleotide diversity')

pdens_MAF <- densplot(filter(sitestats, pop != 'all'),
                      mystat = 'MAF',
                      fill_var = 'pop',
                      fill_name = 'species',
                      fill_vals = c(popcols, 'black'),
                      fill_labs = poplabs,
                      ylims = c(0, 25),
                      xlab = 'MAF')

pdens <- ggarrange(pdens_MAF, pdens_nuc, pdens_pi, pdens_taj,
                   ncol = 1, nrow = 4, heights = c(1, 1, 1, 1))
figfile <- paste0(figfilebase, 'densPlots.png')
ggsave(figfile, width = 8, height = 10)
system(paste0('xdg-open ', figfile))


################################################################################
#### FST AND MAF ####
################################################################################
p_fst_maf <- ggplot(data = sitestats_wide) +
  geom_point(aes(x = MAF_all, y = fst)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = 'MAF', y = expression(F[ST])) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.margin = margin(0.4, 0.4, 0.1, 0.1, 'cm'))
figfile <- paste0(figfilebase, 'FstVsMAF.png')
ggsave(figfile, width = 6, height = 5)
system(paste0('xdg-open ', figfile))

ggplot(data = sitestats_wide) +
  geom_point(aes(x = MAF_gri, y = fst))

ggplot(data = sitestats_wide) +
  geom_point(aes(x = MAF_mur, y = fst))

##
p_maf_all <- ggplot(data = sitestats_wide) +
  geom_point(aes(x = MAF_all, y = fstByMaf_gri.mur)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
  labs(x = 'MAF', y = expression(paste(F[ST], ' by MAF'))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.margin = margin(0.4, 0.4, 0.1, 0.1, 'cm'))
p_maf_all

p_maf_gri <- ggplot(data = sitestats_wide) +
  geom_point(aes(x = MAF_gri, y = fstByMaf_gri.mur)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
  labs(x = 'griseorufus MAF', y = expression(paste(F[ST], ' by MAF'))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.margin = margin(0.4, 0.4, 0.1, 0.1, 'cm'))
p_maf_gri

p_maf_mur <- ggplot(data = sitestats_wide) +
  geom_point(aes(x = MAF_mur, y = fstByMaf_gri.mur)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
  labs(x = 'murinus MAF', y = expression(paste(F[ST], ' by MAF'))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.margin = margin(0.4, 0.4, 0.1, 0.1, 'cm'))
p_maf_mur

## Arrange and save:
p_maf <- ggarrange(p_maf_all, p_maf_gri, p_maf_mur,
               ncol = 1, nrow = 3, heights = c(1, 1, 1))
figfile <- paste0(figfilebase, 'FstByMaf.png')
ggsave(figfile, width = 8, height = 8)
system(paste0('xdg-open ', figfile))
#draw_plot_label(label = c('A', 'B'), size = 25, x = c(0, 0), y = c(1, 0.5))


################################################################################
#### OTHER PLOTS ####
################################################################################
## Number of windows by scaffold length:
# scaf_win <- scaf_all %>%
#   filter(scaffold %in% unique(winstats$scaffold)) %>%
#   arrange(scaffold) %>%
#   mutate(nwin = as.numeric(table(winstats$scaffold)))
#
# plot(scaf_win$len, scaf_win$nwin)

################################################################################
#### PLOTS - KARYPLOTE-R ####
################################################################################
# library(karyoploteR)
# win_gr <- toGRanges(winstats_wide)
# names(mcols(win_gr))
# #seqlevelsStyle(win_gr) <- "UCSC"
# head(win_gr)
#
# kp <- plotKaryotype(plot.type = 4)
# kpPoints(kp, data = win_gr, y = win_gr$fst_gri.mur,
#          cex = 1, r0 = 0.5, r1 = 1)
#kpPoints(kp, data=snp.data, y=snp.data$LRR, cex=0.3, r0=0, r1=0.5, ymax=2, ymin=-3)
