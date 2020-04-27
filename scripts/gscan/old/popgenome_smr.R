################################################################################
##### SET-UP #####
################################################################################
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')
library(tidyverse)
source('/home/jelmer/Dropbox/scripts/genomics/gscan/manplots_fun.R')

## FileID:
fileID <- 'r03.all.mac1.FS6'
fileID_win <- 'r03.all.mac1.FS6_win10_step5'

## Files:
indir_pg <- 'analyses/gscan/popgenome/output/'
infile_scafs <- '/home/jelmer/Dropbox/sc_lemurs/other/seqdata_misc/reference/mmur/scaffolds_withLength.txt'
infile_winstats <- paste0(indir_pg, '/', fileID_win, '_winstats.txt')

## Read data:
winstats <- read.delim(infile_winstats, as.is = TRUE)
head(winstats)


################################################################################
#### PLOTS ####
################################################################################
for(fpop in unique(fst$pop)) {
  #fpop <- unique(fst$pop)[1]
  ggman(fst, yvar = 'fst', scaffold = 'all', pop = fpop,
        drawline = FALSE, drawpoints = TRUE,
        my.ymin = 0, my.ymax = 0.8, plot.title = fpop,
        saveplot = FALSE)
}

yticks <- seq(0, 0.05, by = 0.05/10)
ggman(fst, yvar = 'fst', colvar.lines = 'pop', scaffold = 'all',
      drawline = TRUE, cols.lines = c('red', 'blue', 'green'), drawpoints = FALSE,
      my.ymin = 0, my.ymax = 0.05, my.yticks = yticks, xlab = 'window index',
      plot.title = FALSE, saveplot = FALSE)



## Fish plots:
manplot(pi, var = 'pi', my.ylim = c(0, 0.3), plot.id = 'pi', win.pos = win.pos)
manplot(fst.sl, 'single-pop Fst', my.ylim = c(-1, 1), plot.id = 'fst.singlepop', win.pos = win.pos)

manplot(fst.pr, 'between-pop Fst', plot.id = 'fst', my.ylim = c(-1, 1), win.pos = win.pos)
manplot(dxy, var = 'dxy', plot.id = 'dxy', my.ylim = c(0, 0.2), win.pos = win.pos)

manplot(BDF, plot.id = 'BDF', var = 'BDF', win.pos = win.pos)
manplot(introgstats, plot.id = 'introgressionStats',
          legnames = colnames(introgstats), var = 'introgression', win.pos = win.pos)

manplot(pi[, c('Cdec', 'Ceja', 'Cfus')], var = 'pi', plot.title = 'pi.eja', plot.id = 'pi.eja', my.ylim = c(0, 0.3), win.pos = win.pos)
manplot(fst.pr[, mam.columns], 'Fst', plot.title = 'fst.mam', plot.id = 'fst.mam', my.ylim = c(-1, 1), win.pos = win.pos)
manplot(dxy[, mam.columns], var = 'dxy', plot.title = 'dxy.mam', plot.id = 'dxy.mam', my.ylim = c(0, 0.2), win.pos = win.pos)

plot(win.pos, D, pch = 19, ylab = "D", xlab = "genomic position", main = triplet.name, ylim = c(-1, 1))

## SFS:
# ggplot(sfs_pop1, aes(x = MAF, y = freq, width = 0.02)) +
#   #geom_bar(stat = 'identity', position = 'identity') +
#   geom_density(stat = 'identity', position = 'identity') +
#   scale_x_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   theme_bw()
