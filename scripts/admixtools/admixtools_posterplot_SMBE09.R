################################################################################
##### SET-UP #####
################################################################################
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')
source('/home/jelmer/Dropbox/sc_lemurs/scripts/admixtools/admixtools_plot_fun.R')
library(tidyverse)
library(ggpubr)
library(cowplot)

## File ID:
#file.id.short <- 'r03.wOut'
#file_id <- 'r03.wOutgroups.indsel.mac1.FS6.dstat_r03.wOut'
file.id.short <- file_id <- 'hzproj1.mac1.FS6.bySupersite2'
plotdir <- '/home/jelmer/Dropbox/sc_lemurs/other/posters/2019_07_SMBE/'
figfile_eps <- paste0(plotdir, '/Dstats.eps')

popcombs <- c('(mur-W,mur-C),gri-C', '(mur-W,mur-C),gri-W','(mur-W,mur-E),gri-W',
              '(mur-E,mur-C),gri-W', '(gri-W,gri-C),mur-C')

## Prep df with d-stats:
(d.df <- return.dfmode(file.id = file_id) %>%
    filter(popcomb %in% popcombs))

d.df$sigcol <- cut(abs(d.df$Z), breaks = c(0, 3, Inf),
                labels = c('black', 'red'))
d.df$siglab <- gsub('black', '|Z|<3', d.df$sig)
d.df$siglab <- gsub('red', '|Z|>3', d.df$siglab)

col.labs.all <- c('|Z|<3', '|Z|<2', '3>|Z|>2', '|Z|>3')
col.labs <- col.labs.all[col.labs.all %in% d.df$siglab]


################################################################################
##### PLOT ####
################################################################################
p <- ggplot(d.df, aes(x = factor(popcomb), y = -D)) +
  geom_pointrange(aes(ymax = -d.df$D + d.df$se,
                      ymin = -d.df$D - d.df$se,
                      colour = sigcol), size = 1.6) +
  geom_hline(yintercept = 0, colour = 'grey40') +
  coord_flip() +
  scale_colour_identity(guide = 'legend', labels = col.labs, name = '') +
  scale_y_continuous(breaks = c(0, 0.05, 0.1)) +
  labs(y = "D value\n[D>0: p2-p3 admixture]") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 22),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 26, margin = margin(12, 0, 0, 0)),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(1, 1, 1.2, 0.8, "cm"),
    panel.border = element_rect(colour = 'grey20', size = 1),
    legend.position = 'top',
    legend.text = element_text(size = 22))
p

ggsave(figfile_eps, p, width = 5, height = 9)
system(paste('xdg-open', figfile.eps))
