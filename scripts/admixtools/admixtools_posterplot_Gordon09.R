################################################################################
##### SET-UP #####
################################################################################
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')
source('/home/jelmer/Dropbox/scripts/genomics/admixtools/admixtools_plot_fun.R')
library(tidyverse)
library(ggpubr)
library(cowplot)

## File ID:
## output should be: analyses/admixtools/output/$file.id.dmode.out
file.id.short <- 'r03.wOut'
file.id.d <- 'r03.wOutgroups.indsel.mac1.FS6.dstat_r03.wOut'
file.id.f4r <- 'r03.wOutgroups.indsel.mac1.FS6.f4ratio_r03.wOut'
plotdir <- '/home/jelmer/Dropbox/sc_lemurs/presentations/posters/2019_03_GordonConf/plots/'
figfile.eps <- paste0(plotdir, '/Dstats.eps')

popcombs <- c('murW-murHZ-griHZ', 'murW-murSE-griSW','murW-murGan-griSW',
              'murGan-murHZ-griSW', 'griSW-griHZ-murHZ')

## Prep df with d-stats:
(d.df <- return.dfmode(file.id = file.id.d) %>%
    mutate(popcomb = gsub('^m', '', popcomb)) %>%
    mutate(popcomb = gsub('--m', '-', popcomb)) %>%
    mutate(popcomb = gsub('hybridzone', 'hz', popcomb)) %>%
    mutate(popcomb = gsub('_w', 'W', popcomb)) %>%
    mutate(popcomb = gsub('_se', 'SE', popcomb)) %>%
    mutate(popcomb = gsub('_hz', 'HZ', popcomb)) %>%
    mutate(popcomb = gsub('_gan', 'Gan', popcomb)) %>%
    mutate(popcomb = gsub('_sw', 'SW', popcomb)) %>%
    filter(popcomb %in% popcombs))
d.df$sig <- gsub('b_nearsig', 'a_nonsig', d.df$sig)

## Plot:
dplot.sp <- plot.dstats(d.df, figfile = figfile)

################################################################################
##### PLOT ####
################################################################################
p <- ggplot(d.df, aes(x = factor(popcomb), y = -D)) +
  geom_pointrange(aes(ymax = -d.df$D + d.df$se,
                      ymin = -d.df$D - d.df$se,
                      colour = sig), size = 1.6) +
  geom_hline(yintercept = 0, colour = 'grey40') +
  coord_flip() +
  scale_colour_manual(values = c('black', 'red'),
                      labels = c('|Z| < 3', '|Z| > 3'),
                      name = '') +
  scale_y_continuous(breaks = c(0, 0.05, 0.1)) +
  labs(y = "D value\n[D>0: p2-p3 admixture]") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 26,
                                    margin = margin(t = 12, r = 0, b = 0, l = 0)),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = margin(1, 1, 1.2, 0.8, "cm"),
        panel.border = element_rect(colour = 'grey20', size = 1),
        legend.position = 'top',
        legend.text = element_text(size = 22))
p

ggsave(figfile.eps, p, width = 5, height = 9)
system(paste('xdg-open', figfile.eps))


################################################################################
##### F4-RATIO ####
################################################################################
#return.f4ratio(file.id = file.id.f4r) %>% arrange(popA)
