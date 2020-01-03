################################################################################
##### SET-UP #####
################################################################################
## Libraries and scripts:
library(tidyverse)
library(ggpubr)
library(cowplot)
library(here)

script_atools_fun <- here('scripts/generalscripts_lemurs_link/admixtools_plot_fun.R')
source(script_atools_fun)

## File ID:
## output should be: analyses/admixtools/output/$file.id.dmode.out
file.id.short <- 'r03.wOut'
file.id.d <- 'r03.wOutgroups.indsel.mac1.FS6.dstat_r03.wOut'
file.id.f4r <- 'r03.wOutgroups.indsel.mac1.FS6.f4ratio_r03.wOut'


################################################################################
##### D-STATS BY SPECIES ####
################################################################################
#figfile <- paste0('analyses/admixtools/figures/', file.id.d, '.png')
figfile <- paste0('analyses/admixtools/figures/', file.id.d, '.poster.png')

popcombs <- c('murW-murHZ-griHZ', 'murW-murSE-griSW','murW-murGan-griSW',
              'murGan-murHZ-griSW',
              'griSW-griHZ-murHZ')

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

## Plot:
dplot.sp <- plot.dstats(d.df, figfile = figfile)


################################################################################
##### POSTER PLOT ####
################################################################################
figfile <- paste0('analyses/admixtools/figures/', file.id.d, '.poster.png')

d.df$sig <- gsub('b_nearsig', 'a_nonsig', d.df$sig)

p <- ggplot(d.df, aes(x = factor(popcomb), y = -D)) +
  geom_pointrange(aes(ymax = -d.df$D + d.df$se, ymin = -d.df$D - d.df$se,
                             colour = sig), size = 2) +
  geom_hline(yintercept = 0, colour = 'grey40') +
  coord_flip() +
  scale_colour_manual(values = c('black', 'red'),
                      labels = c('|Z| < 3', '|Z| > 3'), name = '') +
  scale_y_continuous(breaks = c(0, 0.05, 0.1)) +
  labs(y = "D\n[D>0: p2p3 admixture]") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 22,
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 18, hjust = 2),
        plot.margin = margin(1, 1, 1.2, 0.8, "cm"),
        legend.position = 'top',
        legend.text = element_text(size = 18))
p
ggsave(figfile, p, width = 7, height = 7)
system(paste('xdg-open', figfile))


################################################################################
##### F4-RATIO ####
################################################################################
return.f4ratio(file.id = file.id.f4r) %>%
  arrange(popA)


################################################################################
##### COMBINE PLOTS ####
################################################################################
# ## Arrange A and B panels:
# plots <- ggarrange(dplot.sp, dplot.pop,
#                    ncol = 2, nrow = 1, widths = c(1, 1.2))
# plots <- plots + draw_plot_label(label = c("A", "B"), size = 24,
#                                  x = c(0, 0.45), y = c(1, 1))
#
# ## Save as png:
# figfile <- paste0('analyses/admixtools/figures/Dstats_', file.id.short, '.png')
# ggexport(plots, filename = figfile, width = 1000, height = 650)
# system(paste0('xdg-open ', figfile))
