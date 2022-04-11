################################################################################
##### SET-UP #####
################################################################################
## Libraries and scripts:
library(tidyverse)
library(ggpubr)
library(cowplot)
library(here)

script_atools_fun <- here('scripts/admixtools/admixtools_plot_fun.R')
source(script_atools_fun)

## Output dirs and files:
plotdir <- 'results/admixtools/figures/'
figfile_png <- paste0(plotdir, '/admixstats.png')
figfile_eps <- paste0(plotdir, '/admixstats.eps')


################################################################################
##### D-STATS BY SUPERSITE2 ####
################################################################################
file.id.short <- file.id.d <- 'hzproj1.mac1.FS6.bySupersite2'

## Prep df with d-stats:
(d.df <- return.dfmode(file.id = file.id.d) %>%
    mutate(popcomb = gsub('mgri', 'gri', popcomb)) %>%
    mutate(popcomb = gsub('mmur', 'mur', popcomb)) %>%
    mutate(popcomb = gsub('hybridzone', 'hz', popcomb)) %>%
    mutate(popcomb = gsub('_w', '-W', popcomb)) %>%
    mutate(popcomb = gsub('_se', '-SE', popcomb)) %>%
    mutate(popcomb = gsub('_hz', '-C', popcomb)) %>%
    mutate(popcomb = gsub('_gan', '-E', popcomb)) %>%
    mutate(popcomb = gsub('_sw', '-W', popcomb)))
d.df$sig <- gsub('b_nearsig', 'a_nonsig', d.df$sig)

## Plot:
(dplot <- plot.dstats(d.df, fig.save = FALSE) +
  labs(y = "D") +
  theme(axis.text.x = element_text(size = 18),
       axis.title.x = element_text(size = 20),
       plot.margin = margin(0.2, 0.5, 0.1, 0.1, "cm")))


################################################################################
##### F4-RATIO BY SUPERSITE2 ####
################################################################################
## Prep df:
file.id.f4r <- 'hzproj1.mac1.FS6.bySupersite2'
(f.df <- return.f4ratio(file.id = file.id.f4r) %>% arrange(popA))

f.df1 <- f.df[7:8, ]
f.df2 <- f.df[1:6, ]

## Plot:
(fplot1 <- plot.f4r(f.df1, ylims = c(0.94, 1.02), alpha.breaks = c(0.95, 1)) +
  guides(colour = FALSE) +
    theme(axis.title.x = element_blank(),
          plot.margin = margin(t = 0.8, r = 0.5, b = 0.1, l = 0.5, "cm")))

(fplot2 <- plot.f4r(f.df2, ylims = c(-0.01, 0.07), alpha.breaks = c(0, 0.05)) +
    guides(colour = FALSE) +
    labs(y = expression(paste(alpha))) +
    theme(axis.title.x = element_text(size = 20),
          plot.margin = margin(t = 0.1, r = 0.5, b = 0.1, l = 0.5, "cm")))

(fplot <- ggarrange(fplot1, fplot2, nrow = 2, heights = c(0.25, 0.75)))


################################################################################
##### COMBINE PLOTS ####
################################################################################
dfstats <- ggarrange(dplot, fplot, ncol = 2, widths = c(0.5, 0.5))
dfstats <- dfstats +
  draw_plot_label(label = c('A', 'B'), size = 24, x = c(0, 0.5), y = c(1, 1))

## Save plot:
ggsave(dfstats, filename = figfile_eps, width = 12, height = 8)
system(paste0('xdg-open ', figfile_eps))
