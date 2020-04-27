################################################################################
##### SET-UP #####
################################################################################
library(tidyverse)
source('/home/jelmer/Dropbox/scripts/genomics/PCA/PCA_R_fun.R')
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')

## Vars:
file.ID <- 'r03.all.mac3.FS6'
ID.type <- 'ID.short'
pca.df.dir <- 'analyses/PCA/dfs/'

## Files:
infile_lookup <- '/home/jelmer/Dropbox/sc_lemurs/radseq/metadata/lookup_IDshort.txt'
infile_pca.df <- paste0('analyses/PCA/dfs/', file.ID, '_all.txt')
infile_eigenvalues <- paste0('analyses/PCA/dfs/', file.ID, '_all_eigenvalues.txt')

plotdir <- '/home/jelmer/Dropbox/sc_lemurs/presentations/posters/2019_03_GordonConf/plots/'
outfile_pcplot.eps <- paste0(plotdir, '/PCA.eps')

## Metadata:
lookup <- read.delim(infile_lookup, as.is = TRUE) %>%
  select(-ID)
col.mur <- '#FF00FF'
col.gri <- '#D69B12'
col.hyb <- 'black'

## Get eigenvalues:
eigenvalues <- read.delim(infile_eigenvalues, as.is = TRUE)
var.PC1 <- round(eigenvalues$eig[1] / sum(eigenvalues$eig) * 100, 1)
var.PC2 <- round(eigenvalues$eig[2] / sum(eigenvalues$eig) * 100, 1)

## Prep df:
pca.df <- read.table(infile_pca.df, header = TRUE, as.is = TRUE) %>%
  select(ID, PC1, PC2, PC3, PC4) %>%
  merge(., lookup, by.x = 'ID', by.y = 'ID.short') %>%
  select(ID, PC1, PC2, PC3, PC4, site, species.short, species.cor, supersite) %>%
  mutate(supersite = gsub('mmur.hz', 'contact zone', supersite)) %>%
  mutate(supersite = gsub('mgri.hz', 'contact zone', supersite)) %>%
  mutate(supersite = gsub('mhyb.hz', 'contact zone', supersite)) %>%
  mutate(supersite = gsub('mmur.se', 'nearby allopatric', supersite)) %>%
  mutate(supersite = gsub('mgri.se', 'nearby allopatric', supersite)) %>%
  mutate(species.short = factor(species.short, levels = c('mgri', 'mmur', 'mhyb'))) %>%
  arrange(species.short)
#head(pca.df)

## Labels:
sp.cols <- c(col.gri, col.mur, col.hyb)
sp.labs <- c(expression(italic("griseorufus")), expression(italic("murinus")), 'hybrid')
msat.name <- 'previous\nmicrosatellite\nassignment:'
site.labs <- c('contact zone', 'nearby allopatric')


################################################################################
##### PLOT PCA #####
################################################################################
p <- ggplot(data = pca.df) +
  geom_point(aes(x = PC1, y = PC2, color = species.short, shape = supersite),
             size = 6.5, stroke = 3.5) +
  labs(x = paste0('PC1 (', var.PC1, '%)'), y = paste0('PC2 (', var.PC2, '%)')) +
  scale_color_manual(values = sp.cols, labels = sp.labs, name = msat.name) +
  scale_shape_manual(values = c(1, 0), name = '\nsite type:', labels = site.labs) +
  theme_bw() +
  theme(legend.text = element_text(size = 28),
        legend.title = element_text(size = 28, face = 'bold'),
        legend.spacing.x = unit(0, 'cm'),
        legend.spacing.y = unit(0.3, 'cm'),
        legend.key.size = unit(1.3, "cm"),
        legend.text.align = 0,
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        panel.border = element_rect(colour = 'grey20', size = 1),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 30, margin = unit(c(4, 0, 0, 0), 'mm')),
        axis.title.y = element_text(size = 30, margin = unit(c(0, 4, 0, 0), 'mm')),
        plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm")) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))
p

#ggsave(outfile_pcplot.eps, width = 11.5, height = 11.5*(5/7))
#system(paste('xdg-open', outfile_pcplot.eps))
