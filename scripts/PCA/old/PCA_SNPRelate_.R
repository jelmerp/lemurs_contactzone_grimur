################################################################################
##### SET-UP #####
################################################################################
## Variables:
fileID <- 'r03.all.mac3.FS6' # VCF should be [dir_vcf]/[fileID].vcf.gz
ID.type <- 'ID.short' # "ID" (atri001_r01_xxxx) or "ID.short" (atri001)

## Other scripts:
source('/home/jelmer/Dropbox/sc_lemurs/scripts/PCA/PCA_R_fun.R')

## Files and dirs:
dir_base <- '/home/jelmer/Dropbox/sc_lemurs/hybridzone/'
dir_vcf <- 'seqdata/vcf/'
dir_pca_df <- 'analyses/PCA/dfs/'

infile_vcf <- paste0(dir_vcf, '/', fileID, '.vcf.gz')
infile_lookup <- '/home/jelmer/Dropbox/sc_lemurs/radseq/metadata/lookup_IDshort.txt'
#infile_cols <- '/home/jelmer/Dropbox/sc_lemurs/metadata/colors/popcols.txt' # Could use to assign colours automatically

## Base dir and libraries:
setwd(dir_base)
library(gdsfmt)
library(SNPRelate)
library(vcfR)
library(tidyverse)
library(ggpubr)
library(cowplot)

## Species colours:
#cols.df <- read.delim(infile_cols, header = TRUE, as.is = TRUE)
col.mur <- '#FF00FF'
col.gri <- '#D69B12'
col.hyb <- 'black'

## Read metadata:
cat("##### Reading lookup file...\n")
lookup <- read.delim(infile_lookup, as.is = TRUE)


################################################################################
##### READ VCF #####
################################################################################
outfile_snpgds <- paste0('analyses/PCA/gdsFiles/', fileID, '.gds')
snpgdsVCF2GDS(infile_vcf, outfile_snpgds, method = "biallelic.only")
snps <- snpgdsOpen(outfile_snpgds)

snpgdsSummary(outfile_snpgds)


################################################################################
##### SUBSETTING #####
################################################################################
## LD pruning:
#snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
#snpset.id <- unlist(unname(snpset))
#pca <- snpgdsPCA(snpgds.file, snp.id = snpset.id, num.thread = 2)

## List of individuals:
inds_all <- read.gdsn(index.gdsn(snps, "sample.id"))


################################################################################
##### ALL INDS #####
################################################################################
subset.ID = 'all'
outfile_pca.df <- paste0(dir_pca_df, '/', fileID, '_', subset.ID, '.txt')
outfile_eigenvals <- paste0(dir_pca_df, '/', fileID, '_', subset.ID,
                            '_eigenvals.txt')

## Run PCA:
pca_all_raw <- snpgdsPCA(snps, autosome.only = FALSE, num.thread = 1)

## Process results:
pca_all <- pca.process(pca_all_raw,
                      inds.df = lookup,
                      source = 'snpgdsPCA',
                      subset.ID = subset.ID,
                      pca.df.ID = paste0(fileID, '_', subset.ID)) %>%
  mutate(supersite = gsub('mmur.hz', 'sympatric', supersite)) %>%
  mutate(supersite = gsub('mgri.hz', 'sympatric', supersite)) %>%
  mutate(supersite = gsub('mhyb.hz', 'sympatric', supersite)) %>%
  mutate(supersite = gsub('mmur.se', 'parapatric', supersite)) %>%
  mutate(supersite = gsub('mgri.se', 'parapatric', supersite)) %>%
  mutate(species.short = factor(species.short,
                                levels = c('mgri', 'mmur', 'mhyb'))) %>%
  arrange(species.short)

## Eigenvalues:
eig_all <- data.frame(PC = 1:length(pca_all_raw$eigenval),
                        eig = pca_all_raw$eigenval)
eig_all <- eig_all[complete.cases(eig_all), ]
write.table(eig_all, outfile_eigenvals,
            sep = '\t', quote = FALSE, row.names = FALSE)

## Plot prep:
col.by.labs <- c(expression(italic("griseorufus")),
                 expression(italic("murinus")), 'hybrid')
shape.by.labs <- c('contact zone', 'nearby allopatric')

## Plot:
p1 <- pcplot(pca.df,
             eigenvalues = eigenvals$eig,
             col.by = 'species.short',
             col.by.name =  'microsats:',
             col.by.labs = col.by.labs,
             cols = c(col.gri, col.mur, col.hyb),
             shape.by = 'site.type',
             shape.by.name = 'site type:',
             shape.by.labs = shape.by.labs,
             shapes = c(1, 0))


################################################################################
##### PCA - MGRI ONLY #####
################################################################################
## Subset ID:
subset.ID <- 'mgri_cz'

## Select inds:
hybs_mgri <- c('mhyb001', 'mhyb003', 'mhyb004', 'mhyb005', 'mhyb006',
               'mhyb011', 'mhyb015', 'mhyb016')
inds_mgri <- c(inds_all[grepl('mgri', inds_all)],
              inds_all[inds_all %in% hybs_mgri])

## Run PCA:
pca_mgri_raw <- snpgdsPCA(snps, sample.id = inds_mgri,
                      autosome.only = FALSE, num.thread = 1)

## Process results:
pca_mgri <- pca.process(pca_mgri_raw,
                        inds.df = lookup,
                        subset.ID = subset.ID,
                        pca.df.ID = paste0(fileID, '_', subset.ID))

## Eigenvalues:
eig_gri <- data.frame(PC = 1:length(pca_mgri_raw$eigenval),
                        eig = pca_mgri_raw$eigenval)
eig_gri <- eig_gri[complete.cases(eig_gri), ]
write.table(eig_gri, outfile_eigenvals,
            sep = '\t', quote = FALSE, row.names = FALSE)

## Plot prep:
labs_species <- c(expression(italic("griseorufus")), 'hybrid')
title_gri <- expression(paste(italic('griseorufus'), ' only'))

## Plot:
p_mgri <- pcplot(pca_mgri,
                 eigenvalues = eig_gri$eig,
                 col.by = 'site',
                 col.by.name = 'site:',
                 shape.by = 'species.short',
                 shape.by.name = 'microsats:',
                 shape.by.labs = labs_species,
                 shapes = c(1, 0),
                 plot.title = title_gri) +
  theme(plot.margin = margin(0.2, 0.2, 0.6, 0.2, "cm"))


################################################################################
##### PCA- MMUR ONLY #####
################################################################################
## Subset ID:
subset.ID <- 'mmur_cz'

## Select inds:
hybs_mmur <- c('mhyb002', 'mhyb007', 'mhyb008', 'mhyb009',
               'mhyb010', 'mhyb012', 'mhyb013', 'mhyb014')
inds_mmur <- c(inds_all[grepl('mmur', inds_all)],
               inds_all[inds_all %in% hybs_mmur])

## Run PCA:
pca_mmur_raw <- snpgdsPCA(snps, sample.id = inds_mmur,
                          autosome.only = FALSE, num.thread = 1)

## Process results:
pca_mmur <- pca.process(pca_mmur_raw,
                        inds.df = lookup,
                        subset.ID = subset.ID,
                        pca.df.ID = paste0(fileID, '_', subset.ID)) %>%
  mutate(species.short = factor(species.short, levels = c('mmur', 'mhyb')))

## Eigenvalues:
eig_mur <- data.frame(PC = 1:length(pca_mmur_raw$eigenval),
                        eig = pca_mmur_raw$eigenval)
eig_mur <- eig_mur[complete.cases(eig_mur), ]
write.table(eig_mur, outfile_eigenvals,
            sep = '\t', quote = FALSE, row.names = FALSE)

## Plot prep:
labs_species <- c(expression('hybrid', italic("murinus")))
title_mur <- expression(paste(italic('murinus'), ' only'))

## Plot:
p_mmur <- pcplot(pca_mmur,
                 eigenvalues = eig_mur$eig,
                 col.by = 'site',
                 col.by.name = 'site:',
                 shape.by = 'species.short',
                 shape.by.name = 'microsats:',
                 shape.by.labs = labs_species,
                 shapes = c(1, 0),
                 plot.title = title_mur) +
  theme(plot.margin = margin(0.6, 0.2, 0.2, 0.2, "cm"))


################################################################################
##### COMBINE PLOTS #####
################################################################################
## Arrange  panels:
p <- ggarrange(p_mgri, p_mmur,
               ncol = 1, nrow = 2, heights = c(1, 1)) +
  draw_plot_label(label = c('A', 'B'), size = 25,
                  x = c(0, 0), y = c(1, 0.5))

## Save:
figfile_eps <- paste0('analyses/PCA/figures/', fileID, '_singlespecies.eps')
ggsave(filename = figfile_eps, width = 6, height = 10)
system(paste0('xdg-open ', figfile_eps))
