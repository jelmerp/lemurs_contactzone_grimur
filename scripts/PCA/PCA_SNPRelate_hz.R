################################################################################
##### SET-UP #####
################################################################################
## Variables:
basedir <- '/home/jelmer/Dropbox/sc_lemurs/proj/hybridzone/' # Base directory
fileID <- 'r03.all.mac3.FS6' # VCF name
ID.type <- 'ID.short' # "ID" (mmur001_r01_p1b01) or "ID.short" (mmur001)

## Other scripts:
source('/home/jelmer/Dropbox/scripts/genomics/PCA/PCA_R_fun.R')

## Files and dirs:
setwd(basedir)
infile_lookup <- '/home/jelmer/Dropbox/sc_lemurs/radseq/metadata/lookup_IDshort.txt'
infile_vcf <- paste0('seqdata/vcf/', fileID, '.vcf.gz')
# infile_cols <- '/home/jelmer/Dropbox/sc_lemurs/metadata/colors/popcols.txt' # Could use to assign colours automatically

outdir_pca <- 'analyses/PCA/dfs/'
outdir_figs <- 'analyses/PCA/plots/'
outdir_snpgds <- 'analyses/PCA/gdsFiles/'
outfile_snpgds <- paste0(outdir_snpgds, '/', fileID, '.gds')

## Packages:
library(gdsfmt)
library(SNPRelate)
library(vcfR)
library(tidyverse)
library(ggpubr)
library(cowplot)

## Read metadata:
lookup <- read.delim(infile_lookup, as.is = TRUE)

## Colours:
#cols.df <- read.delim(infile_cols, header = TRUE, as.is = TRUE)
col.mur <- '#FF00FF'
col.gri <- '#D69B12'
col.hyb <- 'black'

## Create dirs if necessary:
if(!dir.exists(outdir_snpgds)) dir.create(outdir_snpgds, recursive = TRUE)
if(!dir.exists(outdir_pca)) dir.create(outdir_pca, recursive = TRUE)
if(!dir.exists(outdir_figs)) dir.create(outdir_figs, recursive = TRUE)


################################################################################
##### READ VCF #####
################################################################################
snpgdsVCF2GDS(infile_vcf, outfile_snpgds, method = "biallelic.only") # Convert VCF to snpgds format, which is what SNPrelate uses within R
snps <- snpgdsOpen(outfile_snpgds) # Conversion above was to a file, which is now loaded
snpgdsSummary(outfile_snpgds) # Show summary of data in VCF/snpgds


################################################################################
##### SUBSETTING #####
################################################################################
## LD pruning:
#snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
#snpset.id <- unlist(unname(snpset))
#pca <- snpgdsPCA(snpgds.file, snp.id = snpset.id, num.thread = 2)

## List of individuals:
inds_all <- read.gdsn(index.gdsn(snps, "sample.id"))
inds_all

################################################################################
##### ALL INDS #####
################################################################################
subset.ID = 'all'

## Run PCA:
pca_all_raw <- snpgdsPCA(snps,
                         autosome.only = FALSE,
                         num.thread = 1)

## Process results:
species.factor <- c('mgri', 'mmur', 'mhyb') # For non-alphabetical ordering of species
pca_all <- pca.process(pca_all_raw,
                       lookup = lookup,
                       subset.ID = subset.ID,
                       pca.ID = paste0(fileID, '_', subset.ID)) %>%
  mutate(species.short = factor(species.short, levels = species.factor)) %>%
  arrange(species.short)

## eigenvals:
eig_all <- data.frame(PC = 1:length(pca_all_raw$eigenval), eig = pca_all_raw$eigenval)
eig_all <- eig_all[complete.cases(eig_all), ]

## Plot prep:
col.by.labs <- c(expression(italic("griseorufus")),
                 expression(italic("murinus")), 'hybrid')
shape.by.labs <- c('contact zone', 'nearby allopatric')

## Plot:
p1 <- pcplot(pca_all,
             eigenvals = eig_all$eig,
             col.by = 'species.short',
             col.by.name =  'microsats:',
             col.by.labs = col.by.labs,
             cols = c(col.gri, col.mur, col.hyb),
             shape.by = 'poptype',
             shape.by.name = 'site type:',
             shape.by.labs = shape.by.labs,
             shapes = c(1, 0))


################################################################################
##### PCA - MGRI ONLY #####
################################################################################
subset.ID <- 'mgri_cz'

## Select inds:
hybs_mgri <- c('mhyb001', 'mhyb003', 'mhyb004', 'mhyb005', 'mhyb006',
               'mhyb011', 'mhyb015', 'mhyb016')
inds_mgri <- c(inds_all[grepl('mgri', inds_all)],
              inds_all[inds_all %in% hybs_mgri])

## Run PCA:
pca_mgri_raw <- snpgdsPCA(snps,
                          sample.id = inds_mgri, # Subset
                          autosome.only = FALSE,
                          num.thread = 1)

## Process results:
pca_mgri <- pca.process(pca_mgri_raw,
                        lookup = lookup,
                        subset.ID = subset.ID,
                        pca.ID = paste0(fileID, '_', subset.ID))

## eigenvals:
eig_gri <- data.frame(PC = 1:length(pca_mgri_raw$eigenval),
                        eig = pca_mgri_raw$eigenval)
eig_gri <- eig_gri[complete.cases(eig_gri), ]

## Plot prep:
labs_species <- c(expression(italic("griseorufus")), 'hybrid')
title_gri <- expression(paste(italic('griseorufus'), ' only'))

## Plot:
p_mgri <- pcplot(pca_mgri,
                 eigenvals = eig_gri$eig,
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
subset.ID <- 'mmur_cz'

## Select inds:
hybs_mmur <- c('mhyb002', 'mhyb007', 'mhyb008', 'mhyb009',
               'mhyb010', 'mhyb012', 'mhyb013', 'mhyb014')
inds_mmur <- c(inds_all[grepl('mmur', inds_all)],
               inds_all[inds_all %in% hybs_mmur])

## Run PCA:
pca_mmur_raw <- snpgdsPCA(snps,
                          sample.id = inds_mmur,
                          autosome.only = FALSE,
                          num.thread = 1)

## Process results:
pca_mmur <- pca.process(pca_mmur_raw,
                        lookup = lookup,
                        subset.ID = subset.ID,
                        pca.ID = paste0(fileID, '_', subset.ID)) %>%
  mutate(species.short = factor(species.short, levels = c('mmur', 'mhyb')))

## eigenvals:
eig_mur <- data.frame(PC = 1:length(pca_mmur_raw$eigenval),
                        eig = pca_mmur_raw$eigenval)
eig_mur <- eig_mur[complete.cases(eig_mur), ]

## Plot prep:
labs_species <- c(expression('hybrid', italic("murinus")))
title_mur <- expression(paste(italic('murinus'), ' only'))

## Plot:
p_mmur <- pcplot(pca_mmur,
                 eigenvals = eig_mur$eig,
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
