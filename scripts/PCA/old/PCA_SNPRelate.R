## TO DO: PLOT SAMPLE NAMES

################################################################################
##### SET-UP #####
################################################################################
## Variables:
fileID <- 'r03.all.mac1.FS6'
#fileID <- 'r03.all.rawSNPs.maxmiss75.mac3'
#fileID <- 'r03.all.rawSNPs.goodInds.maxmiss75.mac3'

basedir <- '/home/jelmer/Dropbox/sc_lemurs/hybridzone/'
ID.type <- 'ID.short'
file_lookup.IDshort <- '/home/jelmer/Dropbox/sc_lemurs/radseq/metadata/lookup_IDshort.txt'
file_inds.keep <- 'notany'
file_cols <- '/home/jelmer/Dropbox/sc_lemurs/metadata/colors/colors.species.txt'
keep.set <- 'all'
vcf.dir <- 'seqdata/vcf/'

## Base dir and libraries:
setwd(basedir)
library(gdsfmt)
library(SNPRelate)
library(vcfR)
library(tidyverse)
library(ggpubr)
library(cowplot)

## Other scripts:
source('/home/jelmer/Dropbox/sc_lemurs/scripts/PCA/PCA_R_fun.R')

cat("##### Reading cols.file ...\n")
cols.df <- read.delim(file_cols, header = TRUE, as.is = TRUE)
col.mur <- '#FF00FF'
col.gri <- '#D69B12'
col.hyb <- 'black'

cat("##### Reading file_inds.keep\n")
if(keep.set != 'all') inds.keep <- readLines(file_inds.keep)

cat("##### Reading lookup file...\n")
lookup.IDshort <- read.delim(file_lookup.IDshort, as.is = TRUE)

## Files:
vcf.file <- paste0(vcf.dir, '/', fileID, '.vcf.gz')
pca.df.dir <- 'analyses/PCA/dfs/'


################################################################################
##### READ VCF #####
################################################################################
cat("#######################################################\n")
cat("##### Reading vcf file...\n")

snpgds.file <- paste0('analyses/PCA/gdsFiles/', fileID, '.gds')
snpgdsVCF2GDS(vcf.file, snpgds.file, method = "biallelic.only")
snpgdsSummary(snpgds.file)
snps <- snpgdsOpen(snpgds.file)


################################################################################
##### SUBSETTING #####
################################################################################
subset.ID = 'all'
outfile_pca.df <- paste0(pca.df.dir, '/', fileID, '_', subset.ID, '.txt')
outfile_eigenvalues <- paste0(pca.df.dir, '/', fileID, '_', subset.ID, '_eigenvalues.txt')

## HOW TO SUBSET INDS?
#if(keep.set != 'all') {
#  keep.rows <- rownames(as.matrix(snps)) %in% inds.keep
#  snps <- new('genlight', as.matrix(snps)[keep.rows, ])
#}

## LD pruning:
#snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
#snpset.id <- unlist(unname(snpset))
#pca <- snpgdsPCA(snpgds.file, snp.id = snpset.id, num.thread = 2)


################################################################################
##### PCA FOR ALL INDS #####
################################################################################
cat("#######################################################\n")
cat("##### Running PCA...\n")

## Run PCA:
pca.res <- snpgdsPCA(snps, autosome.only = FALSE, num.thread = 1)
# missing.rate=

## Process results:
cat("##### Editing PCA df...\n")
pca.df <- pca.process(pca.res, inds.df = lookup.IDshort, source = 'snpgdsPCA',
                      subset.ID = subset.ID, pca.df.ID = paste0(fileID, '_', subset.ID)) %>%
  mutate(supersite = gsub('mmur.hz', 'contact zone', supersite)) %>%
  mutate(supersite = gsub('mgri.hz', 'contact zone', supersite)) %>%
  mutate(supersite = gsub('mhyb.hz', 'contact zone', supersite)) %>%
  mutate(supersite = gsub('mmur.se', 'nearby allopatric', supersite)) %>%
  mutate(supersite = gsub('mgri.se', 'nearby allopatric', supersite)) %>%
  mutate(species.short = factor(species.short, levels = c('mgri', 'mmur', 'mhyb'))) %>%
  arrange(species.short)
cat("##### Dimensions of pca.df:", dim(pca.df))
#head(pca.df)

eigenvalues <- data.frame(PC = 1:length(pca.res$eigenval), eig = pca.res$eigenval)

## Labels:
col.by.labs <- c(expression(italic("griseorufus")), expression(italic("murinus")), 'hybrid')
shape.by.labs <- c('contact zone', 'nearby allopatric')
#my.title <- 'Only filtering for missing data\n(poor inds removed)'
#my.title <- 'Only filtering for missing data\n(poor inds kept)'
my.title <- 'Regular filtering pipeline'

p1 <- pcplot(pca.df, pca.ID = paste0(fileID, '_plotAll_PC12'), eigenvalues = eigenvalues$eig,
             col.by = 'species.short', col.by.name =  'previous\nmicrosatellite\nassignment:',
             col.by.labs = col.by.labs, cols = c(col.gri, col.mur, col.hyb),
             shape.by = 'supersite', shape.by.name = 'site type',
             shape.by.labs = shape.by.labs, shapes = c(1, 0),
             plot.title = my.title)

# p34 <- pcplot(pca.df, pca.ID = paste0(fileID, '_plotAll_PC34'),
#               pcX = 'PC3', pcY = 'PC4', col.by = 'species', col.by.name = 'species',
#               shape.by = 'site', shape.by.name = 'site', provide.cols = TRUE, dotsize = 3,
#               plot.title = 'All species: PC3 and PC4')


###########################################################################
##### PCA - MGRI ONLY #####
###########################################################################
mgri.hybs <- c('mhyb001', 'mhyb003', 'mhyb004', 'mhyb005', 'mhyb006',
               'mhyb011', 'mhyb015', 'mhyb016')

IDs.all <- rownames(as.matrix(snps))
IDs.mgri <- c(IDs.all[grepl('mgri', IDs.all)], IDs.all[IDs.all %in% mgri.hybs])

keep.rows <- rownames(as.matrix(snps)) %in% IDs.mgri
snps.mgri <- new('genlight', as.matrix(snps)[keep.rows, ])

pca.df.mgri <- glPca(snps.mgri, center = TRUE, scale = TRUE, nf = 4) %>%
  pca.process(subset.ID = 'mgri', ID.type = ID.type)

pc12.mgri <- pcplot(pca.df.mgri, pca.ID = paste0(fileID, '_mgri'),
                      col.by = 'site', col.by.name = 'site',
                      shape.by = 'species', shape.by.name = 'species',
                      provide.cols = FALSE, dotsize = 3,
                      plot.title = 'M. griseorufus(-like) only')


###########################################################################
##### PCA- MMUR ONLY #####
###########################################################################
mmur.hybs <- c('mhyb002', 'mhyb007', 'mhyb008', 'mhyb009', 'mhyb010',
               'mhyb012', 'mhyb013', 'mhyb014')
IDs.all <- rownames(as.matrix(snps))
IDs.mmur <- c(IDs.all[grepl('mmur', IDs.all)], IDs.all[IDs.all %in% mmur.hybs])

keep.rows <- rownames(as.matrix(snps)) %in% IDs.mmur
snps.mmur <- new('genlight', as.matrix(snps)[keep.rows, ])

pca.df.mmur <- glPca(snps.mmur, center = TRUE, scale = TRUE, nf = 4) %>%
  pca.process(subset.ID = 'mmur', ID.type = ID.type)

pc12.mmur <- pcplot(pca.df.mmur, pca.ID = paste0(fileID, '_mmur'),
                      col.by = 'site', col.by.name = 'site',
                      shape.by = 'species', shape.by.name = 'species',
                      provide.cols = FALSE, dotsize = 3,
                      plot.title = 'M. murinus(-like) only')


###########################################################################
##### PCA - MANGATSIAKA ONLY #####
###########################################################################
Mtk <- lookup.IDshort$ID.short[which(lookup.IDshort$site == 'Mangatsiaka')]
IDs.all <- rownames(as.matrix(snps))
IDs.Mtk <- IDs.all[IDs.all %in% Mtk]
keep.rows <- rownames(as.matrix(snps)) %in% IDs.Mtk
snps.Mtk <- new('genlight', as.matrix(snps)[keep.rows, ])

pca.df.Mtk <- glPca(snps.Mtk, center = TRUE, scale = TRUE, nf = 4) %>%
  pca.process(subset.ID = 'Mtk', ID.type = ID.type)

pc12.Mtk <- pcplot(pca.df.Mtk, pca.ID = paste0(fileID, '_Mtk'),
                     col.by = 'species', col.by.name = 'species',
                     provide.cols = FALSE, dotsize = 3,
                     plot.title = 'Mangatsiaka only')


###########################################################################
##### PCA - MANGATSIAKA ONLY #####
###########################################################################
Tml <- lookup.IDshort$ID.short[which(lookup.IDshort$site == 'Tsimelahy')]
IDs.all <- rownames(as.matrix(snps))
IDs.Tml <- IDs.all[IDs.all %in% Tml]
keep.rows <- rownames(as.matrix(snps)) %in% IDs.Tml
snps.Tml <- new('genlight', as.matrix(snps)[keep.rows, ])

pca.df.Tml <- glPca(snps.Tml, center = TRUE, scale = TRUE, nf = 4) %>%
  pca.process(subset.ID = 'Tml', ID.type = ID.type)

pc12.Tml <- pcplot(pca.df.Tml, pca.ID = paste0(fileID, '_Tml'),
                     col.by = 'species', col.by.name = 'species',
                     provide.cols = FALSE, dotsize = 3,
                     plot.title = 'Tsimelahy only')


###########################################################################
##### COMBINE PLOTS #####
###########################################################################
## Arrange  panels:
p <- ggarrange(pcaplot.all.pc12, pcaplot.all.pc34, pcaplot.mmacmsp3,
               ncol = 3, nrow = 1, widths = c(1, 1, 1.05))
p <- p + draw_plot_label(label = c('A', 'B', 'C'), size = 20,
                         x = c(0, 0.33, 0.67), y = c(1, 1, 1))

## Save as png:
figfile <- paste0('analyses/PCA/msp3/', fileID, '_combined2.png')
ggexport(p, filename = figfile, width = 1100, height = 300)
system(paste0('xdg-open ', figfile))

# figfile.pdf <- paste0('analyses/PCA/msp3/', fileID, '_combined2.eps')
# ggexport(p, filename = figfile.pdf, width = 11, height = 3)
# system(paste0('xdg-open ', figfile.pdf))
