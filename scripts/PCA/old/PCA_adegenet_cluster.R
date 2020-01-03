#!/usr/bin/env Rscript

# Usage: Rscript scripts/PCA/PCA_adegenet_Beza.R  Microcebus.r01.FS8.mac3.griseorufus
# For larger datasets: pca <- snpgdsPCA(snps)

################################################################################
##### SET-UP #####
################################################################################
## Command-line args:
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)

file.ID <- args[1]
basedir <- args[2]
ID.type <- args[3]
keep.set <- args[4]
file_inds.keep <- args[5]
file_lookup.IDshort <- args[6]
file_cols <- args[7]
vcf.dir <- args[8]

## Base dir and libraries:
setwd(basedir)
library(vcfR)
library(adegenet)
library(dplyr)

## Other scripts:
source('/datacommons/yoderlab/users/jelmer/radseq/scripts/PCA/PCA_adegenet_fun.R')

## Report:
cat('##### PCA_adegenet_cluster.R: File ID:', file.ID, '\n')
cat('##### PCA_adegenet_cluster.R: Base dir:', basedir, '\n')
cat('##### PCA_adegenet_cluster.R: ID type:', ID.type, '\n')
cat('##### PCA_adegenet_cluster.R: Set of inds to keep:', keep.set, '\n')
cat('##### PCA_adegenet_cluster.R: File with inds to keep (if any):', file_inds.keep, '\n')
cat('##### PCA_adegenet_cluster.R: File with ID lookup table:', file_lookup.IDshort, '\n')
cat('##### PCA_adegenet_cluster.R: File with species colors:', file_cols, '\n')
cat('##### PCA_adegenet_cluster.R: Vcf dir:', vcf.dir, '\n\n')

## Read metadata:
cat("##### PCA_adegenet_cluster.R: Reading cols.file ...\n")
cols.df <- read.delim(file_cols, header = TRUE, as.is = TRUE)

cat("##### PCA_adegenet_cluster.R: Reading lookup file...\n\n")
lookup.IDshort <- read.delim(file_lookup.IDshort, as.is = TRUE)

#inds.df$col.sp <- cols.sp$colour[match(inds.df$species.short, cols.sp$species.short)]

## Output files:
pca.df.dir <- 'analyses/PCA/dfs/'
file_pca.df <- paste0(pca.df.dir, '/', file.ID, '_', keep.set, '.txt')
file_pca.df.intermed <- paste0(pca.df.dir, '/', file.ID, '_', keep.set, '.raw.txt')


################################################################################
##### IMPORT SNPS #####
################################################################################
cat("################################################################################\n")
cat("##### PCA_adegenet_cluster.R: Reading vcf file and converting to genlight object...\n")
vcf.file <- paste0(vcf.dir, '/', file.ID, '.vcf.gz')
vcf <- read.vcfR(vcf.file)
snps <- vcfR2genlight(vcf)

if(keep.set != 'all') {
  cat("##### PCA_adegenet_cluster.R: Reading file_inds.keep and subsetting inds...\n")
  inds.keep <- readLines(file_inds.keep)
  keep.rows <- rownames(as.matrix(snps)) %in% inds.keep
  snps <- new('genlight', as.matrix(snps)[keep.rows, ])
}


################################################################################
##### RUN PCA #####
################################################################################
cat("#######################################################\n")
cat("##### PCA_adegenet_cluster.R: Running PCA...\n")
pca.res <- glPca(snps, center = TRUE, scale = TRUE, nf = 4)

cat("##### PCA_adegenet_cluster.R: Processing PCA df...\n")
pca.df <- pca.process(pca.res, ID.type = ID.type,
                      subset.ID = keep.set, file_pca.df = file_pca.df)
cat("##### PCA_adegenet_cluster.R: Dimensions of pca.df:", dim(pca.df), '\n')

cat("##### PCA_adegenet_cluster.R: Done with script.')

