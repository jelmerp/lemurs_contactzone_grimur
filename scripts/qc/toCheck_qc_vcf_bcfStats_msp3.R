#!/usr/bin/env Rscript

################################################################################
##### SET-UP #####
################################################################################

setwd('/home/jelmer/Dropbox/sc_lemurs/radseq/')
library(tidyverse); library(reshape2)
source('scripts/qc/qc_vcf_bcfStats_fun.R')

## Command-line arguments:
args <- commandArgs(trailingOnly = TRUE)
software <- args[1]
ref <- args[2]
read.type <- args[3]
vcf.type <- args[4]
whichInds <- args[5]
fileID.prefix <- args[6]
depthDist.ind <- args[7]
basedir <- args[8]

cat('\n##### Running script: qc_vcf.R', '\n')
cat('##### Ref:', ref, '\n')
cat('##### Read type:', read.type, '\n')
cat('##### Vcf type:', vcf.type, '\n')
cat('##### Which inds', whichInds, '\n')
cat('##### File ID prefix:', fileID.prefix, '\n')
cat('##### Plot depth distribution per ind:', depthDist.ind, '\n')
cat('##### Base dir:', basedir, '\n\n')


###### COMMENT OUT #####
software <- 'gatk4'
ref <- 'map2msp3' # mapped2mmur / mapped2cmed
read.type <- 'paired' # R1 / paired
vcf.type <- 'joint' # ind / joint / joint_ind
whichInds <- 'msp3proj' # allInds / Microcebus / Cheirogaleus
fileID.prefix <- 'msp3proj.mac1.FS6'
depthDist.ind <- 'FALSE'
basedir <- 'analyses/qc/vcf/map2msp3.gatk4.paired.joint'
###### COMMENT OUT #####

open.plots <- TRUE

## Inds & read nrs:
inds.df <-  read.delim('metadata/ID.lookupTable.txt', header = TRUE)
#readnrs <- read.delim('analyses/qc/fastq/readnumbers/readnrs_passOnly.txt')
#inds.df <- merge(inds.df, readnrs, by = 'ID', all.x = TRUE)

## Process settings:
setID <- paste0(ref, '.', software, '.', read.type, '.', whichInds)

plotdir <- paste0(basedir, '/figures/')
if(!dir.exists(plotdir)) dir.create(plotdir, recursive = TRUE)

outdir <- paste0(basedir, '/sumtables/')
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

bcftools.input.dir <- paste0(basedir, '/bcftools')

if(vcf.type == 'ind') {
  if(whichInds == 'allInds') IDs <- fileID <- paste0(as.character(inds.df$ID))
  if(whichInds == 'Microcebus') IDs <- fileID <- as.character(inds.df[inds.df$genus == 'Microcebus', ]$ID)
  if(whichInds == 'Cheirogaleus') IDs <- fileID <- as.character(inds.df[inds.df$genus == 'Cheirogaleus', ]$ID)
}
if(vcf.type == 'joint' | vcf.type == 'joint_ind') {
  fileID <- fileID.prefix
  if(whichInds == 'allInds') IDs <- as.character(inds.df$ID)
  if(whichInds == 'Microcebus') IDs <- as.character(inds.df[inds.df$genus == 'Microcebus', ]$ID)
  if(whichInds == 'Cheirogaleus') IDs <- as.character(inds.df[inds.df$genus == 'Cheirogaleus', ]$ID)
  if(whichInds == 'msp3proj') IDs <- as.character(inds.df$ID[which(inds.df$msp3proj == 'msp3proj')])
}


################################################################################
##### PROCESS BCFTOOLS-STATS OUTPUT #####
################################################################################

## Read bcftoolsStats output:
fnames.stats <- paste0(bcftools.input.dir, '/', fileID, '.bcftools.txt')
if(!file.exists(fnames.stats)) cat('ERROR: cant find fnames.stats', fnames.stats, '\n')

if(length(fnames.stats) > 1)
   psc.df <- do.call(rbind, lapply(fnames.stats, get.psc))
if(length(fnames.stats) == 1)
   psc.df <- get.psc(fnames.stats)

psc.df$ID <- gsub('.1.1$', '', psc.df$ID)
psc.df <- psc.df %>% arrange(ID) %>% dplyr::rename(depth = `average depth`)
psc.df$nvar <- psc.df$nRefHom + psc.df$nNonRefHom + psc.df$nHets
psc.df$nvar.nonref <- psc.df$nNonRefHom + psc.df$nHets

psc.df.write <- psc.df %>%
  mutate(TsTv = round(nTransitions / nTransversions, 3)) %>%
  select(ID, nRefHom, nNonRefHom, nHets, TsTv, depth, nSingletons, nvar, nvar.nonref)
psc.df.write$setID <- setID
psc.df.write_filename <- paste0(outdir, '/bcftools_sum_', setID, '.txt')
write.table(psc.df.write, psc.df.write_filename, sep = '\t', quote = FALSE, row.names = FALSE)

## Merge with inds.df:
inds.df.sel <- inds.df %>%
  dplyr::select(ID, ID.short, species.short, species, genus, dnaType, seqType)
psc.df <- merge(psc.df, inds.df.sel, by = 'ID')


################################################################################
##### FAILED AND BEST INDS #####
################################################################################
## Failed inds:
missing <- sort(as.character(IDs[! IDs %in% psc.df$ID]))
psc.df$depth.ok <- ifelse(psc.df$depth < 0.5, 0, 1)
psc.df$nvar.ok <- ifelse(psc.df$nvar < 50000, 0, 1)
failed <- psc.df$ID[which((psc.df$depth.ok + psc.df$nvar.ok) < 2)]

cat('Nr of inds in psc df:', nrow(psc.df), '\n')
cat('Nr of missing inds:', length(missing), '\n')
print(missing)
writeLines(missing, paste0(outdir, '/', fileID, '_inds.missing.txt'))
cat('Nr of failed inds (low depth/nvar):', length(failed), '\n')

## Best inds:
psc.df$nvarByDepth <- round((psc.df$nvar / 10000) * psc.df$depth)
psc.df %>%
  arrange(species, desc(nvarByDepth)) %>%
  select(ID, nvar, depth, nvarByDepth, species.short, seqType)

################################################################################
##### NR OF VARIANTS #####
################################################################################

## Plot histogram of nr of variants:
p <- ggplot(data = psc.df) +
  geom_histogram(aes(nvar / 1000, fill = species.short)) +
  scale_x_continuous(expand = c(0, 1)) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(x = 'Nr of variants (in 1000s)',
       title = labs(title = paste0(setID, ': nvar across samples')))
filename <- paste0(plotdir, fileID.prefix, '_nvar.hist.bySpecies.png')
ggsave(filename, p, width = 8, height = 6)
if(open.plots == TRUE) system(paste('xdg-open', filename))

p <- ggplot(data = psc.df) +
  geom_histogram(aes(nvar / 1000, fill = seqType)) +
  scale_x_continuous(expand = c(0, 1)) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(x = 'Nr of variants (in 1000s)',
       title = labs(title = paste0(setID, ': nvar across samples')))
filename <- paste0(plotdir, fileID.prefix, '_nvar.hist.bySeqType.png')
ggsave(filename, p, width = 8, height = 6)
if(open.plots == TRUE) system(paste('xdg-open', filename))

## Plot histogram of nr of nonref-variants:
p <- ggplot(data = psc.df) +
  geom_histogram(aes(nvar.nonref / 1000, fill = species.short)) +
  scale_x_continuous(expand = c(0, 1)) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(x = 'Nr of variants (in 1000s)',
       title = labs(title = paste0(setID, ': nvar.nonref across samples')))
filename <- paste0(plotdir, fileID.prefix, '_nvar.nonref.hist.bySpecies.png')
ggsave(filename, p, width = 8, height = 6)
if(open.plots == TRUE) system(paste('xdg-open', filename))

p <- ggplot(data = psc.df) +
  geom_histogram(aes(nvar.nonref / 1000, fill = seqType)) +
  scale_x_continuous(expand = c(0, 1)) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(x = 'Nr of variants (in 1000s)',
       title = labs(title = paste0(setID, ': nvar.nonref across samples')))
filename <- paste0(plotdir, fileID.prefix, '_nvar.nonref.hist.bySeqType.png')
ggsave(filename, p, width = 8, height = 6)
if(open.plots == TRUE) system(paste('xdg-open', filename))

## Plot mutations by ind:
plot.muts(psc.df, plotdir, fileID = fileID.prefix,
          ID.start = 1, nr.IDs = 40, legend = TRUE,
          save = TRUE, open.plot = open.plots)
if(nrow(psc.df) > 40) plot.muts(psc.df, plotdir, fileID = fileID.prefix,
                                ID.start = 41, nr.IDs = 40, legend = TRUE,
                                save = TRUE, open.plot = open.plots)


################################################################################
#### DEPTH ####
################################################################################
## Plot depth by ind:
plot.depth(psc.df, plotdir, fileID = fileID.prefix, 1, 40,
           save = TRUE, open.plot = open.plots)
if(nrow(psc.df) > 40) plot.depth(psc.df, plotdir, fileID = fileID.prefix, 41, 40,
                                 save = TRUE, open.plot = open.plots)
if(nrow(psc.df) > 80) plot.depth(psc.df, plotdir, fileID = fileID.prefix, 81, 40,
                                 save = TRUE, open.plot = open.plots)
if(nrow(psc.df) > 120) plot.depth(psc.df, plotdir, fileID = fileID.prefix, 121, 40,
                                  save = TRUE, open.plot = open.plots)
if(nrow(psc.df) > 160) plot.depth(psc.df, plotdir, fileID = fileID.prefix, 161, 44,
                                  save = TRUE, open.plot = open.plots)

## Plot depth histogram:
p <- ggplot(data = psc.df) +
  geom_histogram(aes(depth, fill = species.short)) +
  labs(title = paste0(setID, ': depth across samples'))
filename <- paste0(plotdir, fileID.prefix, '_depth.histogram.png')
ggsave(filename, p, width = 6, height = 6)
if(open.plots == TRUE) system(paste('xdg-open', filename))

## Boxplot of depth by dnaType:
p <- ggplot(data = psc.df) +
  geom_boxplot(aes(x = dnaType, y = depth)) +
  geom_jitter(aes(x = dnaType, y = depth), width = 0.2) +
  labs(x = 'Number of reads (in millions)', y = 'SNP depth',
       title = paste0(setID, ': depth by dnaType'))
filename <- paste0(plotdir, fileID.prefix, '_boxplot_depthVSdnaType.png')
ggsave(filename, p, width = 8, height = 6)
if(open.plots == TRUE) system(paste('xdg-open', filename))

## Plot depth vs nr of variants:
p <- ggplot(data = psc.df) +
  geom_point(aes(x = depth, y = nvar / 1000, colour = species.short)) +
  labs(x = 'Depth', y = 'Number of SNPs (in 1000s)',
       title = paste0(setID, ': nr of reads v nr of variants'))
filename <- paste0(plotdir, fileID.prefix, '_readsVSnvar.bySpecies.png')
ggsave(filename, p, width = 6, height = 6)
if(open.plots == TRUE) system(paste('xdg-open', filename))

p <- ggplot(data = psc.df) +
  geom_point(aes(x = depth, y = nvar / 1000, colour = seqType)) +
  labs(x = 'Depth', y = 'Number of SNPs (in 1000s)',
       title = paste0(setID, ': nr of reads v nr of variants'))
filename <- paste0(plotdir, fileID.prefix, '_readsVSnvar.bySeqType.png')
ggsave(filename, p, width = 6, height = 6)
if(open.plots == TRUE) system(paste('xdg-open', filename))


### Plot depth distribution per ind:
if(vcf.type == 'ind' | depthDist.ind == 'TRUE') {

  if(depthDist.ind == 'TRUE') {
    fileID <- paste0(fileID.prefix, '.', IDs)
    fnames.stats <- paste0(basedir, '/', ref, '/', read.type, '/',
                           vcf.type, '/', fileID, '.bcftoolsStats.txt')
  }

  depth.df <- do.call(rbind, lapply(fnames.stats, get.dp))
  depth.df$ID <- gsub('.*/(.*).bcftoolsStats.txt', '\\1', depth.df$ID)

  plot.depth.dist.wrap <- function(ID.start, nr.IDs = 10) {
    # ID.start = 1; nr.IDs = 10; save.plot = TRUE
    plotdir <- paste0(plotdir, '/depthDist/')
    if(!exists(plotdir)) dir.create(plotdir)

    ID.range <- ID.range <- ID.start:(ID.start + nr.IDs - 1)
    depth.df.sel <- depth.df %>% filter(ID %in% unique(depth.df$ID)[ID.range])
    plot.depth.dist(depth.df.sel, plotdir, fileID = ID.start)
  }
  sapply(seq(1, length(unique(depth.df$ID)), 10), plot.depth.dist.wrap)

}


################################################################################
##### NR OF VARIANTS VS. NR OF READS #####
################################################################################

## Plot nr of reads vs nr of variants:
p <- ggplot(data = psc.df) +
  geom_point(aes(x = reads.passed / 1000000, y = nvar / 1000)) +
  labs(x = 'Number of reads (in millions)', y = 'Number of SNPs (in 1000s)',
       title = paste0(setID, ': nr of reads v nr of variants'))
filename <- paste0(plotdir, fileID.prefix, '_readsVSnvar.png')
ggsave(filename, p, width = 6, height = 6)
if(open.plots == TRUE) system(paste('xdg-open', filename))

## Plot nr of reads vs nr of variants, by genus:
p <- ggplot(data = psc.df) +
  geom_point(aes(x = reads.passed / 1000000, y = nvar / 1000, colour = genus)) +
  labs(x = 'Number of reads (in millions)', y = 'Number of SNPs (in 1000s)',
       title = paste0(setID, ': nr of reads v nr of variants'))
filename <- paste0(plotdir, fileID.prefix, '_readsVSnvar_byGenus.png')
ggsave(filename, p, width = 8, height = 6)
if(open.plots == TRUE) system(paste('xdg-open', filename))

## Plot nr of reads vs nr of variants, by dnaType:
p <- ggplot(data = psc.df) +
  geom_point(aes(x = reads.passed / 1000000, y = nvar / 1000, colour = dnaType)) +
  labs(x = 'Number of reads (in millions)', y = 'Number of SNPs (in 1000s)',
       title = paste0(setID, ': nr of reads v nr of variants'))
filename <- paste0(plotdir, fileID.prefix, '_readsVSnvar_bydnaType.png')
ggsave(filename, p, width = 8, height = 6)
if(open.plots == TRUE) system(paste('xdg-open', filename))


################################################################################
#### DEPTH VS NR OF READS ####
################################################################################
## Plot nr of reads vs depth:
p <- ggplot(data = psc.df) +
  geom_point(aes(x = reads.passed / 1000000, y = depth)) +
  labs(x = 'Number of reads (in millions)', y = 'Read depth',
       title = paste0(setID, ': nr of reads v read depth'))
filename <- paste0(plotdir, fileID.prefix, '_readsVSdepth.png')
ggsave(filename, p, width = 6, height = 6)
if(open.plots == TRUE) system(paste('xdg-open', filename))

## Plot nr of reads vs depth, by genus:
p <- ggplot(data = psc.df) +
  geom_point(aes(x = reads.passed / 1000000, y = depth, colour = genus)) +
  labs(x = 'Number of reads (in millions)', y = 'SNP depth',
       title = paste0(setID, ': nr of reads v, depth by genus'))
filename <- paste0(plotdir, fileID.prefix, '_readsVSdepth_byGenus.png')
ggsave(filename, p, width = 8, height = 6)
if(open.plots == TRUE) system(paste('xdg-open', filename))

## Plot nr of reads vs depth, by dnaType:
p <- ggplot(data = psc.df) +
  geom_point(aes(x = reads.passed / 1000000, y = depth, colour = dnaType)) +
  labs(x = 'Number of reads (in millions)', y = 'SNP depth',
       title = paste0(setID, ': nr of reads v depth, by dnaType'))
filename <- paste0(plotdir, fileID.prefix, '_readsVSdepth_bydnaType.png')
ggsave(filename, p, width = 8, height = 6)
if(open.plots == TRUE) system(paste('xdg-open', filename))
