################################################################################
##### SET-UP #####
################################################################################
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')
library(tidyverse)

## Files:
infile_lookup <- '/home/jelmer/Dropbox/sc_lemurs/radseq/metadata/lookup_IDshort.txt'
infile_focInds <- 'analyses/qc/vcf/map2mmur.gatk4.paired.joint/filtering/r03.wOutgroups.FS6_indlist.txt'

outfile_lookup <- 'metadata/lookup_r03.wOutgroups.FS6.txt'
outfile_indsel1 <- 'metadata/indSel/r03.wOutgroups.indsel1.txt'
outfile_indsel2 <- 'metadata/indSel/r03.wOutgroups.indsel2.txt'
outfile_indsel3 <- 'metadata/indSel/r03.wOutgroups.indsel3.txt'


################################################################################
##### CREATE LOOKUP #####
################################################################################
## Focal inds:
focInds <- readLines(infile_focInds)
focInds.short <- substr(focInds, 1, 7)

## Metadata master file:
lookup <- read.delim(infile_lookup, header = TRUE, as.is = TRUE) %>%
  filter(ID.short %in% focInds.short) %>%
  arrange(ID.short) %>%
  mutate(ID = focInds) %>%
  select(ID, ID.short, Sample_ID, species.short, site, spSite, supersite)

## Supersite ID for all:
lookup$supersite[which(is.na(lookup$supersite))] <- lookup$spSite[which(is.na(lookup$supersite))]
table(lookup$supersite)

## supersite2:
lookup$supersite2 <- lookup$supersite
lookup$supersite2[which(lookup$supersite == 'mmur_sw')] <- NA
lookup$supersite2[which(lookup$species.short == 'mman')] <- NA
lookup$supersite2[which(lookup$supersite == 'mgri_beza')] <- NA
#lookup$supersite2[which(lookup$species.short == 'mhyb')] <- NA
lookup$supersite2 <- gsub('mgan_Mandena', 'mmur_gan', lookup$supersite2)
lookup$supersite2 <- gsub('mruf_Ranomafana', 'mruf', lookup$supersite2)
cat('Table of supersite2 variable:\n')
table(lookup$supersite2)


################################################################################
##### ADD QC #####
################################################################################
infile_imiss <- 'analyses/qc/vcf/map2mmur.gatk4.paired.joint/vcftools/r03.wOutgroups.mac1.FS6.imiss'
imiss <- read.delim(infile_imiss, as.is = TRUE,
                    col.names = c('ID', 'nSNPs', 'nFilt', 'nMiss', 'perc.miss')) %>%
  select(ID, perc.miss) %>%
  mutate(perc.miss = round(perc.miss * 100, 3)) %>%
  filter(ID %in% lookup$ID) # mgri100

infile_idepth <- 'analyses/qc/vcf/map2mmur.gatk4.paired.joint/vcftools/r03.wOutgroups.mac1.FS6.idepth'
idepth <- read.delim(infile_idepth, as.is = TRUE,
                     col.names = c('ID', 'nSites', 'depth')) %>%
  mutate(depth = round(depth, 2))

qc <- merge(idepth, imiss, by = 'ID')
lookup.qc <- merge(lookup, qc, by = 'ID')
head(lookup.qc)


################################################################################
##### WRITE FILES #####
################################################################################
write.table(lookup.qc, outfile_lookup, sep = '\t', quote = FALSE, row.names = FALSE)

## Indsel:
indsel1 <- lookup %>% filter(!is.na(supersite2)) %>% pull(ID) # Excludes mman, Beza, mmur-Vohimena
indsel2 <- indsel1[-which(indsel1 == 'mruf007_r01_p3d12')] # no rufus outgroup, for e.g. ADMXITURE
indsel3 <- c(indsel2, 'mgri100', # present in mac1 not in mac3
             'mruf003', 'mruf007', 'mruf008')
#'mtan002', 'mtan005', 'mtan007') # for new VCF
writeLines(indsel1, outfile_indsel1)
writeLines(indsel2, outfile_indsel2)
writeLines(indsel3, outfile_indsel3)
