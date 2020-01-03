################################################################################
##### SET-UP  #####
################################################################################
library(tidyverse)
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')

infile_stats.df <- 'analyses/qc/r03_qc_summary.txt'
infile_lookup <- '/home/jelmer/Dropbox/sc_lemurs/radseq/metadata/lookup_IDshort.txt'
outfile_stats.df <- 'analyses/qc/r03_qc_summary_wSites.txt'

stats.df <- read.delim(infile_stats.df, header = TRUE, as.is = TRUE)
lookup <- read.delim(infile_lookup, header = TRUE, as.is = TRUE) %>%
  select(-Sample_ID, -Sample_ID.org, -dnaType, -genus, -species, -species.org, -msp3proj, -spSite)

ss <- merge(stats.df, lookup, by.x = 'ID', by.y = 'ID.short')
write.table(ss, outfile_stats.df, sep = '\t', quote = FALSE, row.names = FALSE)

################################################################################
##### CHECK  #####
################################################################################
ssp <- ss %>% filter(pass == 1)

(gri <- ssp %>%
    filter(site == 'Tsimelahy' | site == 'Mangatsiaka',
           species.short == 'mgri',
           fmiss.vcf.joint < 0.02))
mean(gri$depth.vcfJoint)


(mur <- ssp %>%
    filter(site == 'Tsimelahy' | site == 'Mangatsiaka',
           species.short == 'mmur',
           fmiss.vcf.joint < 0.02))

mur <- mur %>%
  filter(ID %in% c('mmur052', 'mmur053', 'mmur057', 'mmur061', 'mmur066', 'mmur070'))
mean(mur$depth.vcfJoint)

indSel <- c(mur$ID, gri$ID)
writeLines(indSel, 'metadata/indSel/r03.gphocs1.indsel.txt')
