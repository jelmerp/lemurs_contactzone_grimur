################################################################################
##### SET-UP  #####
################################################################################
library(tidyverse)
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')

infile_samplenames <- '/home/jelmer/Dropbox/sc_lemurs/radseq/metadata/samplenames.txt'
infile_qc <- 'analyses/qc/r03_qc_summary.txt'

outfile_all <- 'analyses/introgress/assignment_files/dummyAssignment_all.txt'
outfile_Mtk <- 'analyses/introgress/assignment_files/dummyAssignment_Mtk.txt'


################################################################################
##### CREATE FILE  #####
################################################################################
## Inds that passed QC:
pass <- read.delim(infile_qc, header = TRUE, as.is = TRUE) %>%
  filter(pass == 1) %>%
  pull(ID)

snames <- read.delim(infile_samplenames, header = TRUE, as.is = TRUE) %>%
  distinct(ID.short, .keep_all = TRUE) %>%
  filter(ID.short %in% pass) %>%
  dplyr::select(ID.short, site) %>%
  mutate(K1 = ifelse(grepl('hyb', ID.short), 0.5, ifelse(grepl('mur', ID.short), 0, 1))) %>%
  mutate(K2 = 1 - K1)

Mtk <- snames %>% filter(site == 'Mangatsiaka')

write.table(snames, outfile_all, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
write.table(Mtk, outfile_Mtk, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
