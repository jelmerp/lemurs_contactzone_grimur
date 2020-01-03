################################################################################
##### SET-UP #####
################################################################################
## Libraries and scripts:
library(tidyverse)
library(here)

## Input files and dirs:
infile_lookup <- here('metadata/radseq_metadata_link/lookup_IDshort.txt')
infile_indlist.wOut <- xx
outfile_indfile <- here('analyses/admixtools/input/indfile_r03.wOut.txt')

## Read metadata:
ID.lookup <- read.delim(infile_lookup, header = TRUE, as.is = TRUE)

## Indlists:
indlist.wOut <- readLines(infile_indlist.wOut)
indlist.wOut.short <- substr(indlist.wOut, 1, 7)


################################################################################
##### INDFILE ####
################################################################################
species <- substr(indlist.wOut.short, 1, 4)
pop <- ID.lookup$supersite2[match(indlist.wOut.short, ID.lookup$ID.short)]
#pop[which(is.na(pop))] <- species[which(is.na(pop))]
sex <- 'U'
r03.wOut <- data.frame(indlist.wOut, sex, pop)
table(r03.wOut$pop)

write.table(r03.wOut, outfile_indfile, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)


################################################################################
##### D4 popfile #####
################################################################################
#outfile_popfile <- 'analyses/admixtools/input/outfile_popfile_dstat_msp3.msp3pops.txt'

#unique(msp3.df$sp.pop)
# pop1 <- rep('mmac', 10)
# pop2 <- c(rep('msp3.Anjiahely', 8), 'msp3.Ambavala', 'msim')
# pop3 <- c(rep(c('msp3.Mananara_Nord', 'msp3.Ambavala', 'msp3.Antanambe', 'msp3.Antsiradrano'), 2), 'msim', 'mmit')
# pop4 <- c(rep('mmur', 4), c(rep('msim', 4)), 'mmur', 'mmur')
#
# outfile_popfile.d.df <- data.frame(pop1, pop2, pop3, pop4)
# write.table(outfile_popfile.d.df, outfile_popfile,  sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)


################################################################################
##### F3 outfile_popfile #####
################################################################################
# pop1 <- rep('mmac', 6)
# pop2 <- c(rep('msp3.Anjiahely', 4), 'msp3.Ambavala', 'msim')
# pop3 <- c('msp3.Mananara_Nord', 'msp3.Ambavala', 'msp3.Antanambe', 'msp3.Antsiradrano', 'msim', 'mmit')
#
# outfile_popfile.f3.df <- data.frame(pop1, pop2, pop3)
# outfile_popfile <- 'analyses/admixtools/input/outfile_popfile_f3stat_msp3.msp3pops.txt'
# write.table(outfile_popfile.f3.df, outfile_popfile,
#             sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

