################################################################################
##### SET-UP #####
################################################################################
## Libraries and scripts
library(tidyverse)
library(here)

## Input files and dirs
infile_lookup <- here('metadata/hzlookup_bysample.txt')
infile_indlist.wOut <- xx
outfile_indfile <- here('results/admixtools/input/indfile_r03.wOut.txt')

## Read metadata
ID.lookup <- read.delim(infile_lookup, header = TRUE, as.is = TRUE)

## Indlists
indlist.wOut <- readLines(infile_indlist.wOut)
indlist.wOut.short <- substr(indlist.wOut, 1, 7)


################################################################################
##### INDFILE ####
################################################################################
species <- substr(indlist.wOut.short, 1, 4)
pop <- ID.lookup$supersite2[match(indlist.wOut.short, ID.lookup$ID.short)]
sex <- 'U'
r03.wOut <- data.frame(indlist.wOut, sex, pop)
table(r03.wOut$pop)

write_tsv(r03.wOut, outfile_indfile)
