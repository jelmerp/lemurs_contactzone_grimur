################################################################################
##### SET-UP #####
################################################################################
## Libraries:
library(gdata)
library(tidyverse)
library(here)
setwd(basedir)

## Input files:
infile_inds <- here('metadata/sampleinfo_hzproj.xls')

## Output files:
outfile_hzproj <- here('metadata/indsel/stacks/hzproj.txt')
outfile_mur3gri2c <- here('metadata/indsel/stacks/hzproj.mur3gri2c.txt')
outfile_mur3gri2ruf <- here('metadata/indsel/stacks/hzproj.mur3gri2ruf.txt')

## Read metadata:
inds <- read.xls(infile_inds, sheet = 'Sheet1')


################################################################################
##### CREATE POPMAPS #####
################################################################################
## All individuals:
(inds_hzproj <- inds %>%
  select(ID, pop2))

## G-phocs selection with outgroup:
(inds_mur3gri2ruf <- inds %>%
  filter(mur3gri2c == 1 | ID == 'mruf008') %>%
  select(ID, pop2))

## Without outgroup:
(inds_mur3gri2c <- inds_mur3gri2ruf %>%
  filter(ID != 'mruf008'))


################################################################################
##### WRITE FILES #####
################################################################################
write.table(inds_hzproj, outfile_hzproj,
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(inds_mur3gri2ruf, outfile_mur3gri2ruf,
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(inds_mur3gri2c, outfile_mur3gri2c,
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
