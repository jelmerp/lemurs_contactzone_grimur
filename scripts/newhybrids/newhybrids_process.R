# SET-UP -----------------------------------------------------------------------
library(tidyverse)

## QC info
infile_qc <- "results/qc/r03_qc_summary_wSites.txt"
qc <- read.delim(infile_qc, header = TRUE, as.is = TRUE) %>%
  filter(pass == 1)
inds.pass <- qc$ID

## Newhybrids output
runID <- "r03.all.test2"
runID <- "r03.all.predef"
runID <- "r03.keepHybs"
infile_newhybrids <- paste0("results/newhybrids/output/", runID, "/aa-PofZ.txt")


# PROCESS NEWHYBRIDS OUTPUT  ---------------------------------------------------
nh <- read.delim(infile_newhybrids, header = TRUE, as.is = TRUE,
                 col.names = c("ind.nr", "ID", "Pure_0", "Pure_1",
                               "F1", "F2", "Bx_0", "Bx_1"))
nh$ID <- inds.pass
nh$sp <- substr(nh$ID, 1, 4)

table(nh$sp, nh$Pure_0)
table(nh$sp, nh$Pure_1)
