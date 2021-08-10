## Packages
library(here)
library(tidyverse)
library(data.table)

## Define input files
infile_gphocs <- here('results/gphocs/output/gphocs_mergedlogs.txt')
infile_bpp <- here('results/bpp/output/bpp_mergedlogs.txt')

## Define output file
outfile_logs <- here('results/gphocs/output/gphocs-bpp_mergedlogs.txt')

## Read GPHOCS and BPP output files
log_gphocs <- as.data.frame(fread(infile_gphocs, stringsAsFactors = FALSE))
log_bpp <- as.data.frame(fread(infile_bpp, stringsAsFactors = FALSE))

## Define populations
poplist <- list(murW = 'A', B = 'A',
                murC = 'B', murE = 'B',
                griC = 'C', griW = 'C',
                A = 'root', C = 'root')
childpops <- names(poplist)
parentpops <- unique(as.character(poplist))
allpops <- unique(c(childpops, parentpops))

migorder <- c('griC_2_murC', 'murC_2_griC', 'C_2_B', 'B_2_C', 'C_2_A', 'A_2_C')

## Combine GPHOCS and BPP output files
logs <- bind_rows(log_gphocs, log_bpp) %>%
  filter(runID %in% c('noMig', 'g2m2anc', 'bpp_6mig'),
         !(runID == 'bpp_6mig' & rep == 4),
         var %in% c('theta', 'tau', '2Nm', 'm', 'mprop', 'phi')) %>%
  mutate(pop = factor(pop, levels = allpops),
         migpattern = factor(migpattern, levels = migorder))

## Remove tau/popsizes at introgression events
logs <- logs[-intersect(which(is.na(logs$pop)), which(is.na(logs$migfrom))), ]

## Write output to file
write.table(logs, outfile_logs,
            row.names = FALSE, quote = FALSE, sep = '\t')
