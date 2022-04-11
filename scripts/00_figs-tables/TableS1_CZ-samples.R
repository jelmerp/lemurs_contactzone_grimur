## Packages
library(here)
library(tidyverse)

## Output file
outfile <- here("results/manuscript/TableS1.txt")

## Input files
infile_lookup <- here("metadata/hzlookup_bysample.txt")
infile_hapkeIDs <- here("metadata/Hapke2011/hzlookup_HapkeIDs.txt")
infile_coords <- here("metadata/coordinates/hz_coords_all.txt")
infile_qc <- here("results/qc/qc_byind_r03.txt")

## Read input files
lookup_init <- read.delim(infile_lookup, as.is = TRUE)
hapkeIDs <- read.delim(infile_hapkeIDs, as.is = TRUE) %>%
  select(-Sample_ID)
coords <- read.delim(infile_coords, as.is = TRUE) %>%
  select(ID, lat, lon)
qc <- read.delim(infile_qc, as.is = TRUE) %>%
  select(ID, filter_pass)

## To add BPP/G-PHOCS TRUE/FALSE column (inds used for BPP/G-PHOCS)
coal_inds <- c("mgan007", "mgan008", "mgan014",
               "mgri044", "mgri045", "mgri051",
               "mgri088", "mgri093", "mgri104",
               "mmur001", "mmur006", "mmur009",
               "mmur052", "mmur056", "mmur066")

## Combine dfs
lookup <- lookup_init %>%
  merge(., coords, by = "ID", all.x = TRUE) %>%
  merge(., qc, by = "ID", all.x = TRUE) %>%
  merge(., hapkeIDs, by = "ID", all.x = TRUE) %>%
  mutate(put_hybrid = ifelse(sp_msat_overall == "hybrid", 1, 0),
         lat = round(lat, 7),
         lon = round(lon, 7),
         coal = ifelse(ID %in% coal_inds, TRUE, FALSE)) %>%
  select(ID, Sample_ID, Hapke_ID, site, pop, poptype, lat, lon,
         sp_radseq = sp, sp_msat = sp_msat_nuc, sp_mtDNA,
         filter_pass, put_hybrid)

## Write to file
write_tsv(lookup, outfile)
