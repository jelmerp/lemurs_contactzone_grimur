## NOTE THAT TWO SAMPLES WERE SEQUENCED TWICE IN R03: mmur045 & mmur052

# SET-UP -----------------------------------------------------------------------

## Packages
library(here)
library(tidyverse)
library(readxl)

## Input files - RADseq metadata:
infile_lookup_short <- here("metadata/links/lookup_IDshort.txt")   # General mouse lemur RADseq lookup table 1
infile_lookup_long <- here("metadata/links/lookup_IDlong.txt")     # General mouse lemur RADseq lookup table 2
infile_sites <- here("metadata/links/sites_pops.txt")
infile_labwork <- here("metadata/labwork/JoergSamples_sequenced.xls")

## Input files - focal samples:
infile_IDs_r03 <- here("metadata/samplelists/sampleIDsShort_r03.txt")
infile_IDs_add <- here("metadata/samplelists/hzproj_add-to-r03_IDs.txt")

## Output files:
outfile_lookup_short <- here("metadata/hzlookup_bysample.txt")
outfile_lookup_long <- here("metadata/hzlookup_bylibrary.txt")

## Read input files:
lookup_short_raw <- read.delim(infile_lookup_short, as.is = TRUE)
lookup_long <- read.delim(infile_lookup_long, as.is = TRUE) %>%
  rename(ID_short = ID.short)
sites <- read.delim(infile_sites, sep = ",", as.is = TRUE)
inds <- sort(c(readLines(infile_IDs_r03), readLines(infile_IDs_add)))
labwork <- read_xls(infile_labwork) %>%
  select(ID, sp_msat_nuc, sp_mtDNA, sp_msat_overall)


# LOOKUP - ID_short ------------------------------------------------------------

## Process and merge input files
lookup_short <- lookup_short_raw %>%
  filter(ID %in% inds) %>%
  merge(., sites, by = c("site", "sp"), all.x = TRUE) %>%
  rename(pop = pop2) %>%
  replace_na(list(poptype = "allopatric")) %>%
  arrange(ID) %>%
  select(ID, Sample_ID, species, sp, seqruns, site, pop, poptype) %>%
  merge(., labwork, by = "ID")

## Write output file
write.table(lookup_short, outfile_lookup_short,
            sep = "\t", quote = FALSE, row.names = FALSE)


# LOOKUP - ID_long -------------------------------------------------------------

lookup_long <- lookup_long %>%
  filter(ID_short %in% inds) %>%
  arrange(ID)

## Write output file
write.table(lookup_long, outfile_lookup_long,
            sep = "\t", quote = FALSE, row.names = FALSE)
