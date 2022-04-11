## Packages
library(here)
library(tidyverse)

## Output file
outfile <- "results/manuscript/TableS2.txt"

## Input files
infile_lookup <- here("metadata/lookup_IDshort.txt")
infile_IDs <- here("metadata/inds_add-to-r03.txt")

## Read and process input files
IDs <- readLines(infile_IDs)

lookup <- read.delim(infile_lookup, as.is = TRUE) %>%
  filter(ID %in% IDs) %>%
  select(ID, Sample_ID, species, site) %>%
  mutate(pop = case_when(species == "ganzhorni" ~ "mur-E",
                         species == "rufus" ~ "NA",
                         species == "murinus" ~ "mur-W",
                         species == "griseorufus" ~ "gri-W"))

## Write to file
write_tsv(lookup, outfile)
