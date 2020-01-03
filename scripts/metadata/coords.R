################################################################################
##### SET-UP #####
################################################################################
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')
source('scripts/metadata/coords_fun.R')
library(tidyverse)
library(gdata)

## Files:
infile_lookup <- '../radseq/metadata/lookup_IDshort.txt'
infile_joergIDs <- 'metadata/joerg/JoergIDs.txt'

infile_Mtk <- 'metadata/gps/Mangatsiaka_coords.txt'
infile_coords_all <- 'metadata/gps/Andohahela_Jacques_coordonnees_Mai2019.xlsx'

outfile_coords_all <- 'metadata/gps/Andohahela_coords_all.txt'

## Metadata:
lookup <- read.delim(infile_lookup, as.is = TRUE) %>%
  select(Sample_ID, ID.short, species.short, species.cor)
jIDs <- read.delim(infile_joergIDs, as.is = TRUE)


################################################################################
##### PROCESS COORDS #####
################################################################################
coords <- read.xls(infile_coords_all, sheet = 'Sheet1') %>%
  rename(site = Site, lat = Lat, lon = Long) %>%
  mutate(lon = gsub(',', '.', lon),
         lat = gsub(',', '.', lat)) %>%
  mutate(lon = sapply(lon, dg2dec),
         lat = -sapply(lat, dg2dec),
         Sample_ID = jIDs$Sample_ID[match(Lab_ID, jIDs$tubeID)]) %>%
  select(Sample_ID, site, lat, lon)
#coords[which(is.na(coords$Sample_ID)), ]

coords <- merge(coords, lookup, by = 'Sample_ID') %>%
  mutate(msat.hybrid = gsub('mmur', 'pure', species.short)) %>%
  mutate(msat.hybrid = gsub('mgri', 'pure', msat.hybrid)) %>%
  mutate(msat.hybrid = gsub('mhyb', 'hybrid', msat.hybrid))
#head(coords)

## Write files:
write.table(coords, outfile_coords_all,
            sep = '\t', quote = FALSE, row.names = FALSE)
