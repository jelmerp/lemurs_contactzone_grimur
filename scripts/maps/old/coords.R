################################################################################
##### SET-UP #####
################################################################################
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')
library(tidyverse)
library(gdata)
source('scripts/maps/coords_fun.R')

## Files:
infile_lookup <- '../radseq/metadata/lookup_IDshort.txt'
infile_joergIDs <- 'metadata/joerg/JoergIDs.txt'

infile_Mtk <- 'metadata/gps/Mangatsiaka_coords.txt'
infile_coords_all <- 'Andohahela_Jacques_coordonnees_Mai2019.xlsx',

outfile_Mtk <- 'metadata/gps/Mangatsiaka_coords_dec.txt'
outfile_Tm <- 'metadata/gps/Tsimelahy_coords_dec.txt'

## Metadata:
lookup <- read.delim(infile_lookup, as.is = TRUE) %>%
  select(Sample_ID, ID.short, site, species.short, species.cor)
jIDs <- read.delim(infile_joergIDs, as.is = TRUE)


################################################################################
##### PROCESS COORDS #####
################################################################################
coords <- read.xls(infile_coords_all, sheet = 'Sheet1')

both <- read.delim(infile_Mtk, as.is = TRUE) %>%
  #mutate(tubeID = gsub(" ", "", tubeID)) %>%
  separate(location, into = c('lat', 'lon'), sep = ', ') %>%
  mutate(lon = sapply(lon, dg2dec),
         lat = -sapply(lat, dg2dec),
         Sample_ID = jIDs$Sample_ID[match(tubeID, jIDs$tubeID)]) %>%
  select(Sample_ID, lat, lon)

both <- merge(both, lookup, by = 'Sample_ID') %>%
  mutate(msat.hybrid = gsub('mmur', 'pure', species.short)) %>%
  mutate(msat.hybrid = gsub('mgri', 'pure', msat.hybrid)) %>%
  mutate(msat.hybrid = gsub('mhyb', 'hybrid', msat.hybrid))

gridloc <- data.frame(Sample_ID = 'gridloc', lat = -24.96458, lon = 46.55468,
                      ID.short = 'gridloc', site = 'Mangatsiaka',
                      species.short = 'mgri', species.cor = 'mgri',
                      msat.hybrid = 'hybrid')
both <- rbind(both, gridloc)


## Write files:
Mtk <- both %>% filter(site == 'Mangatsiaka')
write.table(Mtk, outfile_Mtk, sep = '\t', quote = FALSE, row.names = FALSE)

Tm <- both %>% filter(site == 'Tsimelahy')
write.table(Tm, outfile_Tm, sep = '\t', quote = FALSE, row.names = FALSE)

