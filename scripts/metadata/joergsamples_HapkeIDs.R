################################################################################
##### SET-UP #####
################################################################################
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')
library(gdata)
library(dplyr)

## Files:
joerg.file <- 'metadata/joerg/JoergSamples.xls'
hapke.file <- 'metadata/joerg/Hapke_labID_forJelmer_30.04.19.csv'
gps.hapke.file <- 'metadata/Hapke2011/12862_2011_1900_MOESM1_ESM.XLS'
gps.jacques.file <- 'metadata/gps/Andohahela_Jacques_coordonnees_Mai2019.xlsx'
myIDs.file <- 'metadata/r03/samplenames_r03.txt'
outfile <- 'metadata/joerg/JoergSamples2.csv'


################################################################################
##### READ AND EDIT ALL FILES #####
################################################################################
## Original GPS file:
joerg <- read.xls(joerg.file, sheet = 'main', as.is = TRUE)
joerg <- joerg[, 1:31]

## Missing gps:
gps.missing <- c('2006', '2009', '2019', '2020', '2021', '2026', '2027', '2033',
                 '2034', '886', '887', '889', '890')
gps.jacques <- read.xls(gps.jacques.file, as.is = TRUE) %>%
  filter(Lab_ID %in% gps.missing) %>%
  rename(lat = Lat, long = Long) %>%
  select(Lab_ID, lat, long)
joerg$lat[match(gps.jacques$Lab_ID, joerg$Sample_ID)] <- gps.jacques$lat
joerg$long[match(gps.jacques$Lab_ID, joerg$Sample_ID)] <- gps.jacques$long

## Hapke:
hapke <- read.csv(hapke.file, as.is = TRUE) %>%
  rename(Hapke_ID = Hapke_specimen_voucher, Hapke_mtID = mt_haplotype_species) %>%
  select(-Region, -Site) %>%
  mutate(Hapke_mtID = gsub('Microcebus ', '', Hapke_mtID))
#hapke$Lab_ID %in% joerg$Joerg_ID2
#hapke$Lab_ID[! hapke$Lab_ID %in% joerg$Joerg_ID]

## Hapke GPS:
gps.hapke <- read.xls(gps.hapke.file, as.is = TRUE) %>%
  rename(lat.hapke = Latitude..WGS84., lon.hapke = Longitude..WGS.84.,
         Hapke_ID = Specimen.voucher) %>%
  select(-Site, -mt_haplotype_assigned_to_species)

## IDs:
myIDs <- read.delim(myIDs.file, as.is = TRUE) %>%
  select(ID.short, Sample_ID) %>%
  distinct(ID.short, .keep_all = TRUE)


################################################################################
##### MERGE #####
################################################################################
joerg.hapke <- merge(joerg, hapke, by.x = 'Joerg_ID2', by.y = 'Lab_ID', all.x = TRUE)
joerg.hapke.gps <- merge(joerg.hapke, gps.hapke, by = 'Hapke_ID', all.x = TRUE)
all <- merge(myIDs, joerg.hapke.gps, by = 'Sample_ID')
#myIDs$Sample_ID[!myIDs$Sample_ID %in% joerg.hapke.gps$Sample_ID]
#head(all)


################################################################################
##### WRITE FILE #####
################################################################################
#write.table(all, outfile, sep = '\t', quote = FALSE, row.names = FALSE)
write.csv(all, outfile, row.names = FALSE)
