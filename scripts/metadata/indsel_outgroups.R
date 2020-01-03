library(tidyverse)
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone')

################################################################################
##### R03.WITH-OUTGROUPS #####
################################################################################
## Read metadata:
IDlookup.file <- '/home/jelmer/Dropbox/sc_lemurs/radseq/metadata/lookup_IDshort.txt'
sn <- read.delim(IDlookup.file, header = TRUE, as.is = TRUE)

## Select pops:
supersites.gri <- c('mgri_sw', 'mgri_se', 'hybridzone')

sites.gri <- 'Beza_Mahafaly_site3_P2_Spiny_Forest'
sites.gan <- 'Mandena'
sites.man <- 'Bemanasy'
sites.all <- c(sites.gri, sites.gan, sites.man)

bad.inds <- c('mmur003', 'mmur011', 'mmur019', 'mgri042')

sn.sel <- sn %>%
  filter(site %in% sites.all | ID.short == 'mruf007' | species.short == 'mmur' | supersite %in% supersites.gri,
         ! ID.short %in% bad.inds) %>%
  distinct(ID.short, .keep_all = TRUE) %>%
  arrange(ID.short)

## Check:
sn.sel %>%
  select(ID.short, site, species) %>%
  group_by(species, site) %>%
  tally() %>%
  print(n = 50)

## Write file:
IDs.sel <- sn.sel$ID.short
IDs.sel.file <- 'metadata/indSel/r03.wOutgroups.IDs.txt'
writeLines(IDs.sel, IDs.sel.file)

# sn %>% filter(species == 'griseorufus') %>% group_by(loc) %>% tally()
# sn %>% filter(species == 'murinus') %>% group_by(loc) %>% tally()
# sn %>% filter(species == 'ganzhorni') %>% group_by(loc) %>% tally()
# sn %>% filter(species == 'manitatra') %>% group_by(loc) %>% tally()
# sn %>% filter(species == 'rufus') %>% group_by(loc) %>% tally()


################################################################################
##### KEEP HYBRIDS #####
################################################################################
r03.file <- 'metadata/r03/samplenames_r03.txt'
r03 <- read.delim(r03.file, header = TRUE, as.is = TRUE)

r03 %>% select(Sample_ID, ID.short, species, site) %>%
  filter(Sample_ID %in% c('886', '887', '890') | species == 'hybrid') %>%
  pull(ID.short)
