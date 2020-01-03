setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')
library(tidyverse)

sn.file <- '/home/jelmer/Dropbox/sc_lemurs/hybridzone/metadata/r03/samplenames_r03.txt'
sn <- read.delim(sn.file, header = TRUE, as.is = TRUE)

sn %>%
  group_by(species, site) %>%
  tally()

sn %>%
  group_by(site) %>%
  tally()
