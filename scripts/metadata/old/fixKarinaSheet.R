library(gdata)
library(tidyverse)

## My excel with Joerg's metadata:
jsamples <- read.xls('metadata/JoergSamples.xls', sheet = 'main', as.is = TRUE) %>%
  select(Sample_ID, Joerg_ID, Originalbezeichnung, tubeID, DNAconc, sequenced_181127) %>%
  filter(sequenced_181127 == 1)

## New excel from Karina with faulty sequencing info:
newlist <- read.xls('metadata/fixKarinaSheet/Microcebus_sentToDuke_copy.xls', sheet = 'Sheet1', as.is = TRUE) %>%
  select(Lab_ID, species_2018, Site, Year)

## Check if sample IDs match:
#all(jsamples$Joerg_ID %in% newlist$Lab_ID)

## Sequencing statistics:
seqStats <- read.delim('analyses/qc/r03_qc_summary.txt', header = TRUE, as.is = TRUE)
seqStats$Joerg_ID <- jsamples$Joerg_ID[match(seqStats$Sample_ID, jsamples$Sample_ID)]

## Add sequencing statistics to new excel:
newlist$RADseq_done <- ifelse(newlist$Lab_ID %in% seqStats$Joerg_ID, 1, 0)
newlist$RADseq_pass <- seqStats$pass[match(newlist$Lab_ID, seqStats$Joerg_ID)]
newlist$RADseq_rescue <- seqStats$rescue[match(newlist$Lab_ID, seqStats$Joerg_ID)]
newlist$DnaConc_Duke <- jsamples$DNAconc[match(newlist$Lab_ID, jsamples$Joerg_ID)]

## Write table:
write.table(newlist, 'metadata/fixKarinaSheet/Microcebus_sentToDuke_JP.txt',
            sep = '\t', quote = FALSE, row.names = FALSE)
