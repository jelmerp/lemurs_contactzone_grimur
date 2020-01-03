library(gdata); library(dplyr)
setwd('Dropbox/sc_lemurs/radseq/')

inds.df.file <- 'metadata/hybridzone/JoergSamples.xls'
inds.df <- read.xls(inds.df.file, sheet = 'main') %>%
  select(ID_ind, year, site, species, species_mtDNA, hybrid, DNAconc, seq)

inds.df$DNAconc <- as.numeric(inds.df$DNAconc)

length(which(!is.na(inds.df$DNAconc)))
length(which(inds.df$DNAconc < 10))

length(which(inds.df$DNAconc < 10)) / length(which(!is.na(inds.df$DNAconc)))
length(which(inds.df$DNAconc < 5)) / length(which(!is.na(inds.df$DNAconc)))

inds.df %>%
  filter(!is.na(DNAconc)) %>%
  group_by(site, year) %>%
  summarize(mean = mean(DNAconc, na.rm = TRUE))
