################################################################################
##### SET-UP #####
################################################################################
## Libraries and scripts:
library(here)
script_fun <- here('scripts/genscriptlink_lemurs/gphocs/gphocs_processlogs_fun.R')
source(script_fun)

## Variables:
setID <- 'hz.mur3gri2c' # (put output withhin folder named "raw" folder within logdir)
gentime <- 3.75
mutrate.gen <- 1.64e-8
mutrate.year <- mutrate.gen / gentime
m.scale <- 1000
t.scale <- 0.0001
burnIn <- 100000
subSample <- 50
lastSample <- NULL

## Input files:
infile_popAbbrev <- here('analyses/gphocs/popinfo/ghocs_popAbbrev.txt')
infile_inds <- here('metadata/lemurs_metadata_link/lookup/lookup_IDshort.txt')
logdir <- paste0(here(), '/analyses/gphocs/output/', setID)

## Read metadata:
pop.lookup <- read.delim(infile_popAbbrev, header = TRUE, as.is = TRUE)
inds <- read.delim(infile_inds, header = TRUE, stringsAsFactors = FALSE)
inds.arr <- inds

## Pops:
if(setID %in% c('r03.all.gphocs2', 'r03.all.gphocs2')) {
  pops <- c('mur', 'gri', 'root')
  kidpops <- c('mur', 'gri')
  parentpops <- 'root'
  allpops <- pops
  currentpops  <- kidpops
}

if(setID == 'hz.mur2gri2c') {
  pops.org <- c('mmur_w', 'mmur_hz', 'mgri_hz', 'mgri_sw', 'anc_mgri', 'anc_mmur', 'anc_root')
  pops <- c('mmur.w', 'mmur.hz', 'mgri.hz', 'mgri.sw', 'anc.mgri', 'anc.mmur', 'anc.root')
  kidpops <- c('murW', 'murHZ', 'griHZ', 'mgriSW', 'a.gri', 'a.mur')
  parentpops <- c('a.mur', 'a.mur', 'a.gri', 'a.gri', 'a.root', 'a.root')
  allpops <- unique(c(kidpops, parentpops))
  currentpops <- kidpops
}

if(setID == 'hz.mur3gri2c') {
  pops.org <- c('anc_mmur_all', 'mmur_w', 'mmur_hz', 'mmur_gan', 'mgri_hz', 'mgri_sw', 'anc_mmur_se', 'anc_mmur', 'anc_mgri', 'anc_root')
  pops <- c('anc.mmur', 'mmur.w', 'mmur.hz',  'mmur.gan', 'mgri.hz', 'mgri.sw',  'anc.mmur.se','anc.mmur',  'anc.mgri', 'anc.root')
  cbind(pops.org, pops)
  kidpops <- c('murW', 'murHZ', 'murGan', 'griHZ', 'mgriSW', 'a.murSE', 'a.gri', 'a.mur')
  parentpops <- c('a.mur', 'a.murSE', 'a.murSE', 'a.gri', 'a.gri', 'a.mur', 'a.root', 'a.root')
  cbind(kidpops, parentpops)
  allpops <- unique(c(kidpops, parentpops))
  currentpops <- kidpops
}


################################################################################
##### PROCESS LOGS #####
################################################################################
Log <- getlogs(logdir = logdir, setID = setID,
               burnIn = burnIn, lastSample = lastSample, subSample = subSample)

## Write table:
saveRDS(Log, file = paste0(logdir, '/mergedlog.RDS'))
mergedLogs.file <- paste0(logdir, '/mergedLogs.txt')
write.table(Log, mergedLogs.file, sep = '\t', quote = FALSE, row.names = FALSE)

