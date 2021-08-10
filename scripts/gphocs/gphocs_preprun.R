## Set-up:
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/analyses/gphocs/')
source('/home/jelmer/Dropbox/scripts/genomics/gphocs/gphocs_3_prepRun_fun.R')

## Focal file ID:
fileID <- 'hz.mur3gri2c'

indir <- paste0('controlfiles/master/', fileID, '/')
outdir <- paste0('controlfiles/reps/', fileID, '/')
pattern.focal <- '.ctrl'

## Prepare replicates:
file.remove(list.files(indir, pattern = '~', full.names = TRUE))
#file.remove(list.files(outdir, full.names = TRUE))
(masterfiles <- list.files(indir, pattern = pattern.focal,
                          full.names = FALSE, include.dirs = FALSE))

aap <- sapply(masterfiles, prepReps,
              nreps = 4, rm.master = FALSE, dir.source = indir, dir.target = outdir)

