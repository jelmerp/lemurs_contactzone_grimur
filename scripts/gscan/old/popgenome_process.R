##### SET-UP #####
setwd('/home/jelmer/Dropbox/sc_fish/cichlids/')
library(data.table); library(dplyr)

triplets <- read.table('analyses/windowstats/popgenome/input/triplets.txt', as.is = TRUE)
triplets <- paste0(triplets$V1, '.', triplets$V2, '.', triplets$V3)
scaffolds <- readLines('metadata/scaffolds.txt')[1:500]

##### APPLY #####
## Set variables:
file.id = 'EjaC.Dstat.DP5.GQ20.MAXMISS0.5.MAF0.01'
winsize <- 50000
stepsize <- 5000

## Run:
file.create('analyses/windowstats/popgenome/missingScafs.txt')
bdf <- lapply(scaffolds, collect.scaffold, file.id = file.id, winsize = winsize, stepsize = stepsize)
bdf <- do.call(rbind, bdf)
#bdf <- collect.scaffold(scaffold = 'NC_022214.1')

filename <- paste0('analyses/windowstats/popgenome/bdf.combined.', file.id, '.win', winsize, '.step', stepsize, '.txt')
write.table(bdf, filename, sep = '\t', quote = FALSE, row.names = FALSE)

missing <- write(unique(readLines('analyses/windowstats/popgenome/missingScafs.txt')), 'analyses/windowstats/popgenome/missingScafs2.txt')

##### SUMMARIZE #####
bdf %>% group_by(pop) %>% summarise(bdf.mean = round(mean(BDF, na.rm = TRUE), 2),
                                    fd.mean = round(mean(fd, na.rm = TRUE), 2),
                                    D.mean = round(mean(D, na.rm = TRUE), 2))

