# SET-UP -----------------------------------------------------------------------

## Packages
library(here)

## Other scripts
source(here('scripts/gphocs-bpp/gphocs_processlogs_fun.R'))

## Define variables
set.seed(7643)
setID <- 'hz.mur3gri2c'

mutrate_gen <- 1.52                    # Mutation rate per generation (will be multiplied by 1e-8)
mutrate_var <- 0.000001                # Variance for mutation rate  (will be multiplied by 1e-8)
gentime <- 3.5                         # Generation time
gentime_sd <- 1.16                     # SD for generation time

m_scale <- 1000                        # Gphocs scaling of m (migration rate)
t_scale <- 0.0001                      # Gphocs scaling of theta and tau
burnin <- 70000                        # Size of burn-in, i.e. iterations to remove
subsample <- 50                        # Subsample 1 in x samples (output lines). Default: 50
last_sample <- NULL                    # Stop processing log at sample x

rename_pops_before <- TRUE             # Rename populations
rename_pops_after <- FALSE             # Rename populations

## Define input files
infile_lookup <- here('results/gphocs/metadata/gphocs_poplookup.txt')
logdir <- here('results/gphocs/output/raw/')

## Define output file
outfile <- here('results/gphocs/output/gphocs_mergedlogs.txt')

## Define populations
poplist <- list(murW = 'A', B = 'A',
                murC = 'B', murE = 'B',
                griC = 'C', griW = 'C',
                A = 'root', C = 'root')

## Read metadata
lookup <- read.delim(infile_lookup, header = TRUE, as.is = TRUE)


# PROCESS GPHOCS OUTPUT LOGS ---------------------------------------------------

logs <- getlogs(
  setID = setID, logdir = logdir, poplist = poplist, lookup = lookup,
  burnin = burnin, last_sample = last_sample, subsample = subsample,
  mutrate_var = mutrate_var, gentime_sd = gentime_sd,
  rename_pops_before = rename_pops_before, rename_pops_after = rename_pops_after
)

## Write to file
write.table(logs, outfile,
            sep = '\t', quote = FALSE, row.names = FALSE)
