# SET-UP -----------------------------------------------------------------------

## Packages
library(here)

## Other scripts
source(here('scripts/gphocs-bpp/bpp_processlogs_fun.R'))

## Define settings
set.seed(7643)
setID <- 'mur3gri2c'
runID <- 'bpp_6mig'

mutrate_gen <- 1.52                  # Mutation rate per generation (will be multiplied by 1e-8)
mutrate_var <- 0.000001              # Variance for mutation rate  (will be multiplied by 1e-8)
gentime <- 3.5                       # Generation time
gentime_sd <- 1.16                   # SD for generation time

t_scale <- 1                         # Scaling of theta and tau
burnin <- 100000                     # Size of burn-in, i.e. iterations to remove
subsample_by <- 50                   # Subsample 1 in x samples (output lines). Default: 50

## Define input files
logdir <- here('results/bpp/output/raw/')
bpp_prefix <- 'bpp-mcmc_'
infile_poplookup <-  here('results/bpp/metadata/bpp_poplookup.txt')

## Define output file with combined BPP logs
outfile_log <- here('results/bpp/output/bpp_mergedlogs.txt')

## Read metadata
lookup <- read.table(infile_poplookup, header = TRUE)


# PROCESS BPP LOGS -------------------------------------------------------------

logs <- preplogs(
  logdir = logdir, bpp_prefix = bpp_prefix,
  setID = setID, runID = runID,
  lookup = lookup, t_scale = t_scale,
  burnin = burnin, subsample_by = subsample_by,
  gentime = gentime, gentime_sd = gentime_sd,
  mutrate_gen = mutrate_gen, mutrate_var = mutrate_var
  )

## Make sure all phi have the same direction and are close to 0, not 1 (invert if necessary):
logs <- logs %>%
  mutate(val = ifelse(var == 'phi' & val > 0.5, 1 - val, val))

## Write the combined BPP output logs to file
write.table(logs, outfile_log,
            row.names = FALSE, quote = FALSE, sep = '\t')
