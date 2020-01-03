################################################################################
#### SET-UP #####
################################################################################
## Libraries and scripts:
library(tidyverse)
library(here)

script_plotfun <- here('scripts/generalscripts_lemurs_link/geno/vcf2loci/locusstats_plotfun.R')
source(script_plotfun)

## Parameters:
fileID <- 'mur3gri2c'

## Input files:
infile_stats <- paste0(here(), '/analyses/vcf2fullfasta/', fileID, '.locusstats.txt')
infile_allscafs <- here('seqdata/seqdata_lemurs_link/reference/mmur/scaffolds.txt')

## Output files:
dir_out <- here('analyses/geno/vcf2loci/')
dir_plot <- file.path(dir_out, 'figs')
outfile_xloci <- paste0(dir_out, '/', fileID, '.sexlinkedloci.txt')
outfile_stats <- paste0(dir_out, '/', fileID, '.statsprocessed.txt')
if(!dir.exists(dir_plot)) dir.create(dir_plot, recursive = TRUE)

## Scaffolds - metadata:
allscafs <- readLines(infile_allscafs)
scaf <- allscafs[grepl('NC_', allscafs)]


################################################################################
#### PREP STATS DF #####
################################################################################
## Read and process stats file:
stats <- read.delim(infile_stats, as.is = TRUE) %>%
  select(1:12)
colnames(stats) <-  c('locus.full', 'nInd', 'bp', 'nCells', 'nN', 'pN',
                      'nvar', 'pvar', 'nPars', 'pPars', 'AT', 'GC')
stats <- stats %>%
  mutate(locus = gsub('.fa', '', locus.full)) %>%
  separate(locus, into = c('scaffold', 'pos'), sep = ':') %>%
  separate(pos, into = c('start', 'end'), sep = '-') %>%
  mutate(start = as.integer(start), end = as.integer(end)) %>%
  select(locus.full, scaffold, start, end, nInd, bp, nN, pN, nvar, pvar, nPars,
         pPars, AT, GC) %>%
  #filter(scaffold != 'Super_Scaffold0') %>%
  arrange(scaffold, start) %>%
  group_by(scaffold) %>%
  mutate(distToNext = lead(start) - (start + bp), # Include distance-to-next locus
         calledSites = round((bp * nInd) * (100 - pN)))

write.table(stats, outfile_stats, sep = '\t', quote = FALSE, row.names = FALSE)


################################################################################
#### PRINT BASIC STATS #####
################################################################################
cat("Mean locus size:", round(mean(stats$bp)), 'bp \n')
cat("Locus size - quantiles:\n")
quantile(stats$bp)

cat("Mean nr of variable sites:", round(mean(stats$nvar), 2), '\n')
cat("Nr of variable sites - quantiles:\n")
quantile(stats$nvar)

cat("Mean nr of parsimony-informative sites:", round(mean(stats$nPars), 2), '\n')
cat("Nr of PI sites - quantiles:\n")
quantile(stats$nPars)

cat("Mean % missing data:", round(mean(stats$pN), 2), '% \n')
cat("Nr of PI sites - quantiles:\n")
round(quantile(stats$pN), 2)

cat("Mean distance to next locus:", round(mean(stats$distToNext, na.rm = TRUE)), '\n')
cat("Distance to next locus - quantiles:\n")
round(quantile(stats$distToNext, na.rm = TRUE), 2)
cat("Number of loci with a neighbouring locus at <10 kb:\n")
length(which(stats$distToNext < 10000))
cat("Number of loci with a neighbouring locus at <5 kb:\n")
length(which(stats$distToNext < 5000))
cat("Number of loci with a neighbouring locus at <1 kb:\n")
length(which(stats$distToNext < 1000))

cat("Number of invariant loci:", length(which(stats$nvar == 0)), '\n')

## By scaffold:
#table(stats$scaffold)
cat("Chromosome-level scaffolds that are not present among the loci:\n")
scaf[!scaf %in% stats$scaffold]

# stats %>% group_by(scaffold) %>%
#   summarise(pvar_mean = mean(pvar), pPars_mean = mean(nPars)) %>%
#   arrange(pPars_mean) %>%
#   print(n=50)

xloci <- stats %>% filter(scaffold == 'NC_033692.1') %>% pull(locus.full)
cat("Number of X chromosome (NC_033692.1) loci:", length(xloci))
writeLines(xloci, outfile_xloci)


################################################################################
#### CREATE PLOTS #####
################################################################################
oneVarPlot('bp', xtitle = 'Locus length (bp)', xmax = 1000)
oneVarPlot('nvar', xtitle = 'Number of variable sites per locus', nbins = max(stats$nvar))
oneVarPlot('nPars', xtitle = 'Number of PI sites per locus', nbins = max(stats$nPars))
oneVarPlot('pN', xtitle = 'Percentage Ns per locus')

stats$distToNext <- stats$distToNext/1000
oneVarPlot('distToNext', xtitle = 'Distance to next locus (kbp)', xmax = 1000)

byPlot('scaffold', 'bp', ytitle = 'Locus length (bp)')
byPlot('scaffold', 'nvar', ytitle = 'Number of variable sites')
byPlot('scaffold', 'nPars', ytitle = 'Number of PI sites')
byPlot('scaffold', 'pN', ytitle = 'Percentage missing data')
byPlot('scaffold', 'distToNext', ytitle = 'Distance to next locus (bp)', ymax = 50000)
