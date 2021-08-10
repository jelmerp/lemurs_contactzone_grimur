#!/usr/bin/env Rscript

## Script to prepare a "popfile" - population assignment file - for Treemix.

#### SET-UP --------------------------------------------------------------------
cat("\n######################################################################\n")
cat("#### Script: treemix_makePopfile.R\n")

## Packages:
if(!'pacman' %in% rownames(installed.packages())) install.packages('pacman')
library(pacman)
packages <- c('tidyverse')
p_load(char = packages, install = TRUE)

## Command-line args:
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
infile_focInds <- args[1]
infile_lookup <- args[2]
outfile_popfile <- args[3]
group.by.column <- args[4]

## Report:
cat('#### treemix_makePopfile.R: infile_focInds:', infile_focInds, '\n')
cat('#### treemix_makePopfile.R: infile_lookup:', infile_lookup, '\n')
cat('#### treemix_makePopfile.R: outfile_popfile:', outfile_popfile, '\n')
cat('#### treemix_makePopfile.R: lookup column name to group by:', group.by.column, '\n\n')

## Read input files:
# List of focal inds:
focInds <- readLines(infile_focInds)
cat('#### treemix_makePopfile.R: Number of focal inds:', length(focInds), '\n')

# metadate with population grouping:
inds.df <- read.delim(infile_lookup, as.is = TRUE) %>%
  filter(ID %in% focInds) %>%
  select(ID, group.by.column)
cat('#### treemix_makePopfile.R: Dimensions of metadata file:', dim(inds.df), '\n\n')

## Remove individuals with no population assignments:
notanys <- which(is.na(inds.df[, group.by.column]))
if(length(notanys) > 0) {
  cat('#### treemix_makePopfile.R: Removing',
      length(notanys), 'inds with no group assignment\n')
  inds.df <- inds.df[-notanys, ]
  cat('#### treemix_makePopfile.R: Dimensions of metadata file after filtering:',
      dim(inds.df), '\n')
}

################################################################################
#### CREATE TREEMIX POPFILE #####
################################################################################
getLine <- function(focalGroup) {
  inds.df.filt <- inds.df[which(inds.df[, group.by.column] == focalGroup), ]

  inds <- inds.df.filt %>%
    pull(ID) %>%
    paste0(collapse = ' ')

  return(paste0(focalGroup, ': ', inds))
}

focalGroups <- unique(inds.df[, group.by.column])
lines <- sapply(focalGroups, getLine)
names(lines) <- NULL
writeLines(lines, outfile_popfile)

cat('#### treemix_makePopfile.R: Done with script\n\n')
