#!/usr/bin/env Rscript

#### SET-UP --------------------------------------------------------------------
cat("\n----------------------\n")
cat("## Starting with script admixtools_makeIndfile.R\n")

library(tidyverse)
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)

inds.focal.file <- args[1]
inds.metadata.file <- args[2]
indfile.out <- args[3]
IDcolumn <- args[4]
groupby <- args[5]

cat('## inds.focal.file:', inds.focal.file, '\n')
cat('## inds.metadata.file:', inds.metadata.file, '\n')
cat('## indfile.out:', indfile.out, '\n')
cat('## ID column:', IDcolumn, '\n')
cat('## groupby column:', groupby, '\n')


#### CREATE INPUT FILE ---------------------------------------------------------
## Read files:
inds.focal <- readLines(inds.focal.file)

inds.df <- read.delim(inds.metadata.file, as.is = TRUE)
inds.df <- inds.df[inds.df[, IDcolumn] %in% inds.focal, ]
inds.df <- inds.df[, c(IDcolumn, groupby)]
print(inds.df)
stopifnot(all(inds.focal %in% inds.df[, IDcolumn]))

## Create input file:
getLine <- function(ind.ID) {
  group <- inds.df[, groupby][inds.df[, IDcolumn] == ind.ID]
  line <- data.frame(ind.ID, U = 'U', group)
  return(line)
}

indfile.df <- do.call(rbind, lapply(inds.focal, getLine))
indfile.df <- indfile.df %>% arrange(group)


#### REPORT AND WRITE FILE -----------------------------------------------------
cat('## Resulting indfile.df:\n')
print(indfile.df)

write.table(indfile.df, indfile.out,
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

cat('\n\n## Done with script.\n')
