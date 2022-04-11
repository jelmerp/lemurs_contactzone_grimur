## Packages
library(here)
library(tidyverse)

## Output file
outfile <- here("results/manuscript/TableS5_FST.txt")

## Input dir
indir <- here("results/fst/output")

## Functions
get_fst <- function(file) {
  as.numeric(sub(".*: ", "", grep("weighted", readLines(file), value = TRUE)))
}

## Contact zone inds (r03 VCF)
file_id <- "r03.all.mac1.FS6"
indir_full <- here(indir, file_id)
fst_files <- list.files(indir_full, pattern = "log$", full.names = TRUE)

pop1 <- gsub("(.*)_vs.*", "\\1", basename(fst_files))
pop2 <- gsub(".*_vs_(.*)_win.*", "\\1", basename(fst_files))
FST <- sapply(fst_files, get_fst)

fst_r03 <- data.frame(pop1, pop2, FST, VCF = "contact_zone",
                      row.names = NULL)

## All inds (hzproj1 VCF)
file_id <- "hzproj1.mac1.FS6"
indir_full <- here(indir, file_id)
fst_files <- list.files(indir_full, pattern = "log$", full.names = TRUE)

pop1 <- gsub("(.*)_vs.*", "\\1", basename(fst_files))
pop2 <- gsub(".*_vs_(.*)_win.*", "\\1", basename(fst_files))
FST <- sapply(fst_files, get_fst)

fst_hzproj1 <- data.frame(pop1, pop2, FST, VCF = "all_inds",
                          row.names = NULL)

## Final file with FST results from both VCF files
fst_all <- rbind(fst_r03, fst_hzproj1) %>%
  arrange(pop1, pop2)

write.table(fst_all, outfile,
            sep = "\t", quote = FALSE, row.names = FALSE)
