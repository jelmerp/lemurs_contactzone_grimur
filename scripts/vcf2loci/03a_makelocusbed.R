#!/usr/bin/env Rscript

#### SET-UP --------------------------------------------------------------------
cat("## 03a_makelocusbed.R: Starting script.\n")

## Command-line args:
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
set_id <- args[1]
infile_inds <- args[2]
indir_bed <- args[3]
outfile_bed <- args[4]

## Hardcoded settings:
min_element_ovl <- 0.9
min_element_ovl_trim <- 0.8
min_locus_size <- 100
max_dist_within_ind <- 10
max_dist_between_ind <- 0
min_element_size <- 25

## Packages:
if (!"pacman" %in% rownames(installed.packages())) install.packages("pacman")
library(pacman)
packages <- c("data.table", "tidyverse", "valr", "IRanges")
p_load(char = packages, install = TRUE)

## Process command-line args:
ids <- readLines(infile_inds)
bedfiles <- paste0(indir_bed, "/", ids, ".callable.bed")

min_element_ovl <- length(bedfiles) * min_element_ovl
min_element_ovl_trim <- length(bedfiles) * min_element_ovl_trim

## Report:
cat("\n## Set ID:", set_id, "\n")
cat("## Bedfile dir:", indir_bed, "\n")
cat("## Bedfiles:", bedfiles, "\n")
cat("## Number of bedfiles:", length(bedfiles), "\n")
cat("## Bed output file:", outfile_bed, "\n\n")

cat("## Max distance within inds:", maxdist_within_ind, "\n")
cat("## Max distance between inds:", max_dist_between_ind, "\n")
cat("## Min element overlap - for locus creation:", min_element_ovl, "\n")
cat("## Min element overlap - for locus trimming:", min_element_ovl_trim, "\n")
cat("## Min element size:", min_element_size, "\n")
cat("## Min locus size:", min_locus_size, "\n")

#### FUNCTIONS -----------------------------------------------------------------
get_indloci <- function(bedfile, max_dist, last_row = 0) {
    cat("## get_indloci function: bedfile:", bedfile, "\n")

    bed <- fread(bedfile,
        header = FALSE,
        colClasses = c("character", "integer", "integer", "character"),
        col.names = c("chrom", "start", "end", "status")
    )
    bed$status <- NULL

    if (last_row == 0) {
        cat("## get_indloci(): Using all rows...\n")
    } else {
        cat(
            "## get_indloci(): Selecting rows until row:",
            last_row, "\n"
        )
        bed <- bed[1:last_row, ]
    }

    bed <- as.tbl_interval(bed)
    bed <- bed_merge(bed, max_dist = max_dist)
    bed <- arrange(bed, start)
    return(bed)
}

collect_indloci <- function(bedfiles, last_row,
                            maxdist_within_ind, min_element_size) {
    bedlist <- lapply(bedfiles, get_indloci,
        max_dist = maxdist_within_ind, last_row = last_row
    )
    bed_by_ind <- do.call(rbind, bedlist)

    # Somehow this step is necessary or bed_merge wont work:
    bed_by_ind <- data.frame(
        chrom = bed_by_ind$chrom,
        start = bed_by_ind$start,
        end = bed_by_ind$end
    ) %>%
        as.tbl_interval() %>%
        arrange(start) %>%
        filter(end - start >= min_element_size)

    return(bed_by_ind)
}

trim_locus <- function(row.nr, locus_df, element_df, min_element_ovl_trim) {
    locus <- locus_df[row.nr, ]
    locus_id <- paste0(locus$chrom, "_", locus$start)

    locus_elements <- bed_intersect(bed_by_ind, locus)
    bed_cov_ir <- IRanges(
        start = locus_elements$start.x,
        end = locus_elements$end.x
    )
    cov <- IRanges::coverage(bed_cov_ir)

    ok <- which(cov@values > min_element_ovl_trim)

    if (length(ok >= 1)) {
        first_ok <- ok[1]
        if (first_ok > 1) {
            first_base <- sum(cov@lengths[1:first_ok])
        } else {
            first_base <- locus$start
        }

        last_ok <- ok[length(ok)]
        if (last_ok < length(cov@values)) {
            last_base <- sum(cov@lengths[1:last_ok])
        } else {
            last_base <- locus$end
        }

        locus_length <- locus$end - locus$start
        trimmed_start <- first_base - locus$start
        trimmed_end <- locus$end - last_base
        locus_length_final <- last_base - first_base

        cat(
            row.nr, locus_id, "/ length:", locus_length,
            "/ trimmed start:", trimmed_start, "/ trimmed end:", trimmed_end,
            "/ remaining bases:", locus_length_final, "\n"
        )

        locus_trimmed <- data.frame(
            chrom = locus$chrom,
            start = first_base, end = last_base
        )
        locus_trimmed_length <- locus_trimmed$end - locus_trimmed$start

        if (locus_trimmed_length < min_locus_size) {
            cat(row.nr, "Locus too small...\n")
        } else {
            return(locus_trimmed)
        }
    } else {
        cat(row.nr, "Coverage too low: skipping locus...\n")
    }
}


#### RUN -----------------------------------------------------------------------
## Get per-individual bed file and merge loci:
bed_by_ind <- collect_indloci(bedfiles,
    last_row = last_row,
    maxdist_within_ind = maxdist_within_ind,
    min_element_size = min_element_size
)

## Merge per-individual loci:
bed_merged <- bed_merge(bed_by_ind, max_dist = max_dist_between_ind)
cat("\n## Nr of initial loci:", nrow(bed_merged), "\n")

## Calculate "coverage" (number of elements) overlapping with each locus:
## Retain only those with "min.elements" number of overlapping elements,
## and "min.frac" fraction of overlap (latter is not very important)
bed_cov <- bed_coverage(bed_merged, bed_by_ind) %>%
    filter(.ints >= min_element_ovl)
cat(
    "## Number of loci after filtering by coverage:",
    nrow(bed_cov), "\n\n"
)

## Trim loci:
cat("## Trimming loci...\n")
bed_trim_list <- lapply(1:nrow(bed_cov), trim_locus,
    locus_df = bed_cov,
    element_df = bed_by_ind,
    min_element_ovl_trim = min_element_ovl_trim
)
bed_trim <- do.call(rbind, bed_trim_list)
cat("## Trimming done.\n\n")
cat("## Number of loci retained:", nrow(bed.trim), "\n\n")

## Save bedfile:
cat("## Writing bedfile", outfile_bed, "\n")
write.table(bed.trim, outfile_bed,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

## Report:
cat("\n## Done with script.\n")