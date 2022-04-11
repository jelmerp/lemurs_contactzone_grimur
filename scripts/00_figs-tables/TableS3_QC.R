# SET-UP -----------------------------------------------------------------------
## Packages
library(here)
library(tidyverse)
library(janitor)
# NOTE: This script needs `bcftools` to be available in the $PATH

## Source script with some QC functions
source(here("scripts/qc/qc_fun.R"))

## VCF file IDs
ID_final <- "r03.all.mac1.FS6"
ID_rescue <- "r03.keepHybs.mac1.FS6"
ID_FS7 <- "r03.all.mac1.FS7"

## Define input files
infile_lookup <- here("metadata/hzlookup_bysample.txt")
infile_lookup_bylib <- here("metadata/hzlookup_bylib.txt")

dir_fq <- here("results/qc/fastq/")
dir_bam <- here("results/qc/bam/")
dir_vcf <- here("results/qc/vcf/gatk_joint/")

infile_dupd <- here(dir_fq, "dedupstats_sumr.txt")
infile_trim <- here(dir_fq, "trimstats_sumr.txt")

infile_proppairs <- here(dir_bam, "r03_proppairedreads.txt")
infile_bam_dp <- here(dir_bam, "r03_meanDP.txt")
infile_bamcolnames <- here(dir_bam, "eautils_colnames.txt")
infile_idepth <- here(dir_vcf, paste0(ID_final, ".idepth"))
infile_imiss <- here(dir_vcf, paste0(ID_final, ".imiss"))
infile_idepth_rescue <- here(dir_vcf, paste0(ID_rescue, ".idepth"))
infile_imiss_rescue <- here(dir_vcf, paste0(ID_rescue, ".imiss"))

infile_vcf_FS7 <- here("seqdata/vcf/gatk", paste0(ID_FS7, ".vcf.gz"))
inds_FS7 <- system(paste("bcftools query -l", infile_vcf_FS7), intern = TRUE)

## Define output files
outdir <- here("results/manuscript")
outfile_qc_ind <- here(outdir, "qc_byind_r03.txt")
outfile_qc_smr <- here(outdir, "qc_summary_r03.txt")
outfile_rescue <- here(outdir, "qc_rescue_r03.txt")
outfile_filterbycat <- here(outdir, "qc_filterbycat_r03.txt")

## Read metadata:
lookup <- read_tsv(infile_lookup) %>% select(ID, sp)
lookup_bylib <- read_tsv(infile_lookup_bylib)
longIDs <- lookup_bylib %>% filter(seqrun == "r03") %>% pull(ID)


# FASTQ STATS ------------------------------------------------------------------
## Stacks dedupping stats
dupd_raw <- readLines(infile_dupd)
longID <- gsub("slurm.deduptrim.(.*):.*", "\\1", dupd_raw)
fq_raw <- as.integer(gsub(".*:([0-9]+) pairs of reads input.*", "\\1", dupd_raw)) * 2
fq_dedupd <- as.integer(gsub(".* ([0-9]+) pairs of reads output.*", "\\1", dupd_raw)) * 2
dup_pct <- as.numeric(gsub(".* ([0-9]+.[0-9]+)% clone reads.*", "\\1", dupd_raw))
dupd_df <- data.frame(longID, fq_raw, fq_dedupd, dup_pct)

## FASTQ trimming stats
trim_raw <- readLines(infile_trim)
longID <- gsub("slurm.deduptrim.(.*):Input.*", "\\1", trim_raw)
fq_trimd <- as.integer(gsub(".*Both Surviving: ([0-9]+) .*", "\\1", trim_raw)) * 2 # *2 because these are nr of read pairs
trim_df <- data.frame(longID, fq_trimd)

## Stats on trimmed and dedupped fastqs:
fqstats_R1 <- read_fqstats_all(longIDs, read = "R1",
                               dir_fq = here(dir_fq, "by_ind/"))
fqstats_R2 <- read_fqstats_all(longIDs, read = "R2",
                               dir_fq = here(dir_fq, "by_ind/"))

fqstats_part <- rbind(fqstats_R1, fqstats_R2) %>%
  rename(longID = ID) %>%
  group_by(longID) %>%
  summarise(fq_qual = round(mean(fq_qual), 4),
            fq_len = round(mean(fq_meanlen), 2))

fqstats <- dupd_df %>%
  merge(., trim_df, by = "longID") %>%
  merge(., fqstats_part, by = "longID")


# BAM: NR OF PROP PAIRED READS, DEPTH, MQ --------------------------------------
## Samtools flagstats
flagstats <- read_flagstats_all(longIDs, bamdir = here(dir_bam, "by_ind/")) %>%
  select(ID, bam_map, bam_pp)

## Mean depth:
bam_dp <- read.delim(infile_bam_dp,
                     header = FALSE, as.is = TRUE, sep = " ",
                     col.names = c("ID", "bam_dp")) %>%
  filter(grepl("_", ID)) %>%
  mutate(bam_dp = round(bam_dp, 2))

## eautils-stats - mapping quality, etc:
eaustats_colnames <- readLines(infile_bamcolnames) %>%
  gsub(" ", "_", .) %>%  gsub("%", "pct", .)

eaustats_unfilt <-
  read_bamstats_all(longIDs,
                    mycolnames = eaustats_colnames,
                    suffix = ".rawbam.ea-utils-samstats.txt",
                    bamdir = paste0(dir_bam, "by_ind/")) %>%
  select(ID, mapq_mean) %>%
  rename(bam_mapq_unfilt = mapq_mean)

eaustats_filt <-
  read_bamstats_all(longIDs,
                    mycolnames = eaustats_colnames,
                    suffix = ".ea-utils-samstats.txt",
                    bamdir = paste0(dir_bam, "by_ind/")) %>%
  select(ID, mapq_mean) %>%
  rename(bam_mapq_filt = mapq_mean)

eaustats <- merge(eaustats_unfilt, eaustats_filt, by = "ID")

## Merge dataframes
bamstats <- merge(flagstats, bam_dp, by = "ID", all.x = TRUE) %>%
  merge(., eaustats, by = "ID", all.x = TRUE) %>%
  rename(longID = ID)


# VCF --------------------------------------------------------------------------
## idepth for joint-VCFs
idepth <- read.delim(infile_idepth, as.is = TRUE,
             colClasses = c("character", "integer", "numeric"),
             col.names = c("ID", "vcf_nSNP", "vcf_dp")) %>%
  select(ID, vcf_dp) %>%
  mutate(vcf_dp = round(vcf_dp, 2))

## imiss for joint-VCFs
imiss <-
  read.delim(infile_imiss, as.is = TRUE,
             colClasses = c("character", "integer", "integer", "integer", "numeric"),
             col.names = c("ID", "vcf_nSNP_all", "vcf_nfilt",
                           "vcf_nmiss", "vcf_fmiss")) %>%
  mutate(vcf_nSNP = vcf_nSNP_all - vcf_nmiss,
         vcf_miss_pct = (vcf_fmiss) * 100) %>%
  select(ID, vcf_nSNP, vcf_miss_pct)

vcfstats <- merge(idepth, imiss, by = "ID", all.x = TRUE)


# VCF STATS FOR "RESCUED" INDIVS -----------------------------------------------
idepth_rescue <-
  read.delim(infile_idepth_rescue, as.is = TRUE,
             colClasses = c("character", "integer", "numeric"),
             col.names = c("ID", "vcf_nSNP", "vcf_dp"))

inds_rescue <- idepth_rescue$ID[! idepth_rescue$ID %in% vcfstats$ID]

imiss_rescue <-
  read.delim(infile_imiss_rescue, as.is = TRUE,
             colClasses = c("character", "integer", "integer", "integer", "numeric"),
             col.names = c("ID", "vcf_nSNP_all", "vcf_nfilt", "vcf_nmiss", "vcf_fmiss"))

qc_rescue <- merge(idepth_rescue, imiss_rescue) %>%
  mutate(vcf_fmiss = round(vcf_fmiss, 3)) %>%
  select(ID, vcf_nSNP, vcf_dp, vcf_nSNP_all, vcf_fmiss) %>%
  filter(ID %in% inds_rescue) # rescue inds only


# MERGE DFs --------------------------------------------------------------------
## Merge FASTQ and BAM stats
fq_bam <-
  merge(fqstats, bamstats, by = "longID") %>%
  mutate(ID = substr(longID, 1, 7),
         bam_map_pct = round((bam_map / fq_trimd) * 100, 4),
         bam_pp_pct = round((bam_pp / bam_map) * 100, 4),
         repd = ifelse(duplicated(ID) | duplicated(ID, fromLast = TRUE), TRUE, FALSE))

## Remove 2nd library for mmur045 and mmur052
fq_bam <- fq_bam %>%
  filter(longID %in% c("mmur045_r03_p5g12", "mmur052_r03_p5h12"))

## Merge with VCF stats
qc <-
  merge(fq_bam, vcfstats, by = "ID", all.x = TRUE) %>%
  merge(., lookup, by = "ID") %>%
  mutate(
    filter_pass = ifelse(!is.na(vcf_dp), "FS6",
                         ifelse(ID %in% inds_FS7, "FS7",
                                ifelse(ID %in% inds_rescue, "rescue", "fail"))),
    sp_org = ifelse(ID == "mhyb004", "mgri", substr(ID, 1, 4)),
    score_fq = (fq_trimd * fq_qual) / 1000000,
    score_vcf = vcf_dp * vcf_nSNP,
    score_vcfbyfq = score_vcf / score_fq
    ) %>%
  select("ID", "longID", "sp_org", "sp", "repd", "fq_raw", "fq_dedupd", "dup_pct",
         "fq_trimd", "fq_qual", "fq_len",
         "bam_map", "bam_pp", "bam_dp", "bam_mapq_unfilt", "bam_mapq_filt",
         "bam_map_pct", "bam_pp_pct",
         "vcf_nSNP", "vcf_miss_pct", "vcf_dp",
         "filter_pass")

# SUMMARIZE --------------------------------------------------------------------
qc_smr_raw <- qc %>%
  filter(filter_pass == "FS6") %>%
  select(-ID, -longID, -sp, -filter_pass, -repd) %>%
  group_by(sp_org) %>%
  summarize_all(list(mean = mean, median = median))

qc_smr <- as_tibble(t(qc_smr_raw), rownames = "row_names") %>%
  row_to_names(row_number = 1) %>%
  mutate_at(c("mgri", "mhyb", "mmur"), as.numeric) %>%
  rename(stat = sp_org) %>%
  filter(!grepl("score", stat)) %>%
  arrange(stat)

## Pass by group:
filterbycat <- addmargins(table(qc$sp_org, qc$filter_pass))


# WRITE FILES ------------------------------------------------------------------
write_tsv(qc, outfile_qc_ind)
write_tsv(qc_smr, outfile_qc_smr)
write_tsv(qc_rescue, outfile_rescue)
write.table(filterbycat, outfile_filterbycat,
            sep = "\t", quote = FALSE, row.names = TRUE)
