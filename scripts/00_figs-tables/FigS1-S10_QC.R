# SET-UP -----------------------------------------------------------------------

## Source script with functions
source(here('scripts/genotyping/4_QC/qc_fun.R'))

## Load packages
library(here)
library(tidyverse)
library(janitor)

## Define input files
infile_qc <- here('results/qc/qc_byind_r03.txt')
infile_popcols <- here('metadata/colors.txt')

## Read and prep input files
popcols <- read.delim(infile_popcols, as.is = TRUE) %>%
   filter(pop %in% c('mmur', 'mgri'))

qc <- read.delim(infile_qc, as.is = TRUE) %>%
   filter(repd == FALSE) %>%
   mutate(
      fq_raw = fq_raw / 1000000,
      fq_dedupd = fq_dedupd / 1000000,
      fq_trimd = fq_trimd / 1000000,
      bam_map = bam_map / 1000000,
      bam_pp = bam_pp / 1000000,
      trim_pct = (1 - (fq_trimd / fq_dedupd)) * 100,
      filter_pass = factor(filter_pass, levels = c('FS6', 'FS7', 'rescue', 'fail'))
   )
colnames(qc)


# PLOT PREP --------------------------------------------------------------------

## Function with useful defaults
gg_hzqc <- function(
   my_y, my_df = qc, my_x = 'sp_org',
   colby = 'sp', colby_cols = c('#D69B12', '#FF00FF'),
   colby_labs = c(expression(italic("griseorufus")), expression(italic("murinus"))),
   colby_name = 'RADseq\nassignment:',
   shapeby = 'filter_pass', shapeby_name = 'filter status:',
   figname = NULL,
   ...
   ) {

   ggbox(my_df = my_df, my_x = my_x, my_y = my_y,
         colby = colby, colby_cols = colby_cols,
         colby_labs = colby_labs, colby_name = colby_name,
         shapeby = shapeby, shapeby_name = shapeby_name,
         figname = figname,
         ...)

}

## Set default labels
xlabs <- c(expression(italic("griseorufus")), 'hybrid', expression(italic("murinus")))
colby_labs <- c(expression(italic("griseorufus")), expression(italic("murinus")))
xtitle <- 'microsatellite assignment'


# CREATE PLOTS -----------------------------------------------------------------

## Plot FASTQ stats
gg_hzqc(my_y = 'fq_raw', xlabs = xlabs, xtitle = xtitle,
        ytitle = 'reads in unfiltered FASTQ (in millions)', openfig = TRUE)
gg_hzqc(my_y = 'dup_pct', xlabs = xlabs, xtitle = xtitle,
        ytitle = '% duplicate reads')
gg_hzqc(my_y = 'trim_pct', xlabs = xlabs, xtitle = xtitle,
        ytitle = '% reads removed during FASTQ QC')

## Plot BAM stats
gg_hzqc(my_y = 'bam_map_pct', xlabs = xlabs, xtitle = xtitle,
        ytitle = '% successfully mapped reads')
gg_hzqc(my_y = 'bam_pp_pct', xlabs = xlabs, xtitle = xtitle,
        ytitle = '% properly paired reads in BAM')
gg_hzqc(my_y = 'bam_pp', xlabs = xlabs, xtitle = xtitle,
        ytitle = 'properly paired reads in filtered BAM (in millions)')
gg_hzqc(my_y = 'bam_mapq_unfilt', xlabs = xlabs, xtitle = xtitle,
        ytitle = 'mean mapping quality in unfiltered BAM', ymax = 48)
gg_hzqc(my_y = 'bam_mapq_filt', xlabs = xlabs, xtitle = xtitle,
        ytitle = 'mean mapping quality in filtered BAM')
gg_hzqc(my_y = 'bam_dp', xlabs = xlabs, xtitle = xtitle,
        ytitle = 'mean depth in filtered BAM')

## Plot VCF stats
gg_hzqc(my_y = 'vcf_dp', xlabs = xlabs, xtitle = xtitle,
        shapeby = NULL, ytitle = 'Mean FS6 VCF depth')
gg_hzqc(my_y = 'vcf_nSNP', xlabs = xlabs, xtitle = xtitle,
        shapeby = NULL, ytitle = 'nr of passed SNPs in FS6 VCF')
gg_hzqc(my_y = 'vcf_miss_pct', xlabs = xlabs, xtitle = xtitle,
        shapeby = NULL, ytitle = '% missing data in FS6 VCF')

## Plot filter status
gg_hzqc(my_x = 'filter_pass', my_y = 'bam_pp', shapeby = NULL,
        xtitle = "filtering status",
        ytitle = 'properly paired reads in filtered BAM (in millions)',
        figfile = 'filterstatus')
