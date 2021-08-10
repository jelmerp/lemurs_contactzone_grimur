#!/bin/bash

# PREP ------------------------------------------------------------------------

## Combine the "failed" and 2nd run:
R1_A=r03_20190129/original/by_run/POELSTRA_5372_190116A1/raw/lemurRadseqHybridZone_r03_S25_L002_R1_001.fastq
R1_B=r03_20190129/original/by_run//POELSTRA_5372_190125B1/raw/lemurRadseqHybridZone_r03_S25_L002_R1_001.fastq
R1_COMBINED=r03_20190129/original/murgriHZ.r03_S25_L002_R1_001.fastq
cat $R1_A $R1_B >$R1_COMBINED

R2_A=r03_20190129/original/by_run/POELSTRA_5372_190116A1/raw/lemurRadseqHybridZone_r03_S25_L002_R2_001.fastq
R2_B=r03_20190129/original/by_run/POELSTRA_5372_190125B1/raw/lemurRadseqHybridZone_r03_S25_L002_R2_001.fastq
R2_COMBINED=r03_20190129/original/murgriHZ.r03_S25_L002_R2_001.fastq
cat $R2_A $R2_B >$R2_COMBINED

# RUN SCRIPT TO PROCESS FASTQS -------------------------------------------------

LIBRARY_ID="murgriHZ.r03"
INDIR=r03/original/
OUTDIR=r03/processed/
BARCODES=r03/metadata/barcodes_df_r03.txt #Two-column file with barcodes and corresponding sample IDs
SKIP=""
sbatch scripts/geno/1_fastq_process/01_fastq-process_runner.sh \
    -i $INDIR -o $OUTDIR -l $LIBRARY_ID -b $BARCODES -s "$SKIP"

## Process QC stats:
grep "Input Read Pairs" slurm.dedup* >r03/processed/QC/trimstats_sumr.txt
grep "pairs of reads input" slurm.dedup* >r03/processed/QC/dedupstats_sumr.txt
