#!/bin/bash

## Settings
VCF_DIR=results/geno/vcf/gatk/
QC_DIR=results/qc/vcf/gatk/
RUN_BCF=TRUE
RUN_BCF_BY_IND=FALSE
MEM=4

FILE_IDS=(r03.all.mac1.FS6 r03.all.mac1.FS7 r03.keepHybs.mac1.FS6)

## Run QC script
for file_id in "${FILE_IDS[@]}"; do

    sbatch --mem ${MEM}G -o slurm_qc-vcf_"$file_id" \
        scripts/qc/qc_vcf_run.sh "$file_id" "$VCF_DIR" "$QC_DIR" "$RUN_BCF" "$RUN_BCF_BY_IND"

done
