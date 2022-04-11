#!/bin/bash

## Software
conda activate bcftools-env

## Settings
VCF_DIR=results/gatk/                   # Dir with existing VCF file(s)
PLINK_DIR=results/admixture/input       # Dir for PLINK files (to be produced)
ADMIX_INDIR=results/admixture/input
ADMIX_OUTDIR=results/admixture/output/  # Dir for ADMIXTURE output files (to be produced)
FILE_IDS=(r03.all.mac1.FS6 r03.all.mac1.FS7 r03.keepHybs.mac1.FS6) # File IDs for VCF files, VCF files should be: $VCF_DIR/$FILE_ID.vcf.gz

mkdir -p "$ADMIX_INDIR"

## Run Admixture
for file_id in "${FILE_IDS[@]}"; do
    echo -e "\n## VCF file ID: $file_id"

    indfile="$ADMIX_INDIR"/"$file_id"_indlist.txt                               # List of individuals to include in ADMIXTURE analysis
    bcftools query -l "$VCF_DIR"/"$file_id".vcf.gz | grep -v "mruf" >"$indfile" # Send list of inds minus "mruf" (outgroup) to separate textfile

    sbatch scripts/admixture/01_admixture_runner.sh "$file_id" "$VCF_DIR" \
        "$PLINK_DIR" "$ADMIX_OUTDIR" -i "$indfile" -k 1 -K 9
done
