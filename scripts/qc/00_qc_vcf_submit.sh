#!/bin/bash

## Settings
vcf_dir=results/geno/vcf
qc_dir=results/geno/qc/vcf
run_bcf=TRUE
run_bcf_ind=FALSE
file_ids=(r03.all.mac1.FS6 r03.all.mac1.FS7 r03.keepHybs.mac1.FS6)

## Run QC script
for id in "${file_ids[@]}"; do
    sbatch scripts/qc/qc_vcf_run.sh "$id" "$vcf_dir" "$qc_dir" "$run_bcf" "$run_bcf_ind"
done
