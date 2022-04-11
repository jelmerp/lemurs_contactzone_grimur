#!/bin/bash

# SETUP ------------------------------------------------------------------------
## Scripts
SCRIPT_GENO_INDS=scripts/geno/01_ind-geno_runner.sh
SCRIPT_GENO_JOINT=scripts/geno/02_joint-geno_runner.sh
SCRIPT_FILTER=scripts/03_filterVCF_FS6_runner.sh

## Input files
IDS_FILE_ALL=metadata/inds_r03.txt
IDS_FILE_HZPROJ1=metadata/inds_hzproj1.txt
META_FILE=metadata/hzlookup_bylibrary.txt
REF_DIR=data/refgenome
REF_ID=GCF_000165445.2_Mmur_3.0_genomic_stitched  # Generate "stitched" genome from genome on NCBI by running `scripts/conversion/stitchScaffolds_submit.sh`
REF=$REF_DIR/$REF_ID.fasta
SCAFFOLD_FILE=$REF_DIR/$REF_ID.scaffoldList.txt

## Settings and output file
ids_short=($(cat "$IDS_FILE_ALL")) # Replicates: mmur045 + mmur052
gatk_version=gatk4
fq_dir=data/fastq
bam_dir=results/bam
vcf_dir=results/geno
qc_dir_vcf=results/geno/qc/vcf_ind
qc_dir_bam=results/geno/qc/bam
minmapqual=30
mem=12
ncores=4
use_r2=TRUE
regions_file=notany
exclude_or_include_regions=notany


# FASTQ TO GVCF ----------------------------------------------------------------
bam_suffix=notany
skip_flags=""

for id_short in ${ids_short[@]}; do
    lane=$(grep "$id_short" "$META_FILE" | cut -f 8 | head -n 1)
    library=$(grep "$id_short" "$META_FILE" | cut -f 9 | head -n 1)
    readgroup_string="@RG\tID:${lane}\tSM:${id_short}\tPL:ILLUMINA\tLB:$library"

    id_long_file=metadata/r03/replicates/$id_short.txt
    grep "$id_short" "$META_FILE" | cut -f 1 > "$id_long_file"

    sbatch "$SCRIPT_GENO_INDS" \
        "$id_short" "$id_long_file" "$use_r2" "$REF" \
        "$fq_dir" "$bam_dir" "$vcf_dir" "$qc_dir_vcf" "$qc_dir_bam" \
        "$minmapqual" "$gatk_version" "$bam_suffix" "$readgroup_string" \
        "$regions_file" "$exclude_or_include_regions" "$mem" "$ncores" "$skip_flags"
done

## Process QC results
cat "$qc_dir_bam"/*meandepth* > "$qc_dir_bam"/r03_meanDepth.txt
cat "$qc_dir_vcf"/vcftools/*FS6*idepth | grep -v "INDV" > "$qc_dir_vcf"/vcftools/r03_idepth.txt


# JOINT GENOTYPING -------------------------------------------------------------
Gvcf_dir=results/geno/vcf/ind
vcf_dir=results/geno/vcf/joint
qc_dir_vcf=results/geno/qc/vcf
increment=1
start_at=1
stop_at="none"
add_commands="none"
mem_job=36
mem_gatk=28
ncores=1

## All r03 inds:
file_id=r03.all
ids_short=($(cat "$IDS_FILE_ALL"))
"$SCRIPT_GENO_JOINT" \
    "$file_id" "$SCAFFOLD_FILE" "$increment" "$start_at" "$stop_at" \
    "$Gvcf_dir" "$vcf_dir" "$qc_dir_vcf" "$REF" \
    "$add_commands" "$mem_job" "$mem_gatk" "$ncores" \
    ${ids_short[@]}

## Including allopatric pops and 3 outgroup rufus
file_id=hzproj1
ids_short=($(cat "$IDS_FILE_HZPROJ1"))
"$SCRIPT_GENO_JOINT" \
    "$file_id" "$SCAFFOLD_FILE" "$increment" "$start_at" "$stop_at" \
    "$Gvcf_dir" "$vcf_dir" "$qc_dir_vcf" "$REF" \
    "$add_commands" "$mem_job" "$mem_gatk" "$ncores" \
    ${ids_short[@]}


# GENERATE keepHybs VCF WITH MOST PUTATIVE HYBRIDS -----------------------------
file_id=r03.all
input_name=$file_id.rawSNPs.ABHet
indfile=metadata/indsel/r03.keepHybs.txt # To be generated below
output_name=r03.keepHybs
vcf_dir_main=results/geno/vcf/joint/intermed
vcf_dir_final=results/geno/vcf/joint/final
qc_dir=results/geno/qc/vcf_joint/
dp_mean=15
MAC=3
mem=12
skip_common_steps="-456789tew"
skip_final_steps="-123"
select_inds_by_file=TRUE
filter_inds_by_missing=FALSE
jobname=filterVCF.$output_name
skip_common_steps="-12456789tew"
skip_in_pip=""
cat "$IDS_FILE_ALL" |
    grep -Ev 'mgan015|mhyb005|mgri068|mgri073|mgri074|mgri081|mgri083|mgri089|mgri091|mgri092|mgri096|mgri101|mmur054|mmur062|mmur074|mmur075' >$indfile

$SCRIPT_FILTER \
    "$input_name" "$output_name" "$vcf_dir_main" "$vcf_dir_final" "$qc_dir" "$REF" \
    "$dp_mean" "$MAC" "$filter_inds_by_missing" "$select_inds_by_file" "$indfile" \
    "$mem" "$jobname" "$skip_common_steps" "$skip_final_steps" "$skip_in_pip"
