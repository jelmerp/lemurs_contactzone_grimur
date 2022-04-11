#!/bin/bash
set -euo pipefail

## Software:
# bcftools
# vcftools

## Settings defined as constants:
# Minimum depth for GATK CallableLoci command: 3 (02_process_inds.sh)
# Maximum percentage of missing data:          10 (03f_filterloci.R)
# Minimum distance between loci:               10000 bp  (03f_filterloci.R)
# Minimum locus size:                          100 bp (03a_makelocusbed.R -- + other more specific settings for locus creation)

## Scripts:
SCRIPT01=01_maskbed.sh
SCRIPT02=02_vcf2loci1.sh
SCRIPT03=03_vcf2loci2.sh

## Positional args:
set_id=$1             # Arbitrary ID for set of loci
ref=$3                # INPUT: Reference genome FASTA file
dir_bam=$3            # INPUT: Directory with BAM files
suffix_bam=$4         # INPUT: Common suffix for BAM files
vcf_altref=$5         # INPUT: Raw (pre-filtering) VCF -- Used for producing altref fasta file
vcf_filtered_mask=$6  # INPUT: Partially filtered VCF -- Used for masking bad sites [Produce this by running the filter pipeline without the two filter-missing steps]
vcf_filt_intersect=$7 # INPUT: Final (fully filtered) VCF -- Used to select loci with only high-qual snps
vcf_highdepth=$8      # INPUT: VCF with sites with excessive depth [This is "$vcf_hidepth" in VCF filter script 05_filter_max-dp.sh]
outdir=$9             # Output directory

## Other variables:
JOBNAME=vcf2fasta.$set_id
dir_fasta=$outdir/fasta
dir_indfasta=$dir_fasta/$set_id.byind
dir_bed=$outdir/bed
bed_removed_sites=$dir_bed/"$set_id".sitesinvcfremovedbyfilters.bed

## Create dirs if necessary:
[[ ! -d "$dir_fasta" ]] && mkdir -p "$dir_fasta"
[[ ! -d "$dir_indfasta" ]] && mkdir -p "$dir_indfasta"
[[ ! -d "$dir_bed" ]] && mkdir -p "$dir_bed"

## get individuals present in vcf:
inds=($($bcftools query -l $vcf_filt_intersect))
file_inds=slurm.indfile.$set_id.tmp
printf "%s\n" "${inds[@]}" >"$file_inds"

## Script 1 - create bedfile with sites removed by vcf filtering
cluster_command="sbatch --job-name=$JOBNAME --dependency=singleton -o slurm_vcf2loci_01.$set_id"
maskbed_job=$($cluster_command "$SCRIPT01" "$vcf_altref" "$vcf_filtered_mask" "$bed_removed_sites")
maskbed_dependency="--dependency=afterok:${maskbed_job/submitted batch job //}"

## Script 2 - for each ind, create calllable-loci bed, and masked alt-ref fasta
for ind in "${inds[@]}"; do
    bam="$dir_bam/$ind.$suffix_bam*bam"
    cluster_command="sbatch --job-name=$JOBNAME $maskbed_dependency -o slurm_vcf2loci_02_${ind}_${set_id}"
    $cluster_command "$SCRIPT02" \
        "$ind" "$set_id" "$bam" "$vcf_altref" "$bed_removed_sites" "$dir_indfasta" "$dir_bed" "$ref"
done

# Script 3 - for all inds at once, delimit and extract loci
cluster_command="sbatch --job-name=$JOBNAME --dependency=singleton -o slurm_vcf2loci_03_$set_id"
$cluster_command "$SCRIPT03" \
    "$set_id" "$file_inds" "$dir_bed" "$dir_fasta" "$vcf_filt_intersect" "$vcf_highdepth"

echo "## Done with script $0"
