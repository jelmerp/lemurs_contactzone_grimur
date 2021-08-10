#!/bin/bash
set -euo pipefail

# Software
# Needs R 3.4+

## Scripts:
SCRIPT_A=03a_makelocusbed.R
SCRIPT_B=03b_intersect.sh
SCRIPT_C=03c_mergedfasta.sh
SCRIPT_D=03d_locusfasta.sh
SCRIPT_F=03e_locusstats.sh
SCRIPT_E=03f_filterloci.sh

## Command-line args:
set_id=$1
file_inds=$2
dir_bed=$3
dir_fasta=$4
vcf_filtered_intersect=$5
vcf_highdepth=$6

## Process args:
jobname=vcf2loci_$set_id

dir_indfasta=$dir_fasta/"$set_id"_byind/
dir_locusfasta_intermed=$dir_fasta/"$set_id"_bylocus_intermed/
dir_locusfasta_final=$dir_fasta/"$set_id"_bylocus_final/

locusbed_intermed=$dir_bed/"$set_id"_loci_intermed.bed
locusbed_final=$dir_bed/"$set_id"_loci_all.bed
locuslist=$dir_bed/locuslist.txt
fasta_merged=$dir_indfasta/"$set_id"_merged.fasta

locusstats_intermed=$dir_bed/"$set_id"/"$set_id"_locusstats_all.txt
locusstats_final=$dir_bed/"$set_id"/"$set_id"_locusstats_filt.txt

[[ ! -d $dir_locusfasta_intermed ]] && mkdir "$dir_locusfasta_intermed"
[[ ! -d $dir_locusfasta_final ]] && mkdir "$dir_locusfasta_final"
[[ ! -d $dir_bed/"$set_id" ]] && mkdir "$dir_bed"/"$set_id"

## Report:
date
echo "##  Starting script $0"
echo "##  Set ID:                              $set_id"
echo "##  Dir for bedfiles (etc):              $dir_bed"
echo "##  Dir for fasta files:                 $dir_fasta"
echo "##  [input] file with individual ids:    $file_inds"
echo "##  [input] vcf - filtered:              $vcf_filtered_intersect"
echo "##  [input] vcf - excessive coverage:    $vcf_highdepth"
echo "##  [output] bedfile w/ loci - intermed: $locusbed_intermed"
echo "##  [output] bedfile w/ loci - final:    $locusbed_final"

echo "##  contents of file id list:" && cat "$file_inds"

## A. Create locus-BEDfile:
Rscript "$SCRIPT_A" "$set_id" "$file_inds" "$dir_bed" "$locusbed_intermed"

## B. Intersect locus-bedfile with high-depth VCF:
"$SCRIPT_B" "$locusbed_intermed" "$locusbed_final" "$vcf_highdepth" "$vcf_filtered_intersect"

## C. Get fasta with all loci and inds:
"$SCRIPT_C" "$file_inds" "$dir_indfasta" "$locusbed_final" "$locuslist" "$fasta_merged"

## D. Create by-locus FASTA files:
cluster_command="sbatch --job-name=$jobname -o slurm_vcf2loci_03d.$set_id"
$cluster_command "$SCRIPT_D" "$locuslist" "$dir_locusfasta_intermed" "$fasta_merged"

## R. Compute locus-stats for intermediate loci:
cluster_command="sbatch --job-name=$jobname --dependency=singleton -o slurm_vcf2loci-03e_$set_id"
$cluster_command "$SCRIPT_E" "$dir_locusfasta_intermed" "$locusstats_intermed"

## F. Filter loci:
cluster_command="sbatch --job-name=$jobname --dependency=singleton -o slurm_vcf2loci_03e_$set_id"
$cluster_command "$SCRIPT_F" "$locusstats_intermed" "$dir_locusfasta_intermed" "$dir_locusfasta_final"

## G. Compute locus-stats for final loci:
cluster_command="sbatch --job-name=$jobname --dependency=singleton -o slurm_vcf2loci-03f_$set_id"
$cluster_command "$SCRIPT_E" "$dir_locusfasta_final" "$locusstats_final"

# Report:
echo "## Done with script $0."
date
