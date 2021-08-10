#!/bin/bash
set -euo pipefail

## Setup:
locusbed_intermed=$1
locusbed_final=$2
vcf_highdepth=$3
vcf_filtered_intersect=$4

echo "## Starting script $0"

## Run bedtools:
bedtools intersect -v -a "$locusbed_intermed" -b "$vcf_highdepth" >"$locusbed_final"

## Report:
echo "## Nr loci after removing too-high-depth variants: $(wc -l <"$locusbed_final")"

nr_snps_in_loci=$(bedtools intersect -u -a "$vcf_filtered_intersect" -b "$locusbed_final" | grep -cv "##")
nr_snps_in_vcf=$(grep -cv "##" "$vcf_filtered_intersect")
nr_snps_lost=$(("$nr_snps_in_vcf" - "$nr_snps_in_loci"))

echo "## Nr lost SNPs (in VCF but not in loci): $nr_snps_lost"
echo "## Total nr SNPs in VCF:                  $nr_snps_in_vcf"
echo "## Total nr SNPs in loci:                 $nr_snps_in_loci"

echo "## Done with script $0"
