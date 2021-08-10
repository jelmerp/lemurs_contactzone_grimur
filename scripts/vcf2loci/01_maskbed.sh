#!/bin/bash
set -euo pipefail

## create a mask bedfile from vcf files

## Software:
# uses vcf2bed from bedops
# uses bedtools

## Command-line args:
vcf_altref=$1
vcf_filtered_mask=$2
bed_removed_sites=$3

## Report:
date
echo "## starting script maskbed.sh"
echo "## input: altref vcf: $vcf_altref"
echo "## input: filtered vcf (for mask): $vcf_filtered_mask"
echo "## output: bedfile with removed sites: $bed_removed_sites"

## Run vcf2bed and bedtools-intersect
bedtools intersect -v \
    -a <(vcf2bed --sort-tmpdir=tmpdir <"$vcf_altref" | cut -f 1,2,3) \
    -b <(vcf2bed --sort-tmpdir=tmpdir <"$vcf_filtered_mask" | cut -f 1,2,3) \
    >"$bed_removed_sites"

## Report:
echo "## linecount removed-sites: $(wc -l <"$bed_removed_sites")"
echo "## done with script."
date
