#!/bin/bash
set -euo pipefail

## Setup:
file_inds=$1
dir_indfasta=$2
locusbed_final=$3
locuslist=$4
fasta_merged=$5

## For each individual, extract loci in locus-bedfile from altrefmasked fasta
while read -r ind; do
    fasta_in="$dir_indfasta"/"$ind"_altrefmasked.fasta
    fasta_out="$dir_indfasta"/"$ind"_allloci.fasta
    bedtools getfasta -fi "$fasta_in" -bed "$locusbed_final" >"$fasta_out"
done <"$file_inds"

## make list with loci
fasta_1=$(find "$dir_indfasta"/*allloci.fasta | head -1)
grep ">" "$fasta_1" | sed 's/>//' >"$locuslist"

## merge by-individual fasta files
>"$fasta_merged"

while read -r ind; do
    fasta=$dir_indfasta/"$ind"_allloci.fasta
    # replacing ":" by "," for compatibility with faidx:
    sed "s/>/>${ind}__/g" "$fasta" | sed 's/:/,/g' >>"$fasta_merged"
done <"$file_inds"

## index merged fasta file:
samtools faidx "$fasta_merged"

## report:
echo "Done with script $0"
