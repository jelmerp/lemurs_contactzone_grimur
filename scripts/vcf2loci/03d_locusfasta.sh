#!/bin/bash
set -euo pipefail

## Software:
# faidx

## Command-line args:
locuslist=$1
dir_locusfasta_intermed=$2
fasta_merged=$3

## Report:
date
echo "## Starting script $0."
echo "## Locus list:            $locuslist"
echo "## Fasta dir - by locus:  $dir_locusfasta_intermed"
echo "## Fasta - merged:        $fasta_merged"

## Extract single-locus FASTA files from merged FASTA
nlines=$(wc -l <"$locuslist")

echo -e "## Cycling through locus-list lines... \n"
for line in $(seq 1 "$nlines"); do

    locus_realname=$(head -n "$line" "$locuslist" | tail -n 1)
    locus_faidxname=$(echo "$locus_realname" | sed 's/:/,/g')
    fasta=$dir_locusfasta_intermed/$locus_realname.fa

    $FAIDX --regex "$locus_faidxname" "$fasta_merged" | sed 's/,/:/g' >"$fasta"

    echo "## Line: $line    Fasta: $fasta"
done

## Report:
echo -e "\n## Done with script."
date
