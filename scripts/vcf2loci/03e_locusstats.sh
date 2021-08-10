#!/bin/bash
set -euo pipefail

## Software:
AMAS=AMAS.py # https://github.com/marekborowiec/AMAS

## Command-line args:
dir_fasta=$1
file_stats_all=$2

## Process args:
dir_stats=$(dirname "$file_stats_all")

## Report:
echo "##  Starting script $0."
echo "##  Fasta dir:           $dir_fasta"
echo "##  Output file:         $file_stats_all"

## Initiate output file:
>"$file_stats_all"

## Get stats for each fasta file:
echo "##  Cycling through fasta files..."
for fasta in "$dir_fasta"/*; do
    echo "##  Fasta file: $fasta"

    fasta_id=$(basename "$fasta")
    file_stats=$dir_stats/tmp."$fasta_id".stats.txt

    $AMAS summary -f fasta -d dna -i "$fasta" -o "$file_stats"

    grep -v "Alignment_name" "file_stats" >>"$file_stats_all"
done

## Include header with column names:
header=$(head -n 1 "$file_stats")
(echo "$header" && cat "$file_stats_all") >"$dir_stats"/tmp.txt && mv "$dir_stats"/tmp.txt "$file_stats_all"

## Housekeeping:
echo "##  Removing temporary files:"
find "$dir_stats" -name 'tmp*txt' -delete

echo -e "##  Done with script $0."
date
