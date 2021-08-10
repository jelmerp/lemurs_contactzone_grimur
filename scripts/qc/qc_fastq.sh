#!/bin/bash

set -euo pipefail

module load fastqc

echo "## Starting script fastqc.sh"
date

## Command-line args:
if [ "$#" -ne 2 ]; then
    "Usage: fastqc.sh <input-file.fastq[.gz]> <output-dir>"
    "Exiting."
    exit 1
fi

input="$1"
outdir="$2"

## Report:
echo "## Input file:         $input"
echo "## Output directory:   $outdir"

[[ ! -d $outdir ]] && mkdir -p "$outdir"

echo "## Running fastqc..."
fastqc --outdir="$outdir" "$input"

echo -e "\n Done with script."
date
