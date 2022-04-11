#!/bin/bash

set -euo pipefail

## Software
EAUTILS_FASTQSTATS=software/ExpressionAnalysis-ea-utils-bd148d4/clipper/fastq-stats

## Command-line args
FASTQ=$1
OUTDIR=$2

## Process args
FASTQ_ID=$(basename -s .fastq.gz "$FASTQ")
OUTPUT="$OUTDIR"/"$FASTQ_ID".fastqstats.txt

## Report
echo "## Starting script."
date
echo
echo "## Input fastq file:    $FASTQ"
echo "## Outdir:              $OUTDIR"
echo "## Output file:         $OUTPUT"
echo

## Run
$EAUTILS_FASTQSTATS "$FASTQ" > "$OUTPUT"

## Report
echo "## Output file:"
ls -lh "$OUTPUT"
echo -e "Done with script fastq.stats.sh"
date
echo
