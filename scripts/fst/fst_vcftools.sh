#!/bin/bash

set -euo pipefail

## Software
eval "$(conda shell.bash hook)"
conda activate vcftools-env

## Command-line arguments:
INDIR=$1
FILE_ID=$2
POPFILE1=$3
POPFILE2=$4
WINSIZE=$5
WINSTEP=$6
OUTDIR=$7
POPFILEDIR=$8

## Process arguments:
VCF=$INDIR/$FILE_ID.vcf.gz

POP1=$(basename "$POPFILE1" .txt)
POP2=$(basename "$POPFILE2" .txt)

OUTPUT=$OUTDIR/"$POP1"_vs_"$POP2"_win$WINSIZE.step$WINSTEP.txt

[[ ! -d "$OUTDIR" ]] && mkdir -p "$OUTDIR"

## Report:
echo "## Starting script fst_vcftools.sh"
date
echo "## VCF source: $VCF"
echo "## Population file 1: $POPFILE1"
echo "## Population file 2: $POPFILE2"
echo "## Population 1: $POP1"
echo "## Population 2: $POP2"
echo "## Window size: $WINSIZE"
echo "## Window stepsize: $WINSTEP"
echo "## Output dir: $OUTDIR"
echo "## Popfile dir: $POPFILEDIR"
echo "## Output: $OUTPUT"

## Run vcftools:
vcftools --gzvcf "$VCF" \
    --fst-window-size "$WINSIZE" --fst-window-step "$WINSTEP" \
    --weir-fst-pop "$POPFILE1" --weir-fst-pop "$POPFILE2" --out "$OUTPUT"

## Report:
echo "Done with script."
date
