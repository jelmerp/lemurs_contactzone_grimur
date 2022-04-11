#!/bin/bash

set -euo pipefail

echo -e "\n## Starting script demultiplex.sh"
date
echo

# SET-UP -----------------------------------------------------------------------
## Software
STACKS_PROCESS_RADTAGS=software/stacks-2.0b/process_radtags

## Command-line args:
R1=$1
R2=$2
OUTDIR=$3
DIR_STATS=$4
BARCODE_FILE=$5

## Prep:
[[ $R1 =~ \.gz$ ]] && INPUT_COMMAND="-i gzfastq" && ID=$(basename $R1 .fastq.gz)
[[ $R1 =~ \.fastq$ ]] && INPUT_COMMAND="-i fastq" && ID=$(basename $R1 .fastq)
STATSFILE=$DIR_STATS/$ID.demultiplexStats.txt

## Report:
echo "## Input fastq R1: $R1"
echo "## Input fastq R2: $R2"
echo "## Output fastq dir: $OUTDIR"
echo "## Output stats dir: $DIR_STATS"
echo "## Output stats file: $STATSFILE"
echo "## File with barcodes: $BARCODE_FILE"
echo
echo "## ID: $ID"
echo "## Stacks input format command: $INPUT_COMMAND"
echo

## Create output dir if not there yet:
mkdir -p "$OUTDIR"


################################################################################
#### RUN STACKS PROCESS_RADTAGS ####
################################################################################
$STACKS_PROCESS_RADTAGS --paired -r --inline_null $INPUT_COMMAND -e sbfI \
	-1 $R1 -2 $R2 -o $OUTDIR -b $BARCODE_FILE -y gzfastq > $STATSFILE
# -r option: rescue barcodes and radtags
# -e option: specify enzyme (=sbfI)
# --inline_null option: barcode is inline with sequence, occurs only on single-end read (default).

## Remove "rem" files:
echo -e "\n## Removing rem files..."
rm -f $OUTDIR/*rem.1.fq*
rm -f $OUTDIR/*rem.2.fq*

## Report:
echo -e "\n## Showing stats file contents:"
cat "$STATSFILE"

echo -e "\n## Done with script."
date
echo
