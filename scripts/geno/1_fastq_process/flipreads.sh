#!/bin/bash
set -e
set -o pipefail
set -u

################################################################################
#### SET-UP ####
################################################################################
## Command-line args:
LIBRARY_ID=$1
INDIR=$2
OUTDIR=$3
STATSDIR=$4
BARCODE_FILE=$5
FLIP_SCRIPT=$6

## Prep:
NLINES_FILE=$STATSDIR/$LIBRARY_ID.flipstats.linenrs.txt # Output stats file
[[ ! -d $OUTDIR ]] && mkdir -p $OUTDIR
[[ ! -d $STATSDIR ]] && mkdir -p $STATSDIR

## Report:
date
echo "## flipreads.sh: Library ID: $LIBRARY_ID"
echo "## flipreads.sh: Indir: $INDIR"
echo "## flipreads.sh: Outdir: $OUTDIR"
echo "## flipreads.sh: Statsdir: $STATSDIR"
echo "## flipreads.sh: File with barcodes: $BARCODE_FILE"
echo "## flipreads.sh: Flipping script: $FLIP_SCRIPT"
echo "## flipreads.sh: File with line numbers: $NLINES_FILE"
echo


################################################################################
#### FLIP READS #####
################################################################################
## Get files:
FILES=$(ls $INDIR/*${LIBRARY_ID}*)
echo -e "## flipreads.sh: Fastq files:\n$FILES"

## Define read-files and prefix:
READ1_IN=$(echo $FILES | cut -d" " -f1)
READ2_IN=$(echo $FILES | cut -d" " -f2)
echo -e "\n## flipreads.sh: File with R1 reads:" && ls -lh $READ1_IN
echo -e "\n## flipreads.sh: File with R2 reads:" && ls -lh $READ2_IN

[[ $READ1_IN =~ \.gz$ ]] && \
  PREFIX=$(basename -s "_R1_001.fastq.gz" $READ1_IN) || \
  PREFIX=$(basename -s "_R1_001.fastq"  $READ1_IN)
echo -e "\n## flipreads.sh: Prefix: $PREFIX"

READ1_OUT=$OUTDIR/${PREFIX}_R1_flipped.fastq
READ2_OUT=$OUTDIR/${PREFIX}_R2_flipped.fastq
echo "## flipreads.sh: Read 1 outfile: $READ1_OUT"
echo "## flipreads.sh: Read 2 outfile: $READ2_OUT"

## Unzip:
[[ $READ1_IN =~ \.gz$ ]] && echo "## flipreads.sh: Unzipping READ1_IN..." && gunzip $READ1_IN
[[ $READ2_IN =~ \.gz$ ]] && echo "## flipreads.sh: Unzipping READ2_IN..." && gunzip $READ2_IN

## Names to use for flipping script:
READ1_IN=$(echo $READ1_IN | sed 's/.gz//')
READ2_IN=$(echo $READ2_IN | sed 's/.gz//')


################################################################################
#### RUN FLIPPING SCRIPT ####
################################################################################
echo -e "\n## flipreads.sh: Running Perl flipping script..."
$FLIP_SCRIPT $BARCODE_FILE $READ1_IN $READ2_IN $READ1_OUT $READ2_OUT > $STATSDIR/flipStats.$PREFIX


################################################################################
#### HOUSEKEEPING ####
################################################################################
## Report:
echo -e "\n## flipreads.sh: Reporting line numbers:"
wc -l $READ1_IN
wc -l $READ2_IN
wc -l $READ1_OUT
wc -l $READ2_OUT

wc -l $READ1_IN >> $NLINES_FILE
wc -l $READ2_IN >> $NLINES_FILE
wc -l $READ2_OUT >> $NLINES_FILE
wc -l $READ2_OUT >> $NLINES_FILE

## Gzip:
echo -e "\n## flipreads.sh: Gzipping fastqs..."
gzip $READ1_OUT
gzip $READ2_OUT

## Report:
echo -e "\n## flipreads.sh: Done with script."
date
