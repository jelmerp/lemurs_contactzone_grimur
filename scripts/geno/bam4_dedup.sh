#!/bin/bash
set -e
set -o pipefail
set -u

date
echo "## bam4_dedup.sh: Starting script."
echo

################################################################################
#### SET-UP ####
################################################################################
## Command-line args:
ID=$1
INDIR=$2
OUTDIR=$3
PREFIX_IN=$4
PREFIX_OUT=$5
BAMSTATS_DIR=$6

## Process args:
BAM_IN=$INDIR/$ID.$PREFIX_IN.bam
BAM_OUT=$OUTDIR/$ID.$PREFIX_OUT.bam
STATSFILE=$BAMSTATS_DIR/$ID.bamfilterstats.txt

## Report:
echo "## bam4_dedup.sh: ID: $ID"
echo "## bam4_dedup.sh: Indir: $INDIR"
echo "## bam4_dedup.sh: Outdir: $OUTDIR"
echo "## bam4_dedup.sh: Prefix_in: $PREFIX_IN"
echo "## bam4_dedup.sh: Prefix_out: $PREFIX_OUT"
echo "## bam4_dedup.sh: Bamstats dir: $BAMSTATS_DIR"
echo "## bam4_dedup.sh: Bamstats file: $STATSFILE"

[[ ! -d tmpdir ]] && mkdir tmpdir

################################################################################
#### DEDUP ####
################################################################################
echo "## bam4_dedup.sh: Deduplicating bamfile..."
# First sort by name, then run fixmate (-r: remove secondary&unmapped reads; -m: add mate score tags), then markdup
samtools sort -T tmpdir -n $BAM_IN | \
  samtools fixmate -r -m - - | \
  samtools sort -T tmpdir - | \
  samtools markdup -r - $BAM_OUT                    # -r: remove duplicates

echo "## bam4_dedup.sh: Indexing bamfile..."
samtools index -b $BAM_OUT                          # index bamfile

################################################################################
#### HOUSEKEEPING ####
################################################################################
NRSEQS_IN=$(samtools view -c $BAM_IN)              # Nr of sequences in input file
NRSEQS_OUT=$(samtools view -c $BAM_OUT)            # Nr of sequences in output file
echo -e "\n## bam4_dedup.sh: Nr of sequences before dedupping: $NRSEQS_IN"
echo -e "## bam4_dedup.sh: Nr of sequences after dedupping: $NRSEQS_OUT \n"

## Write to statsfile:
echo "Nr of sequences before dedupping: $NRSEQS_IN" >> $STATSFILE 
echo "Nr of sequences after dedupping: $NRSEQS_OUT" >> $STATSFILE

echo -e "## bam4_dedup.sh: Done with script. \n"
date