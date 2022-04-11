#!/bin/bash

set -euo pipefail

## Report
echo -e "\n## Starting with script."
date
echo

# SET-UP -----------------------------------------------------------------------
## Scripts & software
TRIMMOMATIC=software/Trimmomatic-0.36/trimmomatic-0.36.jar
STACKS_DEDUP=software/stacks-2.0b/clone_filter

## Command-line args:
R1=$1
R2=$2
PREFIX=$3
OUTDIR=$4
DIR_STATS=$5
ADAPTER_FILE=$6
NCORES=$7
SKIP_DEDUP=$8
SKIP_TRIM=$9

## Prep:
[[ ! -d $OUTDIR/discarded ]] && mkdir -p $OUTDIR/discarded
DEDUPSTATS_FILE=$DIR_STATS/$PREFIX.dedupstats.txt
TRIMSTATS_FILE=$DIR_STATS/$PREFIX.trimstats.txt

## Report:
echo "## Input R1: $R1"
echo "## Input R2: $R2"
echo "## Prefix (ID): $PREFIX"
echo "## Outdir: $OUTDIR"
echo "## Adapter file: $ADAPTER_FILE"
echo "## Number of cores: $NCORES"
echo
echo "## Dedup-stats file: $DEDUPSTATS_FILE"
echo "## Trim-stats file: $TRIMSTATS_FILE"
echo
echo "## Skip dedup (TRUE/FALSE): $SKIP_DEDUP"
echo "## Skip trimming (TRUE/FALSE): $SKIP_TRIM"


################################################################################
#### 1. DEDUP ####
################################################################################
echo -e "\n###################################################################"
if [ "$R2" != "NONE" ] && [ "$SKIP_DEDUP" = "FALSE" ]; then
	
    echo "## Dedup step..."
	echo "## Running stacks..."
	$STACKS_DEDUP -1 "$R1" -2 "$R2" -o "$OUTDIR" -i gzfastq 2>&1 | tee "$DEDUPSTATS_FILE"
	
	R1_IN=$OUTDIR/$PREFIX*1.1.fq.gz
	R2_IN=$OUTDIR/$PREFIX*2.2.fq.gz

elif [ $R2 = "NONE" ]; then

	echo "## Single-end sequences - skipping dedup step."
	R1_IN=$R1
	R2_IN=""

elif [ "$SKIP_DEDUP" = "TRUE" ]; then

	echo "## Skipping dedup step."
	R1_IN=$R1
	R2_IN=$R2

fi


################################################################################
#### 2. TRIM #####
################################################################################
echo -e "\n#####################################################################"
echo "## Trimming step..."

if [ "$R2" != "NONE" ] && [ "$SKIP_TRIM" = "FALSE" ]; then
	echo "## Trimming step - paired-end sequences."
	
	R1_OUT=$OUTDIR/$PREFIX.R1.fastq.gz
	R2_OUT=$OUTDIR/$PREFIX.R2.fastq.gz
	R1_DISCARD=$OUTDIR/discarded/$PREFIX.U1.fastq.gz
	R2_DISCARD=$OUTDIR/discarded/$PREFIX.U2.fastq.gz
	
	echo "## Trimming step - input R1: $R1_IN"
	echo "## Trimming step - input R2: $R2_IN"
	echo "## Trimming step - output R1: $R1_OUT"
	echo "## Trimming step - output R2: $R2_OUT"
	echo "## Trimming step - discarded R1 seqs: $R1_DISCARD"
	echo "## Trimming step - discarded R2 seqs: $R2_DISCARD"
	
	echo -e "\n## Running trimmomatic..."
	java -jar $TRIMMOMATIC PE -threads $NCORES -phred33 $R1_IN $R2_IN \
		$R1_OUT $R1_DISCARD $R2_OUT $R2_DISCARD \
		ILLUMINACLIP:$ADAPTER_FILE:2:30:10 \
		AVGQUAL:20 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:60 2>&1 | tee $TRIMSTATS_FILE

elif [ $R2 = "NONE" ] && [ $SKIP_TRIM = "FALSE" ]; then

	echo "## Trimming step - single-end sequences."
	
	R1_OUT=$OUTDIR/$PREFIX.R0.fastq.gz
	
	echo "## Trimming step - input R0: $R1_IN"
	echo "## Trimming step - output R0: $R1_OUT"
	
	echo -e "\n#### Running trimmomatic..."
	java -jar $TRIMMOMATIC SE -threads $NCORES -phred33 $R1_IN $R1_OUT \
		ILLUMINACLIP:$ADAPTER_FILE:2:30:10 \
		AVGQUAL:20 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:60 2>&1 | tee $TRIMSTATS_FILE

else
	echo "## Skipping trimming step."
fi

echo -e "\n## Removing intermediate files..."
[[ -s $R1_OUT ]] && rm -fv "$R1_IN"
[[ -s $R2_OUT ]] && rm -fv "$R2_IN"

echo -e "\n## Done with script fqdeduptrim.sh"
date
echo
