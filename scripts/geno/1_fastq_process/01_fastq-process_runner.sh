#!/bin/bash

set -euo pipefail

echo -e "## Starting with script fastq-process_runner.sh"
date
echo -e "## Args:\n$@\n"

Help() {
    # Display Help
    echo "## Runner script to process RADseq FASTQ files."
    echo
    echo "## Syntax: fqprocess_pip.sh -b barcode-file -i input-dir -l library-ID -o output-dir -s steps-to-skip [-h]"
    echo "## Options:"
    echo "## -h     Print help."
    echo "## -s     Steps-to-skip, provide one or more letters:"
    echo "##        A: Skip all steps, dry run."
    echo "##        F: Skip read flipping step."
    echo "##        D: Skip demultiplexing step."
    echo "##        T: Skip read trimming."
    echo "##        M: Skip deduplication (removing PCR duplicates)."
    echo "##        Q: Skip QC."
    echo "##        S: Skip stats."
    echo
}

# SET-UP -----------------------------------------------------------------------

## Scripts:
SCRIPT_FLIP_PERL=scripts/geno/1_fastq_proces/flip_trim_sbfI_170601.pl
SCRIPT_FLIP_BASH=scripts/geno/1_fastq_proces/flipreads.sh
SCRIPT_DEMULT=scripts/geno/1_fastq_proces/demultiplex.sh
SCRIPT_DEDUP_TRIM=scripts/geno/1_fastq_proces/fqdeduptrim.sh
SCRIPT_CHECKBARCODES=scripts/geno/1_fastq_proces/checkbarcodes.sh
SCRIPT_QC_FASTQ=scripts/qc/qc_fastq.sh
SCRIPT_STATS_FASTQ=scripts/qc/qc_fastq_stats.sh

## Hardcoded variables:
ADAPTER_FILE=/datacommons/yoderlab/programs/Trimmomatic-0.36/adapters/all.fa # For trimmomatic
CUTSITE="TGCAGG"                                                             # sbfI enzyme
NCORES_DEDUP=4

## Command-line args:
BARCODE_TABLE=""
INDIR=""
LIBRARY_ID=""
OUTDIR=""
SKIP=""

while getopts ':b:i:l:o:s:h' flag; do
    case "${flag}" in
    b) BARCODE_TABLE="$OPTARG" ;;
    i) INDIR="$OPTARG" ;;
    l) LIBRARY_ID="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    s) SKIP="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## ERROR: Invalid option" && exit 1 ;;
    :) echo "## ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done
[[ -z $BARCODE_TABLE ]] && echo "## ERROR: Please provide barcode file with -b flag" && exit 1
[[ -z $INDIR ]] && echo "## ERROR: Please provide input dir with -i flag" && exit 1
[[ -z $LIBRARY_ID ]] && echo "## ERROR: Please provide library ID with -l flag" && exit 1
[[ -z $OUTDIR ]] && echo "## ERROR: Please provide output dir with -o flag" && exit 1

echo "## Library ID: $LIBRARY_ID"
echo "## Input dir: $INDIR"
echo "## Output dir: $OUTDIR"
echo "## Barcode file: $BARCODE_TABLE"
echo "## To-skip-string: $SKIP"
echo

# To skip:
[[ $SKIP =~ "Q" ]] && SKIP_QC="TRUE" || SKIP_QC="FALSE"
[[ $SKIP =~ "F" ]] && SKIP_FLIP="TRUE" || SKIP_FLIP="FALSE"
[[ $SKIP =~ "M" ]] && SKIP_DEMULT="TRUE" || SKIP_DEMULT="FALSE"
[[ $SKIP =~ "D" ]] && SKIP_DEDUP="TRUE" || SKIP_DEDUP="FALSE"
[[ $SKIP =~ "T" ]] && SKIP_TRIM="TRUE" || SKIP_TRIM="FALSE"
[[ $SKIP =~ "S" ]] && SKIP_STATS="TRUE" || SKIP_STATS="FALSE"
[[ $SKIP =~ "A" ]] && SKIP_QC="TRUE" && SKIP_FLIP="TRUE" && SKIP_DEMULT="TRUE" &&
    SKIP_DEDUP="TRUE" && SKIP_TRIM="TRUE" && SKIP_STATS="TRUE"

echo "## Skip QC step: $SKIP_QC"
echo "## Skip flipping step: $SKIP_FLIP"
echo "## Skip demultiplexing step: $SKIP_DEMULT"
echo "## Skip dedupping: $SKIP_DEDUP"
echo "## Skip trimming: $SKIP_TRIM"
echo

## Process args:
DIR_RAW=$INDIR
DIR_FLIPPED=$OUTDIR/flipped
DIR_DEMULT=$OUTDIR/demult
DIR_FINAL=$OUTDIR/final
DIR_STATS=$OUTDIR/QC
DIR_STATS_IND=$DIR_STATS/by_ind

[[ ! -d $DIR_FLIPPED ]] && mkdir -p $DIR_FLIPPED
[[ ! -d $DIR_DEMULT ]] && mkdir -p $DIR_DEMULT
[[ ! -d $DIR_FINAL ]] && mkdir -p $DIR_FINAL
[[ ! -d $DIR_STATS_IND ]] && mkdir -p $DIR_STATS_IND

BARCODE_LIST=$DIR_STATS/barcodelist.$LIBRARY_ID
cut -f 1 "$BARCODE_TABLE" >"$BARCODE_LIST"
echo "## Barcode list:" && ls -lh "$BARCODE_LIST" && echo

## Report:
echo "## Dir with raw fastqs: $DIR_RAW"
echo "## Dir with flipped fastqs: $DIR_FLIPPED"
echo "## Dir with demultiplexed fastqs: $DIR_DEMULT"
echo "## Dir with final/processed fastqs: $DIR_FINAL"
echo "## Dir with qc stats: $DIR_STATS"
echo "## Dir with qc stats by ind: $DIR_STATS_IND"

## Perform checks:
[[ ! -s $SCRIPT_FLIP_PERL ]] && echo "ERROR: $SCRIPT_FLIP_PERL not found" && exit 1
[[ ! -s $SCRIPT_FLIP_BASH ]] && echo "ERROR: $SCRIPT_FLIP_BASH not found" && exit 1
[[ ! -s $SCRIPT_DEMULT ]] && echo "ERROR: $SCRIPT_DEMULT not found" && exit 1
[[ ! -s $SCRIPT_DEDUP_TRIM ]] && echo "ERROR: $SCRIPT_DEDUP_TRIM not found" && exit 1
[[ ! -s $SCRIPT_QC_FASTQ ]] && echo "ERROR: $SCRIPT_QC_FASTQ not found" && exit 1
[[ ! -s $SCRIPT_STATS_FASTQ ]] && echo "ERROR: $SCRIPT_STATS_FASTQ not found" && exit 1
[[ ! -s $SCRIPT_CHECKBARCODES ]] && echo "ERROR: $SCRIPT_CHECKBARCODES not found" && exit 1
[[ ! -s $ADAPTER_FILE ]] && echo "ERROR: $ADAPTER_FILE not found" && exit 1

# 0. QC RAW FILES --------------------------------------------------------------

echo -e "\n------------------------"

if [ $SKIP_QC = "FALSE" ]; then

    for FASTQ in "$DIR_RAW"/*"$LIBRARY_ID"*; do

        ID=$(basename "$FASTQ" .fastq.gz)
        echo -e "\n## Running qc on raw file $ID..."
        ls -lh "$FASTQ"

        echo -e "## Submitting fastqc script:..."
        sbatch -o slurm.fastqc."$ID" "$SCRIPT_QC_FASTQ" "$FASTQ" "$DIR_STATS"

        echo "## Submitting check-barcodes script:..."
        OUTFILE=$DIR_STATS/barcodecounts_$ID.txt
        sbatch -o slurm.checkbarcodes."$ID" "$SCRIPT_CHECKBARCODES" "$FASTQ" "$CUTSITE" "$OUTFILE"
    done

else
    echo "## Skipping QC step...."
fi

# 1. FLIP READS ----------------------------------------------------------------

echo -e "\n------------------------"

if [ "$SKIP_FLIP" = FALSE ]; then

    echo -e "## Step 1 -- flip reads...\n"

    "$SCRIPT_FLIP_BASH" "$LIBRARY_ID" "$DIR_RAW" "$DIR_FLIPPED" "$DIR_STATS" \
        "$BARCODE_LIST" "$SCRIPT_FLIP_PERL"
else
    echo "## Skipping flipping step...."
fi

## 2. DEMULTIPLEX ---------------------------------------------------------------

echo -e "\n------------------------"

if [ "$SKIP_DEMULT" = FALSE ]; then

    echo -e "## Step 2 -- demultiplex fastqs...\n"

    R1=$(ls "$DIR_FLIPPED"/*"$LIBRARY_ID"*_R1_flipped*fastq*)
    R2=$(ls "$DIR_FLIPPED"/*"$LIBRARY_ID"*_R2_flipped*fastq*)
    echo "## R1: $R1"
    echo "## R2: $R2"

    "$SCRIPT_DEMULT" "$R1" "$R2" "$DIR_DEMULT" "$DIR_STATS" "$BARCODE_TABLE"
else
    echo "## Skipping demupltiplexing step...."
fi

# 3. REMOVE PCR DUPS & TRIM ----------------------------------------------------

echo -e "\n------------------------"

if [ "$SKIP_DEDUP" = FALSE ] || [ "$SKIP_TRIM" = FALSE ]; then

    echo -e "## Step 3 -- dedup & trim fastqs...\n"

    ## Paired-end reads:
    for R1 in "$DIR_DEMULT"/*1.f*q.gz; do

        if [ -e "$R1" ]; then
            PREFIX=$(basename "$R1" ".1.f*q.gz" | sed 's/.R1.fastq.gz//' | sed 's/.R1.fq.gz//' | sed 's/.1.fq.gz//' | sed 's/.1.fastq.gz//')
            R2="$DIR_DEMULT"/$PREFIX*2.f*q.gz

            echo -e "\n## Paired-end sequences detected."
            echo "## Prefix: $PREFIX"
            echo "## Listing FASTQ files:"
            ls -lh "$R1"
            ls -lh "$R2"

            echo "## Running dedup/trim script..."
            sbatch --job-name=fqprocess."$PREFIX" --mem=16G -o slurm.deduptrim."$PREFIX" \
                "$SCRIPT_DEDUP_TRIM" "$R1" "$R2" "$PREFIX" "$DIR_FINAL" "$DIR_STATS_IND" "$ADAPTER_FILE" \
                "$NCORES_DEDUP" "$SKIP_DEDUP" "$SKIP_TRIM"

            ## Fastqc on final files
            if [ "$SKIP_STATS" = FALSE ]; then
                for FASTQ in "$DIR_FINAL"/*"$PREFIX"*fastq.gz; do
                    echo "## Running fastqc on final file $FASTQ..."
                    sbatch --dependency=singleton --job-name=fqprocess."$PREFIX" --mem=16G -o slurm.fastqc."$PREFIX" \
                        "$SCRIPT_QC_FASTQ" "$FASTQ" "$DIR_STATS_IND"
                done
            fi

        fi

    done

    ## Single-end reads:
    for R1 in "$DIR_DEMULT"/*0.f*q.gz; do

        if [ -e "$R1" ]; then
            PREFIX=$(basename $R1 ".0.fq.gz" | sed 's/.R0.fastq.gz//')
            R2="NONE"

            echo -e "\n## Single-end sequences detected."
            echo -e "## Prefix: $PREFIX"
            echo -e "## Listing FASTQ file:"
            ls -lh "$R1"

            sbatch --job-name=fqprocess."$PREFIX" --mem=16G -o slurm.deduptrim."$PREFIX" \
                "$SCRIPT_DEDUP_TRIM" "$R1" "$R2" "$PREFIX" "$DIR_FINAL" "$DIR_STATS_IND" $ADAPTER_FILE \
                "$NCORES_DEDUP" "$SKIP_DEDUP" "$SKIP_TRIM"

            ## Fastqc on final files:
            if [ "$SKIP_STATS" = FALSE ]; then
                for FASTQ in "$DIR_FINAL"/*"$PREFIX"*fastq.gz; do
                    echo -e "## Running fastqc on final file $FASTQ..."
                    sbatch --dependency=singleton --job-name=fqprocess."$PREFIX" --mem=16G -o slurm.fastqc."$PREFIX" \
                        $SCRIPT_QC_FASTQ "$FASTQ" "$DIR_STATS_IND"
                done
            fi

        fi

    done
else
    echo "## Skipping dedup & trim step...."
fi

echo -e "\n## Done with script."
date
