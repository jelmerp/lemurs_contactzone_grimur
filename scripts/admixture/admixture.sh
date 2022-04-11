#!/bin/bash

#SBATCH --mem=8G
#SBATCH --ntasks=1

set -euo pipefail

# SET-UP -----------------------------------------------------------------------
## Software
eval "$(conda shell.bash hook)"
conda activate admixture-env

## Commmand-line arguments
FILE_ID=$1
INDIR=$2
OUTDIR=$3
K=$4

## Process arguments
INPUT=$INDIR/$FILE_ID.bed
OUTLOG=$OUTDIR/$FILE_ID.$K.admixtureOutLog.txt

## Other variables
n_cores=1
#n_cores="$SLURM_CPUS_PER_TASK"

## Report
echo -e "\n## Starting script admixture.sh"
date
echo
echo "## File ID:                 $FILE_ID"
echo "## Value of K:              $K"
echo "## Number of cores:         $n_cores"
echo "## Plink input file dir:    $INDIR"
echo "## Output file dir:         $OUTDIR"
echo
echo "## Input bed file:          $INPUT"
echo "## Output file:             $OUTLOG"
echo "## Listing input file:"
ls -lh "$INPUT"
echo -e "--------------------------------\n"

## Make output dir if it doesn't exist
[[ ! -d "$OUTDIR"/pfiles ]] && mkdir -p "$OUTDIR"/pfiles

# CREATE INDIV FILE ------------------------------------------------------------
## Create a list of individuals used in the analysis directly from one of the
## files with the sequence data. The ADMIXTURE output doesn't contain this info.

echo -e "\n-----------------------"
INDIV_FILE=$OUTDIR/$FILE_ID.indivs.txt
echo "## Indiv file: $INDIV_FILE"

if [ ! -e "$OUTDIR"/"$FILE_ID".indivs.txt ]; then
    echo "## Creating indiv file from pedfile..."

    PEDFILE="$INDIR"/"$FILE_ID".ped
    echo "## Ped file: $PEDFILE"

    cut -f1 "$PEDFILE" >"$INDIV_FILE"
    echo "## Listing indiv file:"
    ls -lh "$INDIV_FILE"
fi

echo "Showing indiv file..."
cat "$INDIV_FILE"

# RUN ADMIXTURE ----------------------------------------------------------------
echo -e "\n-----------------------"
echo -e "## Running ADMIXTURE..."
admixture --cv -j"$n_cores" "$INPUT" "$K" >"$OUTLOG"


# HOUSEKEEPING -----------------------------------------------------------------
## Show outlog
echo -e "\n-----------------------"
echo -e "#### Showing contents of ADMIXTURE logfile..."
cat "$OUTLOG"

## Move files
echo -e "\n## Moving files to output dir..."
ls -lh "$FILE_ID"."$K".Q
ls -lh "$FILE_ID"."$K".P

mv "$FILE_ID"."$K".Q "$OUTDIR"/
mv "$FILE_ID"."$K".P "$OUTDIR"/pfiles/

## Report
echo -e "\n## Done with script admixture_run.sh"
date
echo
