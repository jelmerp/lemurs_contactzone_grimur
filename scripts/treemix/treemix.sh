#!/bin/bash

set -euo pipefail

# SET-UP -----------------------------------------------------------------------
## Software
TREEMIX=software/treemix-1.13/src/treemix

## Command-line args
FILE_ID=$1         # File ID for Treemix input file
NMIG=$2            # Nr of migration edges to be added
K=$3               # Blocksize
ROOT=$4            # Root taxon
TRMX_INDIR=$5      # Treemix input dir
TRMX_OUTDIR=$6     # Treemix output dir

## Define Treemix input and output files:
TRMX_INPUT=$TRMX_INDIR/$FILE_ID.tmix
TRMX_OUTPUT=$TRMX_OUTDIR/$FILE_ID.treemixOutput.k$K.mig$NMIG

## Report:
date
echo "## Starting with script."
echo "## file ID:                      $FILE_ID"
echo "## Number of migration events:   $NMIG"
echo "## Value of K:                   $K"
echo "## Root:                         $ROOT"
echo "## Treemix input:                $TRMX_INPUT"
echo "## Treemix output:               $TRMX_OUTPUT"
echo

# RUN TREEMIX ------------------------------------------------------------------
if [ "$ROOT" != 'none' ]; then
	
    echo "## Running treemix with root..."
	$TREEMIX -i "$TRMX_INPUT".gz -root "$ROOT" -m "$NMIG" -k "$K" \
        -o "$TRMX_OUTPUT".root"$ROOT"

else
	
    echo "## Running treemix without root..."
	$TREEMIX -i "$TRMX_INPUT".gz -m "$NMIG" -k "$K" -o "$TRMX_OUTPUT".rootNone

fi


echo -e "\n## Done with script treemix.sh"
date
echo
