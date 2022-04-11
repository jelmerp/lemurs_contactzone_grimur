#!/bin/bash

set -euo pipefail

## SET-UP ------------------------------------------------------------------------
## Scripts and software
SCRIPT_VCF2TREEMIX=scripts/conversion/vcf2treemix.sh
SCRIPT_TREEMIX=scripts/treemix/treemix.sh

## Command-line args
FILE_ID=$1              # VCF file ID (File should be $VCF_DIR/$FILE_ID.vcf.gz
VCF_DIR=$2              # VCF dir 
PREP_INPUT=$3           # TRUE/FALSE, whether or not to create the input files from a VCF
MINMIG=$4               # Minimum nr of migration edges
MAXMIG=$5               # Maximum nr of migration edges
ROOT=$6                 # Root taxon
TREEMIX_DIR=$7          # Treemix base dir
INDS_METADATA=$8        # Metadata file
GROUP_BY_COLUMN=$9      # Column in metadata file to group individuals by (e.g. the column with the species name)

## Define Treemix input and output dirs:
TREEMIX_INDIR=$TREEMIX_DIR/input
TREEMIX_OUTDIR=$TREEMIX_DIR/output

## Make dirs if needed
mkdir -p "$TREEMIX_INDIR" "$TREEMIX_OUTDIR"

## Report
echo -e "\n## Starting with script 01_treemix_runner.sh"
date
echo
echo "## File ID:                $FILE_ID"
echo "## VCF dir:                $VCF_DIR"
echo "## Prep input:             $PREP_INPUT"
echo "## Min nr of mig events:   $MINMIG"
echo "## Max nr of mig events:   $MAXMIG"
echo "## Root:                   $ROOT"
echo "## Treemix input dir:      $TREEMIX_INDIR"
echo "## Treemix output dir:     $TREEMIX_OUTDIR"
echo
echo "## Metadata file:          $INDS_METADATA"
echo "## Column in metadata file to group individuals by: $GROUP_BY_COLUMN"
echo


# PREPARE TREEMIX INPUT --------------------------------------------------------
if [ "$PREP_INPUT" == "TRUE" ]; then
	echo -e "\n## Preparing Treemix input - calling vcf2treemix.sh...\n"
	
    $SCRIPT_VCF2TREEMIX "$FILE_ID" "$VCF_DIR" "$TREEMIX_DIR" "$INDS_METADATA" "$GROUP_BY_COLUMN"
else
	echo -e "\n## Skipping input preparation...\n"
fi


# RUN TREEMIX ------------------------------------------------------------------
echo -e "\n## Running Treemix...\n"

K=1000 # Set blocksize to 1000

for NMIG in $(seq "$MINMIG" "$MAXMIG"); do
	sbatch --mem 50G $SCRIPT_TREEMIX
        "$FILE_ID" "$NMIG" "$K" "$ROOT" "$TREEMIX_INDIR" "$TREEMIX_OUTDIR"
done


## Report:
echo -e "\n## Done with script 01_treemix_runner.sh"
date
echo
