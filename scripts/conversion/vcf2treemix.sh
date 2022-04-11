#!/bin/bash

set -euo pipefail

# SET-UP -----------------------------------------------------------------------
## Software and scripts
SCRIPT_POPFILE=scripts/treemix/treemix_makePopfile.R
SCRIPT_VCF2TREEMIX=scripts/conversion/vcf2treemix.py
# NOTE: SCRIPT ASSUMES THAT `bcftools` IS IN THE $PATH

## Command-line arguments:
FILE_ID=$1
VCF_DIR=$2
TREEMIX_DIR=$3
INDS_METADATA=$4
GROUP_BY_COLUMN=$5

## Process:
VCF=$VCF_DIR/$FILE_ID.vcf
TREEMIX_INFILE=$TREEMIX_DIR/input/$FILE_ID.tmix
INDFILE=$TREEMIX_DIR/popfiles/$FILE_ID.inds.txt
POPFILE=$TREEMIX_DIR/popfiles/$FILE_ID.popfile.txt

## Report:
echo -e "\n## Starting with script vcf2treemix.sh"
date
echo
echo "## File ID: $FILE_ID"
echo "## VCF dir: $VCF"
echo "## Treemix base dir: $TREEMIX_DIR"
echo "## Metadata file: $INDS_METADATA"
echo "## Column in metadata file to group individuals by: $GROUP_BY_COLUMN"
echo
echo "## VCF file: $VCF"
echo "## Indfile (to create): $INDFILE"
echo "## Popfile (to create): $POPFILE"
echo "## Treemix input file (to create): $TREEMIX_INFILE"
echo

## Check for files:
[[ -f "$VCF" ]] && echo -e "## Unzipped VCF found.\n"
[[ ! -f "$VCF" ]] && echo -e "## Unzipped VCF not found, unzipping...\n" && gunzip -c "$VCF".gz > "$VCF"

echo "## Listing VCF file:"
ls -lh "$VCF"
echo

## Create output dirs if needed
mkdir -p "$TREEMIX_DIR"/input "$TREEMIX_DIR"/popfiles


# CREATE TREEMIX POPFILE -------------------------------------------------------
echo "## Creating Treemix popfile..."
bcftools query -l "$VCF" > "$INDFILE"
$SCRIPT_POPFILE "$INDFILE" "$INDS_METADATA" "$POPFILE" "$GROUP_BY_COLUMN"

echo -e "\n## Showing contents of Treemix popfile:"
cat "$POPFILE"
echo


# CONVERT VCF TO TREEMIX FORMAT ------------------------------------------------
echo -e "\n## Converting vcf to Treemix format..."
python3 "$SCRIPT_VCF2TREEMIX" -vcf "$VCF" -pop "$POPFILE"
	
## Move and zip treemix input
echo -e "\n## Moving and zipping treemix input..."
mv "$VCF_DIR"/"$FILE_ID"*tmix "$TREEMIX_INFILE"
gzip -f "$TREEMIX_INFILE"


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing treemix input file:"
ls -lh "$TREEMIX_INFILE".gz
echo -e "\n## Done with script vcf2treemix.sh"
date
echo
