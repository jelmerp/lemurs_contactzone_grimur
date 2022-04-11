#!/bin/bash

set -euo pipefail

# SET-UP -----------------------------------------------------------------------
## Software
BCFTOOLS=software/bcftools-1.6/bcftools
TABIX=software/htslib-1.6/tabix
BGZIP=software/htslib-1.6/bgzip
VCFTOOLS=software/vcfrools/vcftools-master/bin/vcftools

## Command-line args
FILE_ID=$1                # File ID (basename) for VCF
VCF_DIR=$2                # Dir with VCF
QC_DIR=$3                 # Dir for QC outpur
RUN_BCF=$4                # FALSE/TRUE
RUN_BCF_BY_IND=$5         # FALSE/TRUE
RUN_VCF_SITESTATS=$6      # FALSE/TRUE

## Process:
VCF="$VCF_DIR"/"$FILE_ID".vcf.gz

## Report:
echo "## Starting script qc_vcf.sh"
date
echo
echo "## Indiv ID: $FILE_ID"
echo "## VCF dir: $VCF_DIR"
echo "## Input VCF file: $VCF"
echo "## Base QC dir: $QC_DIR"
echo "## Run bcftoolsStats: $RUN_BCF"
echo "## Run bcftoolsStats by ind: $RUN_BCF_BY_IND"
echo
echo "## Run vcftools sitestats modules: $RUN_VCF_SITESTATS"
echo

## Index VCF:
if [ ! -e "$VCF".tbi ]; then
	echo "## No .tbi file found..."
	
	if [ -e "$VCF_DIR"/"$FILE_ID".vcf ]; then
	
    	echo "## Rezipping vcf using bgzip..."
		$BGZIP -f "$VCF_DIR"/"$FILE_ID".vcf
	
    elif [ -e "$VCF" ]; then
	
    	echo "## Rezipping vcf using bgzip..."
		gunzip -f "$VCF"
		$BGZIP -f "$VCF_DIR"/"$FILE_ID".vcf	
	
    else
		echo -e "\n## ERROR: VCF FILE NOT FOUND!!\n\n" && exit 1
	fi
	
	echo "## Indexing vcf with tabix..." 
	$TABIX -f -p vcf "$VCF"
else
	echo "## .tbi file found."
fi

## List input VCF
echo -e "## Listing input VCF file:"
ls -lh "$VCF"

## Create directories
mkdir -p "$QC_DIR" "$QC_DIR"/bcftools "$QC_DIR"/vcftools

## Defaults if "RUN_BCF" and "RUN_BCF_BY_IND" are not assigned
[[ -z $RUN_BCF ]] && RUN_BCF=TRUE
[[ -z $RUN_BCF_BY_IND ]] && RUN_BCF_BY_IND=FALSE
[[ -z $RUN_VCF_SITESTATS ]] && RUN_VCF_SITESTATS=TRUE


# RUN BCFTOOLS-STATS ON ENTIRE VCF ####
if [ $RUN_BCF = TRUE ]; then
	echo -e "\n------------------"
	echo "## Running bcftools-stats on entire vcf..."
	
	OUTFILE=$QC_DIR/bcftools/$FILE_ID.bcftools.txt
	
	$BCFTOOLS stats --samples - "$VCF" > "$OUTFILE"
	
	echo -e "\n## bcftools-stats output file: $OUTFILE"
	ls -lh "$OUTFILE"
fi


# RUN BCFTOOLS-STATS SAMPLE-BY-SAMPLE ------------------------------------------
IDs_SINGLE=( $($BCFTOOLS query -l $VCF) )

if [ $RUN_BCF_BY_IND == TRUE ] && [ ${#IDs_SINGLE[@]} -gt 1 ]; then
	echo -e "\n------------------"
	echo "## Running bcftools-stats sample-by-sample..."
	echo "## Number of samples in vcf: ${#IDs_SINGLE[@]}"
	echo
	
	for ID_SINGLE in ${IDs_SINGLE[@]}; do
		echo "## ID: $ID_SINGLE"
		OUTFILE_IND=$QC_DIR/bcftools/$FILE_ID.$ID_SINGLE.bcftools.txt
		$BCFTOOLS stats --samples "$ID_SINGLE" "$VCF" > "$OUTFILE_IND"
	done
fi


# VCFTOOLS: DEPTH AND MISSINGNESS STATS ----------------------------------------
echo -e "\n## Running vcftools to get missing-indv stats..."
$VCFTOOLS --gzvcf "$VCF" --missing-indv --stdout > "$QC_DIR"/vcftools/"$FILE_ID".imiss

echo -e "\n## Running vcftools to get depth stats..."
$VCFTOOLS --gzvcf "$VCF" --depth --stdout > "$QC_DIR"/vcftools/'$FILE_ID'.idepth

if [ $RUN_VCF_SITESTATS == TRUE ]; then
	echo -e "\n## Running vcftools to get missing-site stats..."
	$VCFTOOLS --gzvcf "$VCF" --missing-site --stdout > "$QC_DIR"/vcftools/"$FILE_ID".smiss
	
	echo -e "\n## Running vcftools to get site-mean-depth stats..."
	$VCFTOOLS --gzvcf "$VCF" --site-mean-depth --stdout > "$QC_DIR"/vcftools/"$FILE_ID".sdepth
fi


# REPORT -----------------------------------------------------------------------
echo -e "\n## QC Output files:"
ls -lh "$QC_DIR"/vcftools/"$FILE_ID"*
[[ $RUN_BCF == TRUE ]] && ls -lh "$QC_DIR"/bcftools/"$FILE_ID"*

echo -e "\n#### idepth file:"
cat "$QC_DIR"/vcftools/"$FILE_ID".idepth

echo -e "\n## Done with script qc_vcf.sh"
date
echo
