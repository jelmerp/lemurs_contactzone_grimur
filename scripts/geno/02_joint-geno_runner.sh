#!/bin/bash

set -euo pipefail

# SET-UP -----------------------------------------------------------------------
## Scripts:
SCR_GENO=scripts/geno/2_gatk/gatk2_jointgeno.sh
SCR_MERGEVCFS=scripts/geno/2_gatk/gatk3_mergevcfs.sh
SCR_FILTER=scripts/geno/3_filter/01_filterVCF_FS6_pip.sh

## Command-line args:
FILE_ID=$1
shift
scaf_file=$1
shift
GVCF_DIR=$1
shift
VCF_DIR=$1
shift
QC_DIR=$1
shift
REF=$1
shift
ADD_COMMANDS=$1
shift
MEM_JOB=$1
shift
MEM_GATK=$1
shift
NCORES=$1
shift
SKIP_GENO=$1
shift
DP_MEAN=$1
shift

count=0
while [ "$*" != "" ]
    do INDS[$count]=$1
    shift
    count=$(expr $count + 1)
done

## Process args:
VCF_DIR_MAIN=$VCF_DIR/intermed/
VCF_DIR_FINAL=$VCF_DIR/final/
VCF_DIR_SCAFFOLD=$VCF_DIR/intermed_byscaffold
SCAFFOLDLIST_DIR=$VCF_DIR/scaffoldlists

## Report:
echo "## Starting script 02_joint-geno_runner.sh"
date
echo
echo "## File ID:                   $FILE_ID"
echo "## Scaffold file:             $scaf_file"
echo "## Gvcf dir:                  $GVCF_DIR"
echo "## Vcf dir - by scaffold:     $VCF_DIR_SCAFFOLD"
echo "## Vcf dir - main:            $VCF_DIR_MAIN"
echo "## Vcf dir - final:           $VCF_DIR_FINAL"
echo "## QC dir:                    $QC_DIR"
echo "## Ref:                       $REF"
echo "## Additonal GATK commands:   $ADD_COMMANDS"
echo "## Min Mean DP:               $DP_MEAN"
echo
echo "## Memory allotted to jobs:   $MEM_JOB"
echo "## Memory for GATK:           $MEM_GATK"
echo "## Nr of cores:               $NCORES"
echo "## Skip genotyping:           $SKIP_GENO"
echo
echo "## Individuals: ${INDS[@]}"
echo

## Make dirs if needed:
mkdir -p "$VCF_DIR_SCAFFOLD"
mkdir -p "$VCF_DIR_FINAL"
mkdir -p "$SCAFFOLDLIST_DIR"


# JOINT GENOTYPING OF GVCF FILES - BY SCAFFOLD ---------------------------------
echo -e "\n-------------------------------"
if [ "$SKIP_GENO" == TRUE ]; then
	echo -e "## Skipping genotyping...\n"
else
	echo -e "## Calling joint genotyping script per (set of) scaffold(s)...\n"
	echo "## Checking for presence of VCFs ..."
	for IND in ${INDS[@]}; do
		echo -e "\n## Ind: $IND"
		ls -lh "$GVCF_DIR"/*"$IND"*vcf
	done

	MULTI_IND=TRUE
	LASTLINE=$(wc -l < "$scaf_file")
	echo -e "\n## Nr of scaffolds: $LASTLINE \n"
	
	for scaf_nr in $(seq 1 1 "$LASTLINE"); do
		INTERVAL_ID=$(head -n "$scaf_nr" "$scaf_file" | tail -n 1)
		INTERVAL_FILE=$SCAFFOLDLIST_DIR/$INTERVAL_ID.list
		head -n "$scaf_nr" "$scaf_file" | tail -n 1 > "$INTERVAL_FILE"
		
		echo "## Interval ID (scaffold nr): $INTERVAL_ID"
		echo "## Scaffold list:"
		cat "$INTERVAL_FILE"
				
		sbatch --mem "$MEM_JOB"G -N 1-1 --ntasks "$NCORES" --job-name=jgeno.pip."$FILE_ID" \
		$SCR_GENO "$FILE_ID" "$MULTI_IND" "$INTERVAL_FILE" "$INTERVAL_ID" \
		"$GVCF_DIR" "$VCF_DIR_SCAFFOLD" "$REF" "$ADD_COMMANDS" "$MEM_GATK" "$NCORES" ${INDS[@]}
	done
	
fi


# MERGE BY-SCAFFOLD VCFS -------------------------------------------------------
echo -e "\n-------------------------------"
echo -e "## Calling script to merge scaffolds...\n"

sbatch --mem "$MEM_JOB"G --job-name=jgeno.pip."$FILE_ID" --dependency=singleton \
	$SCR_MERGEVCFS "$FILE_ID" "$VCF_DIR_SCAFFOLD" "$VCF_DIR_MAIN"


# FILTER VCF FILES -------------------------------------------------------------
echo -e "\n-------------------------------"
echo -e "## Calling vcf filtering script...\n"

INPUT_NAME="$FILE_ID".rawvariants
OUTPUT_NAME="$FILE_ID"
MAC=3
FILTER_INDS_BY_MISSING=TRUE
SELECT_INDS_BY_FILE=FALSE
SAMPLE_ID_FILE=notany
INDSEL_ID=notany
JOBNAME="$FILE_ID"
SKIP_COMMON_STEPS="-456789tew"
SKIP_FINAL_STEPS="-123"
SKIP_IN_PIP=""

sbatch --mem "$MEM_JOB"G --job-name=jgeno.pip."$FILE_ID" --dependency=singleton \
	$SCR_FILTER "$INPUT_NAME" "$OUTPUT_NAME" "$VCF_DIR_MAIN" "$VCF_DIR_FINAL" \
    "$QC_DIR" "$REF" "$DP_MEAN" "$MAC" "$FILTER_INDS_BY_MISSING" "$SELECT_INDS_BY_FILE" \
    "$SAMPLE_ID_FILE" "$MEM_GATK" "$JOBNAME" "$INDSEL_ID" \
    "$SKIP_COMMON_STEPS" "$SKIP_FINAL_STEPS" "$SKIP_IN_PIP"

