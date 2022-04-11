#!/bin/bash

set -euo pipefail

# SET-UP -----------------------------------------------------------------------
## Software
JAVA=software/java_1.8.0/jre1.8.0_144/bin/java
PICARD=software/picard_2.13.2/picard.jar

## Command-line arguments
SETNAME=$1
VCF_DIR_IN=$2
VCF_DIR_OUT=$3

## Process args
VCF_LIST=$VCF_DIR_IN/$SETNAME.vcflist.txt
ls "$VCF_DIR_IN"/"$SETNAME".*rawvariants.vcf > "$VCF_LIST"
VCF_OUT="$VCF_DIR_OUT"/"$SETNAME".rawvariants.vcf

## Make output dir
mkdir -p "$VCF_DIR_OUT" 

## Report
echo "## Starting script gatk3_mergevcf.sh"
date
echo
echo "## Set name:       $SETNAME"
echo "## VCF dir in:     $VCF_DIR_IN"
echo "## VCF dir out:    $VCF_DIR_OUT"
echo
echo "## VCF list:       $VCF_LIST"
echo "## Output VCF:     $VCF_OUT"


# RUN GATK MERGEVCFS -----------------------------------------------------------
echo -e "\n## Running mergevcfs...\n"
$JAVA -jar $PICARD MergeVcfs I="$VCF_LIST" O="$VCF_OUT"


# HOUSEKEEPING -----------------------------------------------------------------
echo -e "\n## Output VCF $VCF_OUT"
ls -lh "$VCF_OUT"
NVAR=$(grep -vc "##" "$VCF_OUT")
echo -e "\n## Number of variants: $NVAR"
echo -e "\n## Done with script."
date
echo