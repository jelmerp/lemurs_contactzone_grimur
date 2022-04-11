#!/bin/bash

set -euo pipefail

# SET-UP -----------------------------------------------------------------------
## Software
SAMTOOLS=software/samtools-1.6/samtools
GATK3=software/gatk-3.8-0/GenomeAnalysisTK.jar
JAVA=software/java_1.8.0/jre1.8.0_144/bin/java
GATK4=software/gatk-4.0.7.0/gatk

## Command line args
ref=$1
input=$2
output=$3
mem=$4
ncores=$5
gatk_version=$6

## Report:
echo "## Starting script gatk1_vardisc.sh"
date
echo
echo "## Reference sequence:   $ref"
echo "## Input file name:      $input"
echo "## Output file name:     $output"
echo "## GATK version:         $gatk_version"
echo "## Memory:               $mem"
echo "## Number of cores:      $ncores"
echo

## Index bamfile:
if [ ! -e "$input".bai ]
then
	echo -e "\n## Indexing bamfile..." 
	$SAMTOOLS index "$input"
fi


# RUN GATK HAPLOTYPECALLER -----------------------------------------------------
echo -e "\n## Starting variant discovery...\n"

if [ "$gatk_version" = "gatk3" ]; then
	
    echo -e "## Running with GATK version 3... \n"
	"$JAVA" -Xmx"$mem"G -jar "$GATK3" -T HaplotypeCaller \
        -R "$ref" -I "$input" -o "$output" \
        --genotyping_mode DISCOVERY \
	    --emitRefConfidence GVCF -mmq 10 -nct "$ncores" 

elif [ "$gatk_version" = "gatk4" ]; then
	
    echo -e "## Running with GATK version 4...\n"
	"$GATK4" --java-options "-Xmx"$mem"g" HaplotypeCaller \
        -R "$ref" -I "$input" -O "$output" \
	    -ERC GVCF --pairHMM AVX_LOGLESS_CACHING_OMP \
        --native-pair-hmm-threads "$ncores"
fi


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Output file:"
ls -lh "$output"

NVAR=$(grep -vc "##" "$output")
echo -e "\n## Number of variants in gvcf: $NVAR \n"

echo -e "## Done with script gatk1_vardish.sh"
date
echo
