#!/bin/bash

set -euo pipefail

echo "## Starting with script."
date
echo

# SET-UP -----------------------------------------------------------------------
## Software:
JAVA=software/java_1.8.0/jre1.8.0_144/bin/java
PICARD=software/picard_2.13.2/picard.jar
SAMTOOLS=software/samtools-1.6/samtools

## Scripts:
SCR_INDEXGENOME=scripts/conversionfasta_index.sh
SCR_REGIONSEL=scripts/conversion/bam_selectregions.sh

SCR_ALIGN=scripts/geno/bam1_align.sh
SCR_SORT=scripts/geno/bam2_sort.sh
SCR_FILTERBAM=scripts/geno/bam3_filter.sh
SCR_DEDUP=scripts/geno/bam4_dedup.sh
SCR_REALIGN=scripts/geno/bam5_realign.sh

SCR_VARDISC=scripts/geno/gatk1_vardisc.sh
SCR_GENO=scripts/geno/gatk2_jointgeno.sh

SCR_QCVCF=scripts/qc/qc_vcf.sh
SCR_QCBAM=scripts/qc/qc_bam.sh

## Command-line args:
ID_SHORT="$1"
shift
ID_LONG_FILE="$1"
shift
USE_R2="$1"
shift
REF="$1"
shift
FASTQ_DIR="$1"
shift
BAM_DIR="$1"
shift
VCF_DIR="$1"
shift
QC_DIR_VCF="$1"
shift
QC_DIR_BAM="$1"
shift
MINMAPQUAL="$1"
shift
DP_MEAN="$1"
shift
BAM_SUFFIX="$1" # When not running bam processing steps: Bam suffix of files to use for next step
shift
READGROUP_STRING="$1"
shift
REGION_FILE="$1" # bed file with regions to include or exclude from bam files
shift
REGION_SEL="$1"
shift
MEM="$1"
shift
NCORES="$1"
shift

MAP='true'
PROCESS_BAM_GENERAL='true'
SORT_BAM='true'
REALIGN_BAM='false' # Default for realignment is false!
DEDUP_BAM='true'
FILTER_BAM='true'
MERGE_BAMS='true'
USE_REGION_FILE_BAM='true'
QC_BAM='true'
VARDISC='true'
GENO='true'
FILTER_VCF='true'

while getopts 'ZAMPSRDQVGFfm' flag; do
  case "${flag}" in
    M) MAP='false' ;;
    m) MERGE_BAMS='false' ;;
    P) PROCESS_BAM_GENERAL='false' ;;
    S) SORT_BAM='false' ;;
    R) REALIGN_BAM='true' ;;
    D) DEDUP_BAM='false' ;;
	f) FILTER_BAM='false' ;;
    A) USE_REGION_FILE_BAM='false' ;;
    Q) QC_BAM='false' ;;
    V) VARDISC='false' ;;
    G) GENO='false' ;;
    F) FILTER_VCF='false' ;;
  esac
done

if [ $PROCESS_BAM_GENERAL == false ]; then
	REALIGN_BAM='false'
	DEDUP_BAM='false'
	SORT_BAM='false'
	USE_REGION_FILE_BAM='false'
	[[ $MAP == false ]] && QC_BAM='false'
fi

## Process:
if [ -s $ID_LONG_FILE ]; then
	echo -e "\n## Found IDS_LONG file, therefore reading IDS_LONG from file..."
	IDS_LONG=( $(cat $ID_LONG_FILE) )
else
	echo -e "\n## Did not find ID_LONG_FILE, therefore IDS_LONG = ID_SHORT"
	IDS_LONG=$ID_SHORT
fi

GVCF_DIR=$VCF_DIR/gvcf/
VCF_DIR_MAIN=$VCF_DIR/intermed/
VCF_DIR_FINAL=$VCF_DIR/final/

## Hardcoded:
INTERVAL_FILE=notany
INTERVAL_ID=allIntervals

## Report:
echo -e "\n## Indiv ID: $ID_SHORT"
echo -e "## Sample ID(s): ${IDS_LONG[@]}"
echo
echo "## Use R2: $USE_R2"
echo "## Ref: $REF"
echo "## Fastq dir: $FASTQ_DIR"
echo "## Bam dir: $BAM_DIR"
echo "## Gvcf dir: $GVCF_DIR"
echo "## Vcf dir main: $VCF_DIR_MAIN"
echo "## Vcf dir final: $VCF_DIR_FINAL"
echo "## Vcf QC dir: $QC_DIR_VCF"
echo "## Bam QC dir: $QC_DIR_BAM"
echo "## Minmapqual for samtools: $MINMAPQUAL"
echo "## Min-mean DP for filtering: $DP_MEAN"
echo "## Bam suffix: $BAM_SUFFIX"
echo "## Bam readgroup string: $READGROUP_STRING"
echo "## Memory: $MEM"
echo "## Nr of cores: $NCORES"
echo
echo "## To map: $MAP"
echo "## To process bam: $PROCESS_BAM_GENERAL"
echo "## To sort bam: $SORT_BAM"
echo "## To locally realign bam: $REALIGN_BAM"
echo "## To deduplicate bam: $DEDUP_BAM"
echo "## To use regions file to exclude non-autosomal seqs: $USE_REGION_FILE_BAM"
echo "## To QC bam: $QC_BAM"
echo "## To merge bams: $MERGE_BAMS"
echo "## To perform variant discovery: $VARDISC"
echo "## To perform (joint) genotyping: $GENO"
echo "## To filter and QC VCF: $FILTER_VCF"
echo

## Create directories, if needed:
mkdir -p $BAM_DIR/intermed $BAM_DIR/final $BAM_DIR/final_merged \
    $GVCF_DIR $VCF_DIR_MAIN $VCF_DIR_FINAL $QC_DIR_BAM $QC_DIR_VCF

## Fastq setting:
if [ $MAP == true ] && [ $PROCESS_BAM_GENERAL == true ]; then
    for ID_LONG in ${IDS_LONG[@]}; do
	    if [ -n "$(find $FASTQ_DIR -name "*$ID_LONG.*R0.*q.gz")" ]; then
	        echo -e "\n## NOTE: Single-end sequenced sample, setting USE_R2 to FALSE\n"
	        USE_R2=FALSE
	    fi
    done
fi

## Remove old bam index (bai) files:
echo -e "## Removing old bam index files...\n"
rm -f $BAM_DIR/intermed/$ID_SHORT*bai
rm -f $BAM_DIR/final/$ID_SHORT*bai
rm -f $BAM_DIR/final_merged/$ID_SHORT*bai


# PREP REFERENCE GENOME --------------------------------------------------------
## Index genome with bwa [Takes 1-2h, necessary for bwa mem mapping]:
if [ ! -f $REF.bwt ]; then
	echo -e "\n## Indexing ref genome with bwa..."
	$SCR_INDEXGENOME $REF
fi

## Create dictionary with Picard: [Takes 1 min, necessary for GATK genotyping]
if [ ! -f $REF.dict ]; then
	echo -e "\n## Creating ref genome dictionary with Picard..."
	$JAVA -Xmx4g -jar $PICARD CreateSequenceDictionary R=$REF O=$REF.dict
	OLDNAME=${REF}.dict
	NEWNAME="${OLDNAME/fasta.dict/dict}"
	NEWNAME="${NEWNAME/fna.dict/dict}"
	cp $OLDNAME $NEWNAME
fi

## Index fasta with samtools:
if [ ! -f $REF.fai ]; then
	echo -e "\n## Indexing ref genome with samtools..."
	$SAMTOOLS faidx $REF
fi

################################################################################
#### STEPS DONE SEPARATELY FOR EACH LIBRARY, PRIOR TO MERGING BY IND ####
################################################################################
for ID_LONG in ${IDS_LONG[@]}; do
	echo -e "\n###############################################################" 
	echo "## Processing $ID_LONG..."
		
	############################################################################
	#### ALIGN SEQUENCES ####
	############################################################################
	if [ $MAP == true ]; then
		echo -e "\n###########################################################"
		echo "## Calling mapping script for $ID_LONG..."
		
		## Get fastq files:
		FASTQ1=$(ls $FASTQ_DIR/*$ID_LONG.*q.gz | grep "\.R[0-1]\..*q.gz")
		[[ $USE_R2 == "TRUE" ]] && \
		  echo -e "## NOTE: Using R2 reads \n" && \
		  FASTQ2=$(ls $FASTQ_DIR/*$ID_LONG.R2.*q.gz) || \
		  FASTQ2=""
		
		## Run mapping script:
		$SCR_ALIGN $ID_LONG $READGROUP_STRING $REF $FASTQ_DIR $BAM_DIR/intermed $FASTQ1 $FASTQ2
		
	fi
	
	############################################################################
	#### PROCESS BAMFILES ####
	############################################################################
	SUFFIX_OUT=""
	
	if [ $PROCESS_BAM_GENERAL == true ]; then
		[[ $BAM_SUFFIX = "notany" ]] && SUFFIX_IN="" || SUFFIX_IN=$BAM_SUFFIX
		[[ $BAM_SUFFIX = "notany" ]] && SUFFIX_OUT="" || SUFFIX_OUT=$BAM_SUFFIX
		
		## Sort:
		if [ $SORT_BAM == true ]; then
			echo -e "\n#######################################################"
			echo "## Calling bam sorting script for $ID_LONG..."
			[[ -z $SUFFIX_OUT ]] && SUFFIX_OUT=sort || SUFFIX_OUT=$SUFFIX_OUT.sort
			echo "## Using bam file suffix: $SUFFIX_OUT"
			
			[[ -z $SUFFIX_IN ]] && INFILE=$BAM_DIR/intermed/$ID_LONG.bam || INFILE=$BAM_DIR/intermed/$ID_LONG.$SUFFIX_IN.bam
			OUTFILE=$BAM_DIR/intermed/$ID_LONG.$SUFFIX_OUT.bam
			MEM_BWA=$(( $MEM / $NCORES ))
			$SCR_SORT -c $NCORES -m $MEM_BWA $INFILE $OUTFILE
			
			SUFFIX_IN=$SUFFIX_OUT
		fi

		## Filter:
		if [ $FILTER_BAM == true ]; then
			echo -e "\n\n#######################################################"
			echo "## Calling bam filtering script for $ID_LONG..."
			SUFFIX_OUT=$SUFFIX_OUT.MQ$MINMAPQUAL
			echo "## Using bam file suffix: $SUFFIX_OUT"

			[[ $USE_R2 == "TRUE" ]] && echo -e "\n## Since USE_R2=TRUE, setting FILTER_PROPPAIRS to TRUE..." && FILTER_PROPPAIRS=TRUE
			[[ $USE_R2 == "FALSE" ]] && echo -e "\n## Since USE_R2=FALSE, setting FILTER_PROPPAIRS to FALSE..." && FILTER_PROPPAIRS=FALSE
			
			INFILE=$BAM_DIR/intermed/$ID_LONG.$SUFFIX_IN.bam
			OUTFILE=$BAM_DIR/intermed/$ID_LONG.$SUFFIX_OUT.bam
			$SCR_FILTERBAM -c $NCORES -f $FILTER_PROPPAIRS -q $MINMAPQUAL -s $QC_DIR_BAM $INFILE $OUTFILE

			SUFFIX_IN=$SUFFIX_OUT
		fi

		## Deduplicate if using R1 and R2:
		if [ $DEDUP_BAM == true ] && [ $USE_R2 == "TRUE" ]; then
		  echo -e "\n###################################################"
		  echo -e "## Calling bam deduplication script for $ID_LONG..."
		  SUFFIX_OUT=$SUFFIX_OUT.dedup
		  echo -e "## Using bam file SUFFIX: $SUFFIX_OUT \n"
		  
		  $SCR_DEDUP $ID_LONG $BAM_DIR/intermed $BAM_DIR/intermed $SUFFIX_IN $SUFFIX_OUT $QC_DIR_BAM

		  SUFFIX_IN=$SUFFIX_OUT
		fi
		
		## Realign:
		if [ $REALIGN_BAM == true ]; then
			echo -e "\n#######################################################"
			echo -e "## Calling bam realignment script for $ID_LONG..."
			SUFFIX_OUT=$SUFFIX_OUT.real
			echo -e "## Using bam file suffix: $SUFFIX_OUT \n"
			
			$SCR_REALIGN $ID_LONG $BAM_DIR/intermed $BAM_DIR/intermed $SUFFIX_IN $SUFFIX_OUT $REF $MEM
			
			SUFFIX_IN=$SUFFIX_OUT
		fi
		
		## Subset sequences:
		if [ $USE_REGION_FILE_BAM == true ]; then
			echo -e "\n\n#######################################################"
			echo -e "## Subsetting bam to specified regions for $ID_LONG..."
			echo -e "## Current bam file suffix: $SUFFIX_OUT \n"
			
			SUFFIX_OUT=$SUFFIX_OUT.auto # EDIT TO $REGION_ID
			REHEADER=FALSE
			PREFIX_REHEADER=notany
			
			$SCR_REGIONSEL $ID_LONG $SUFFIX_IN $SUFFIX_OUT $BAM_DIR/intermed $QC_DIR_BAM $REGION_FILE $REGION_SEL $REHEADER $PREFIX_REHEADER
		fi
		
		## Get stats on resulting bamfiles:
		if [ $QC_BAM == true ]; then
			echo -e "\n\n#######################################################"
			echo -e "## Calling bamstats script for $ID_LONG..."
			echo -e "## Current bam file suffix: $SUFFIX_OUT \n"
			
			INPUT=$BAM_DIR/intermed/$ID_LONG.$SUFFIX_OUT.bam
			UNSORTED=FALSE
			
			$SCR_QCBAM $ID_LONG $INPUT $QC_DIR_BAM $REF $MEM $UNSORTED
		fi
		
		## Move final bamfile to "final" directory:
		echo "## Moving final bamfile to final directory:"
		mv $BAM_DIR/intermed/$ID_LONG.$SUFFIX_OUT.bam $BAM_DIR/final/$ID_LONG.$SUFFIX_OUT.bam
		
		## List final bamfile:
		echo "## Final bamfile:"
		ls -lh $BAM_DIR/final/$ID_LONG.$SUFFIX_OUT.bam
		
		## Remove temp files:
		rm -r $BAM_DIR/intermed/$ID_LONG*bai
		[[ $SORT_BAM == "true" ]] && rm -f $BAM_DIR/intermed/$ID_LONG.sort.bam
		#[[ $SORT_BAM == "false" ]] && rm -f $BAM_DIR/intermed/$ID_LONG.bam
	fi
done


################################################################################
#### STEP 3b: MERGE BAMFILES ####
################################################################################
[[ -z $SUFFIX_OUT ]] && SUFFIX_OUT=$BAM_SUFFIX
BAM_OUT=$BAM_DIR/final_merged/$ID_SHORT.$SUFFIX_OUT.bam

if [ $MERGE_BAMS == true ]; then
	echo -e "\n\n###############################################################"
	echo -e "## Merging bam files...\n"
	
	## Prep:
	NR_FILES=${#IDS_LONG[@]}
	BAM_IN_COMMAND=$(for ID_LONG in ${IDS_LONG[@]}; do printf " $BAM_DIR/final/$ID_LONG.$SUFFIX_OUT.bam"; done)
	
	echo -e "## Current bam file suffix: $SUFFIX_OUT"
	echo -e "## Number of bam files for $ID_SHORT: $NR_FILES"
	echo -e "## Bam outfile: $BAM_OUT"
	echo -e "## Bam infile(s):\n$(echo $BAM_IN_COMMAND | tr " " "\n") \n"
	
	## If ID_SHORT has more than one ID_LONG, merge bams; otherwise, copy bam to final dir:
	if [ $NR_FILES -gt 1 ]
	then
		echo "## Calling samtools merge..."
		$SAMTOOLS merge -rf $BAM_OUT $BAM_IN_COMMAND
	else
		echo "## Individual $ID_SHORT does not have multiple IDs, simply copying bamfile..."
		cp -f $BAM_IN_COMMAND $BAM_OUT
	fi
	
	## Report:
	echo -e "\n## Final/merged bamfile:"
	ls -lh $BAM_OUT
	
	## QC:
	if [ $QC_BAM = 'true' ]; then
	  echo -e "\n## Calling bamstats script on final/merged bam...\n"
	  UNSORTED=FALSE
	  $SCR_QCBAM $ID_SHORT $BAM_OUT $QC_DIR_BAM $REF $MEM $UNSORTED
	fi
fi


################################################################################
#### STEP 4: VARIANT DISCOVERY WITH GATK ####
################################################################################
if [ $VARDISC == true ] ; then
	echo -e "\n\n###############################################################"
	echo "## Calling variant discovery script..."
	
	GATK_VERSION=gatk4
	INPUT=$BAM_OUT
	OUTPUT=$GVCF_DIR/$ID_SHORT.rawvariants.g.vcf
	
	echo -e "## Input: $INPUT"
	echo -e "## Output: $OUTPUT \n"
	
	$SCR_VARDISC $REF $INPUT $OUTPUT $MEM $NCORES $GATK_VERSION
fi


################################################################################
#### STEP 5: GENOTYPING OF GVCF FILES ####
################################################################################
if [ $GENO == true ] ; then
	echo -e "\n\n###############################################################"
	echo -e "## Calling (joint) genotyping script...\n"
	
	ADD_COMMANDS="none"
	MULTI_IND=FALSE
	SETNAME=$ID_SHORT
	
	$SCR_GENO $SETNAME $MULTI_IND $INTERVAL_FILE $INTERVAL_ID $GVCF_DIR $VCF_DIR_MAIN \
		$REF "$ADD_COMMANDS" $MEM $NCORES $ID_SHORT
	
	echo -e "\n## QC of raw VCF...\n"
	$SCR_QCVCF $ID_SHORT.$INTERVAL_ID.rawvariants $VCF_DIR_MAIN $QC_DIR_VCF FALSE FALSE FALSE
fi


################################################################################
#### STEP 6: FILTER VCF FILES ####
################################################################################
if [ $FILTER_VCF == true ]; then
	echo -e "\n\n###############################################################"
	echo -e "## Calling vcf filtering script...\n"
	
	INPUT_NAME=$ID_SHORT.$INTERVAL_ID.rawvariants
	OUTPUT_NAME=$ID_SHORT
	
	## Hardcoded options:
	MAC=1
	INDFILE=notany
	FILTER_INDS_BY_MISSING=TRUE
	SELECT_INDS_BY_FILE=FALSE
	
	## Run filtering script:
	$SCR_FILTER $INPUT_NAME $OUTPUT_NAME $VCF_DIR_MAIN $VCF_DIR_FINAL $QC_DIR_VCF \
		$MEM $REF $INDFILE $DP_MEAN $MAC $FILTER_INDS_BY_MISSING $SELECT_INDS_BY_FILE
	
	echo -e "\n## Removing intermediate VCF files...\n"
	ls -lh $VCF_DIR_MAIN/$ID_SHORT*
	rm -f $VCF_DIR_MAIN/$ID_SHORT*
fi

echo -e "\n## Done with script."
date
echo
