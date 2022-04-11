#!/bin/bash

set -euo pipefail

# SET-UP -----------------------------------------------------------------------
## Software and hard-coded dirs
GATK4_EXC=software/gatk-4.0.7.0/gatk
TMP_DIR=/work/jwp37/javaTmpDir

## Command-line args:
SETNAME="$1"
shift
MULTI_IND="$1"
shift
INTERVAL_FILE="$1"
shift
INTERVAL_ID="$1"
shift
GVCF_DIR="$1"
shift
VCF_DIR_SCAFFOLD="$1"
shift
REF="$1"
shift
ADD_COMMANDS="$1"
shift
MEM="$1"
shift
NCORES="$1"
shift
count=0

while [ "$*" != "" ]; do
    INDS[$count]=$1
    shift
    count=4(expr $count + 1)
done

## Process command-line args:
VCF_OUT=$VCF_DIR_SCAFFOLD/$SETNAME.$INTERVAL_ID.rawvariants.vcf
[[ $ADD_COMMANDS == "none" ]] && ADD_COMMANDS=""
INDS_COMMAND=$(for IND in ${INDS[@]}; do printf " --variant $GVCF_DIR/$IND.rawvariants.g.vcf"; done)

mkdir -p $GVCF_DIR/byScaffold
mkdir -p $TMP_DIR

## Report:
echo -e "\n## Starting script gatk2_jointgeno.sh"
date
echo
echo "## Set name: $SETNAME"
echo "## Multi-ind TRUE/FALSE: $MULTI_IND"
echo "## Interval file: $INTERVAL_FILE"
echo "## Interval ID: $INTERVAL_ID"
echo "## Input dir: $GVCF_DIR"
echo "## Output dir: $VCF_DIR_SCAFFOLD"
echo "## Reference genome file: $REF"
echo "## Additional GATK commands: $ADD_COMMANDS"
echo "## Memory: $MEM"
echo "## Number of cores: $NCORES"
echo
echo "## Number of individuals: ${#INDS[@]}"
echo "## Individuals: ${INDS[@]}"
echo
echo "## Output VCF: $VCF_OUT"


################################################################################
#### MERGE GVCFS BY IND (FOR ONE SCAFFOLD) ####
###############################################################################
echo -e "\n\n###################################################################"

if [ $MULTI_IND = TRUE ]; then
	
    echo -e "## Step 1: running GenomicsDBImport...\n"
	
	DB_DIR=$VCF_DIR_SCAFFOLD/DBs/$SETNAME.$INTERVAL_ID
	[[ -d $DB_DIR ]] && echo "## Removing directory $DB_DIR" && rm -r $DB_DIR # GATK throws error if dir already exists.
	[[ ! -d $VCF_DIR_SCAFFOLD/DBs ]] && echo "## Creating directory $VCF_DIR_SCAFFOLD/DBs" && mkdir -p $VCF_DIR_SCAFFOLD/DBs
	
	echo -e "## Database dir: $DB_DIR"
	echo -e "## Output VCF file: $VCF_OUT \n"
	
	$GATK4_EXC --java-options "-Xmx${MEM}g" GenomicsDBImport \
		$INDS_COMMAND \
		--genomicsdb-workspace-path $DB_DIR \
		--batch-size 0 \
		--intervals $INTERVAL_FILE \
		--reader-threads $NCORES \
		--interval-padding 100
	
	GENO_INPUT="gendb://$DB_DIR"
else
	echo "## Analyzing single sample, so skipping GenomicsDBImport..."
	GENO_INPUT=$GVCF_DIR/$SETNAME.rawvariants.g.vcf
fi


################################################################################
#### RUN GATK JOINT GENOTYPING ####
################################################################################
echo -e "\n\n#################################################################"
echo -e "## Step 2: genotyping GVCF with all inds combined..."
echo -e "## GVCF input: $GENO_INPUT"
echo -e "## Output VCF file: $VCF_OUT \n"

## Run GATK:
$GATK4_EXC --java-options "-Xmx${MEM}g" GenotypeGVCFs \
	-R $REF -V $GENO_INPUT --use-new-qual-calculator \
	-O $VCF_OUT --TMP_DIR=$TMP_DIR

#-G StandardAnnotation


################################################################################
#### HOUSEKEEPING ####
################################################################################
echo -e "\n## Output file:"
ls -lh $VCF_OUT

NVAR=$(grep -v "##" $VCF_OUT | wc -l)
echo -e "\n## Number of variants in jointly-genotyped-vcf: $NVAR \n"

echo "## Done with script gatk2_jointgeno.sh"
date