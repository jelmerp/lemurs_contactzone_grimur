#!/bin/bash

set -euo pipefail

# SETUP ------------------------------------------------------------------------
## Software
BCFTOOLS=software/bcftools-1.10.2/bcftools
VCF2PLINK_SCRIPT=scripts/conversion/vcf2plink.sh
MAKE_INDFILE_SCRIPT=scripts/admixtools/admixtools_makeIndfile.R

## Command-line args
FILE_ID=$1
shift
VCF_DIR=$1
shift
PLINK_DIR=$1
shift
VCF2PLINK=$1
shift
CREATE_INDFILE=$1
shift
SUBSET_INDFILE=$1
shift
INDFILE=$1
shift
POPFILE=$1
shift
PARFILE=$1
shift
ATOOLS_MODE=$1
shift
INDS_METADATA=$1
shift
ID_COLUMN=$1
shift
GROUPBY=$1
shift

## Process parameters
VCF=$VCF_DIR/$FILE_ID.vcf.gz # input file
PEDFILE=$PLINK_DIR/$FILE_ID.ped 
MAPFILE=$PLINK_DIR/$FILE_ID.map

## Report
echo "## Starting with script."
date
echo
echo "## File ID: $FILE_ID"
echo "## VCF dir: $VCF_DIR"
echo "## PLINK dir: $PLINK_DIR"
echo "## Create indfile (TRUE/FALSE): $CREATE_INDFILE"
echo "## Subset indfile (TRUE/FALSE): $SUBSET_INDFILE"
printf "\n"
echo "## VCF file (input): $VCF"
echo "## Popfile (input): $POPFILE"
echo "## Indfile (output): $INDFILE"
echo "## Parfile (output): $PARFILE"
echo "## PED file (output): $PEDFILE"
echo "## MAP file (output): $MAPFILE"
printf "\n"
echo "## Admixtools mode: $ATOOLS_MODE"
printf "\n"
echo "## Inds metadata (optional): $INDS_METADATA"
echo "## ID column: $ID_COLUMN"
echo "## Group-by column: $GROUPBY"
printf "\n"


################################################################################
#### CONVERT VCF TO PLINK FORMAT ####
################################################################################
if [ $VCF2PLINK == TRUE ]
then
	echo "## Converting vcf to plink..."
	MAF=0
	LD_MAX=1
	SELECT_INDS=FALSE
	INDFILE_PLINK="NA"
	ID_OUT="NA"
	$VCF2PLINK_SCRIPT $FILE_ID $VCF_DIR $PLINK_DIR $MAF $LD_MAX $SELECT_INDS $INDFILE_PLINK $ID_OUT
else
	echo -e "\n## Not converting vcf to plink.\n"
fi


################################################################################
#### CREATE EIGENSTAT INDFILE ####
################################################################################
if [ $CREATE_INDFILE == TRUE ]
then
	echo "## Creating eigenstat indfile..."
	echo "## Inds metadata: $INDS_METADATA"
	echo "## Indfile: $INDFILE"
	
	INDLIST=indlist.$FILE_ID.tmp
	$BCFTOOLS query -l $VCF > $INDLIST
	
	Rscript $MAKE_INDFILE_SCRIPT $INDLIST $INDS_METADATA $INDFILE $ID_COLUMN $GROUPBY
	
	rm $INDLIST
else
	echo -e "\n## NOT CREATING EIGENSTAT INDFILE.\n"
fi


################################################################################
#### SUBSET INDFILE ####
################################################################################
## (INDFILE CAN ONLY LIST INDS THAT ARE ACTUALLY PRESENT IN THE VCF FILE)

if [ $SUBSET_INDFILE == TRUE ]
then
	echo -e "\n\n## Subsetting Eigenstat indfile..."
	NLINE_IN=$(cat $INDFILE | wc -l)
	
	INDS=( $($BCFTOOLS query -l $VCF) )
	> $INDFILE.tmp
	for IND in ${INDS[@]}; do grep $IND $INDFILE >> $INDFILE.tmp; done
	sort -u $INDFILE.tmp > $INDFILE
	rm $INDFILE.tmp
	
	NLINE_OUT=$(cat $INDFILE | wc -l)
	echo -e "## Nr of lines before subsetting: $NLINE_IN"
	echo -e "## Nr of lines after subsetting: $NLINE_OUT \n"
	[[ $NLINE_OUT == 0 ]] && echo -e "\n\n\n\n## ERROR: NO LINES LEFT IN INDFILE!\n\n\n\n"
fi


################################################################################
#### CREATE EIGENSTAT PARFILES #####
################################################################################
echo -e "\n#####################################################################"
echo "## Creating Eigenstat parfile..."

if [ $ATOOLS_MODE == "D" ]
then
	echo -e "\n## Creating parfile for D-mode...\n"
		
	printf "genotypename:\t$PEDFILE\n" > $PARFILE
	printf "snpname:\t$MAPFILE\n" >> $PARFILE
	printf "indivname:\t$INDFILE\n" >> $PARFILE
	printf "popfilename:\t$POPFILE\n" >> $PARFILE
	printf "printsd:\tYES\n" >> $PARFILE
	
	echo "## Parfile $PARFILE:"
	cat $PARFILE
	printf "\n"
	
elif [ $ATOOLS_MODE == F4 ]
then
	echo -e "\n#### Creating parfile for F4-mode...\n"
		
	#cp $PARFILE_DMODE $PARFILE_FMODE ### EDIT!
	printf "f4mode:\tYES\n" >> $PARFILE_FMODE
	
	echo -e "\n## Parfile $PARFILE:"
	cat $PARFILE
	printf "\n"
	
elif [ $ATOOLS_MODE == F3 ]
then
	echo "\n## Creating parfile for f3-mode...\n"
	
	printf "genotypename:\t$PEDFILE\n" > $PARFILE
	printf "snpname:\t$MAPFILE\n" >> $PARFILE
	printf "indivname:\t$INDFILE\n" >> $PARFILE
	printf "popfilename:\t$POPFILE\n" >> $PARFILE
	#printf "printsd:\tYES\n" >> $PARFILE
	
	echo "## Parfile $PARFILE:"
	cat $PARFILE
	printf "\n"
	
elif [ $ATOOLS_MODE == F4RATIO ]
then
	echo -e "\n## Creating parfile for f4-ratio-mode...\n"
	echo "## Creating parfile $PARFILE..."
	
	printf "genotypename:\t$PEDFILE\n" > $PARFILE
	printf "snpname:\t$MAPFILE\n" >> $PARFILE
	printf "indivname:\t$INDFILE\n" >> $PARFILE
	printf "popfilename:\t$POPFILE\n" >> $PARFILE
	printf "printsd:\tYES\n" >> $PARFILE
	
	echo "## Parfile $PARFILE:"
	cat $PARFILE
	printf "\n"
	
else
	echo -e "\n\n\n## ATOOLS_MODE variable $ATOOLS_MODE does not match any mode..."
	echo -e "## NOT CREATING PARFILE...\n\n\n"
fi


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Done with script.\n"
date
echo
