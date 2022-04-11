#!/bin/bash

set -euo pipefail

# SETUP ------------------------------------------------------------------------
## Software and scripts
VCFTOOLS=software/vcftools/vcftools-master/bin/vcftools
module load Plink/1.90

## Positional args:
FILE_ID=$1
shift
VCF_DIR=$1
shift
PLINK_DIR=$1
shift

## Optional args:
MAF=0                   # Minor allele frequency; give 0 for no filtering
LD_MAX=1                # Give 1 if no LD pruning needed
INDFILE=FALSE
ID_OUT=$FILE_ID

while getopts 'm:l:i:o:' flag; do
  case "${flag}" in
  	m) MAF="${OPTARG}" ;;
	l) LD_MAX="${OPTARG}" ;;
    i) INDFILE="${OPTARG}" ;;
    o) ID_OUT="${OPTARG}" ;;
  esac
done

## Process args:
OUTFILE=$PLINK_DIR/$ID_OUT
VCF=$VCF_DIR/$FILE_ID.vcf.gz
[[ ! -d $PLINK_DIR ]] && echo "#### Creating dir $PLINK_DIR" && mkdir -p $PLINK_DIR

## Report:
echo "## Starting script."
date
echo
echo "## Positional args:"
echo "## File ID:                                $FILE_ID"
echo "## Vcf dir:                                $VCF_DIR"
echo "## Plink dir:                              $PLINK_DIR"
echo
echo "## Optional args:"
echo "## Minor allele frequency (MAF) cut-off:   $MAF"
echo "## LD_MAX:                                 $LD_MAX"
echo "## Indfile (subset inds to list):          $INDFILE"
echo "## Output ID:                              $ID_OUT"
echo
echo "## Processed args:"
echo "## VCF file:                               $VCF"
echo "## Output file:                            $OUTFILE"
echo -e "-----------------------\n"


# PREP COMMAND TO LET VCFTOOLS SUBSET INDIVIDUALS (IF NEEDED) ------------------
## If $INDFILE exists (-e), then assign a "keep command"
## This command will be passed on to the vcftools program when it does the vcf->plink conversion
## If $KEEP_COMMAND is not assigned here, it can still be passed on to vcftools (see below),
## but since the variable is blank, nothing will be processed.  
if [ $INDFILE != FALSE ]; then

	[[ ! -f $INDFILE ]] && echo "Indfile $INDFILE DOES NOT EXIST" && exit 1
	echo "## Selecting individuals from file: $INDFILE"
	KEEP_COMMAND="--keep $INDFILE"
	
	echo -e "## Keep command: $KEEP_COMMAND \n"

else

	echo "## No ind-file specified - keeping all individuals."
	KEEP_COMMAND=""

fi


# CONVERT VCF TO PLINK ---------------------------------------------------------
echo -e "\n\n## Converting vcf to plink..."
$VCFTOOLS --gzvcf $VCF $KEEP_COMMAND --plink --maf $MAF --out $OUTFILE


# EDIT PLINK FILES -------------------------------------------------------------
## Replace chromosome notations
echo "## Replacing chrom 0 by chrom 1 in PLINK map file..."
sed 's/^0/1/g' $OUTFILE.map > $OUTFILE.tmp.map
mv $OUTFILE.tmp.map $OUTFILE.map

## Create PLINK bed file
echo -e "\n\n## Creating binary PLINK files..."
plink --file $OUTFILE --make-bed --out $OUTFILE

echo -e "\n\n## Creating PLINK files for adegenet..."
plink --file $OUTFILE --recodeA --out $OUTFILE.recodeA

## Perform LD pruning in PLINK, if $LD_MAX is not 1:
if [ $LD_MAX != 1 ]; then

	echo -e "\n\n## LD pruning with PLINK..."
	mkdir $ID_OUT.$LD_MAX
	cd $ID_OUT.$LD_MAX
	plink --file ../$OUTFILE --indep-pairwise 50 5 $LD_MAX # calculate LD
	plink --file ../$OUTFILE --extract plink.prune.in --make-bed --out ../$PLINK_DIR/plink/$ID_OUT.LDpruned$LD_MAX # prune high LD sites
	plink --bfile ../$OUTFILE.LDpruned$LD_MAX --recode --out ../$OUTFILE.LDpruned$LD_MAX # convert binary "bed" files back to ped
	cd ..
	rm -r $ID_OUT.$LD_MAX

else

	echo -e "\n\n## Skipping LD pruning..."

fi


# WRAP UP ----------------------------------------------------------------------
echo -e "\n\n## Output files:"
ls -lh $OUTFILE*
echo -e "\n## Done with script vcf2plink.sh"
date
echo
