#!/bin/bash

set -euo pipefail

# SET-UP -----------------------------------------------------------------------

## Scripts:
SCRIPT_VCF2PLINK=scripts/conversion/vcf2plink.sh
SCRIPT_ADMIX_RUN=scripts/admixture/admixture_run.sh

## Positional args:
ID_IN=$1
shift
VCF_DIR=$1
shift
PLINK_DIR=$1
shift
OUTDIR=$1
shift

## Optional args:
MAF=0
LD_MAX=1
NCORES=1
INDFILE=FALSE
ID_OUT=$ID_IN
K_MIN=1
K_MAX=9

while getopts 'm:l:c:i:o:k:K:' flag; do
    case "${flag}" in
    m) MAF="${OPTARG}" ;;
    l) LD_MAX="${OPTARG}" ;;
    c) NCORES="${OPTARG}" ;;
    i) INDFILE="${OPTARG}" ;;
    o) ID_OUT="${OPTARG}" ;;
    k) K_MIN="${OPTARG}" ;;
    K) K_MAX="${OPTARG}" ;;
    \?) echo "Error: Invalid option" && exit 1 ;;
    esac
done

## Report:
date
echo "## Starting script."
printf "\n"
echo "## Positional args:"
echo "## File ID: $ID_IN"
echo "## VCF dir: $VCF_DIR"
echo "## Plink dir: $PLINK_DIR"
echo "## Output dir: $OUTDIR"
printf "\n"
echo "## Optional args:"
echo "## MAF: $MAF"
echo "## Max LD: $LD_MAX"
echo "## Number of cores: $NCORES"
echo "## File with individuals to subset: $INDFILE"
echo "## Output ID: $ID_OUT"
echo "## Min value of K: $K_MIN"
echo "## Max value of K: $K_MAX"
printf "\n"

[[ ! -d $PLINK_DIR ]] && mkdir -p "$PLINK_DIR" # If absent, create dir for PLINK-files
[[ ! -d $OUTDIR ]] && mkdir -p "$OUTDIR"       # If absent, create dir for ADMIXTURE output

# CONVERT VCF TO PLINK ---------------------------------------------------------

echo -e "\n-------------------"
echo "## Running VCF2PLINK script..."
$SCRIPT_VCF2PLINK "$ID_IN" "$VCF_DIR" "$PLINK_DIR" -m "$MAF" -l "$LD_MAX" -i "$INDFILE" -o "$ID_OUT"

# INDIVIDUALS  -----------------------------------------------------------------

echo -e "\n-------------------"
echo "## Creating individuals-file..."

INDIV_FILE=$OUTDIR/$ID_OUT.indivs.txt

PEDFILE=$PLINK_DIR/$ID_OUT.ped

cut -f1 "$PEDFILE" >"$INDIV_FILE"

echo "## Ped file: $PEDFILE"
echo "## Listing indiv file:"
ls -lh "$INDIV_FILE"

# RUN ADMIXTURE FOR EACH K
echo -e "\n-------------------"
echo "## Submitting ADMIXTURE script..."

for K in $(seq "$K_MIN" "$K_MAX"); do
    echo -e "\n## Value of K: $K"
    #sbatch -o slurm.admixture.run."$ID_OUT".K"$K" \
    $SCRIPT_ADMIX_RUN "$ID_OUT" "$PLINK_DIR" "$OUTDIR" "$K"
done

## Report
echo -e "\n-------------------"
echo "## Done with script 01_admixture_runner.sh"
date
