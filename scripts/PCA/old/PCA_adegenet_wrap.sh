#!/bin/bash
set -e
set -o pipefail
set -u

module load R

## Command-line args:
FILE_ID=$1
BASEDIR=$2
ID_TYPE=$3
KEEP_SET=$4
FILE_INDS_KEEP=$5
FILE_LOOKUP=$6
FILE_COLS=$7
VCF_DIR=$8

## Report:
date
echo "##### PCA_adegenet_submit.sh: Starting script."

## Run R-script:
Rscript scripts/PCA/PCA_adegenet_cluster.R $FILE_ID $BASEDIR $ID_TYPE $KEEP_SET $FILE_INDS_KEEP $FILE_LOOKUP $FILE_COLS $VCF_DIR

## Report:
echo "##### PCA_adegenet_submit.sh: Done with script."
date