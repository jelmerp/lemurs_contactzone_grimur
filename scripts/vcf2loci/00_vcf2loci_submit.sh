#!/bin/bash

## Scripts
SCRIPT_LOCUSSTATS=scripts/vcf2loci/03e_locusstats.sh
SCRIPT_RUNNER=scripts/vcf2loci/vcf2loci_runner.sh

## Individuals & IDs:
FILE_ID=hzproj1 # Basic file ID of VCF / dataset 
SUBSET_ID="hz.mur3gri2c"
INDFILE=results/gphocs/metadata/gphocs_inds.txt

## Reference genome:
REF=data/refgenome/GCF_000165445.2_Mmur_3.0_genomic_stitched.fasta

## Settings for defining and filtering loci:
MIN_ELEMENT_OVERLAP=0.9              # Number of bedfile elements that should overlap to create a locus, given as % of nr of samples.
MIN_ELEMENT_OVERLAP_TRIM=0.8         # Number of bedfile elements that should overlap *at each basepair-level position* at the edges of each locus, given as % of nr of samples: lower coverage ends are trimmed.
MIN_LOCUS_SIZE=100                   # Minimum locus size in bp, smaller loci are not retained.
VCF2FULLFASTA_ID=ov$MIN_ELEMENT_OVERLAP.ovt$MIN_ELEMENT_OVERLAP_TRIM.ls$MIN_LOCUS_SIZE # File suffix to indicate locus production parameters
TRESHOLD_MISSING=10                  # Missing data treshold for final locus selection (percentage)

## Settings for VCF filtering:
VCF_RAWFILE_SUFFIX=rawSNPs.ABHet     # Suffix for input VCF file, vcf file should be $VCF_DIR_MAIN/$FILE_ID.$VCF_RAWFILE_SUFFIX.vcf(.gz)
VCF_DP_MEAN=15                       # Minimum mean-depth (across samples) per site -- lower will be filtered
VCF_MAC=1                            # Minimum Minor Allele Count (MAC). A file without and with MAC filtering will be produced.
SELECT_INDS_BY_FILE=TRUE # Whether or not to (pre)select individuals/samples using a file with IDs (TRUE/FALSE)
FILTER_INDS_BY_MISSING=FALSE # Whether or not to remove individuals/samples with high amounts of missing data

## Settings for CallableLoci:
BAMFILE_SUFFIX=sort.MQ30
CALLABLE_COMMAND="--minDepth 3" # (Added) Command for GATK CallableLoci: use to indicate minimum depth (and can also be used for max depth, etc)
CALLABLE_ID=DP3 # File suffix to indicate CallableLoci settings
CALLABLE_BEDFILE_SUFFIX=callable # Suffix for bedfile with callable sites produced by CallableLoci

## Directories:
VCF_DIR_MAIN=results/geno/vcf/joint_intermed/ # Dir with raw and intermediate VCF files
VCF_DIR_FINAL=results/geno/vcf # Dir with final filtered VCF file
VCF_QC_DIR=results/geno/qc/vcf_joint # Dir for VCF QC output
BAM_DIR=results/bam/final_merged/ # Dir with pre-existing bamfiles
FASTA_DIR=results/fasta_full/ # Dir with fasta files to produce
CREATELOCI_DIR=results/fasta_fill/createLoci/ # Dir with e.g. bed files to produce

## Steps to skip in each script:
SKIP_IN_PIP="" #-FB12" #"F" #"-B1" # What to skip in pipeline script: -F: VCF file filtering / -B: skip bedfile creation / -1: skip vcf2fullFasta1.sh / -2: skip vcf2fullFasta1.sh
SKIP_IN_PIP_VCF_FILTER="-H"
SKIP_COMMON_STEPS_VCF_FILTER="-12456789tew"
SKIP_IN_VCF2FULLFASTA1="" # What to skip in vcf2fullFasta1: -C: skip CallableLoci step / -A: skip Altref step / -M: skip mask step
SKIP_IN_VCF2FULLFASTA2="" #"-CIX" # What to skip in vcf2fullFasta2: -C: skip create-loci step / -I: skip intersect-with-VCF step / -X: skip extract-loci step / -F: skip filter-loci step

## Memory:
MEM=14 # Memory (GB) to use for GATK

## Run main script
"$SCRIPT_RUNNER" \
    $INDFILE $FILE_ID "$SUBSET_ID" $VCF2FULLFASTA_ID $REF \
    $VCF_RAWFILE_SUFFIX $VCF_DP_MEAN $VCF_MAC $SELECT_INDS_BY_FILE $FILTER_INDS_BY_MISSING \
    $MIN_ELEMENT_OVERLAP $MIN_ELEMENT_OVERLAP_TRIM $MIN_LOCUS_SIZE $TRESHOLD_MISSING \
    "$CALLABLE_COMMAND" $CALLABLE_ID $CALLABLE_BEDFILE_SUFFIX $BAMFILE_SUFFIX \
    $VCF_DIR_MAIN $VCF_DIR_FINAL $VCF_QC_DIR $BAM_DIR $FASTA_DIR $CREATELOCI_DIR $MEM \
    $SKIP_COMMON_STEPS_VCF_FILTER "$SKIP_IN_PIP_VCF_FILTER" "$SKIP_IN_VCF2FULLFASTA1" \
    "$SKIP_IN_VCF2FULLFASTA2" "$SKIP_IN_PIP"

## Get locus stats
DIR_FASTA=results/fasta_full/byLocus.final.r03.wOutgroups.hz.mur3gri2c.mac1.FS7.callableDP3.ov0.9.ovt0.8.ls100
FILE_STATS_ALL=mur3gri2c.locusstats.txt
sbatch "$SCRIPT_LOCUSSTATS" $DIR_FASTA $FILE_STATS_ALL
