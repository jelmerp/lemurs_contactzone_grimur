#!/bin/bash
set -e
set -o pipefail
#set -u # Need flexibility to have $FASTQ2 unboundq

################################################################################
#### SET-UP ####
################################################################################
## Command-line args:
ID=$1
READGROUP_STRING="$2" #"@RG\tID:group1\tSM:$IND\tPL:illumina\tLB:lib1"
REF=$3
FASTQ_DIR=$4
BAM_DIR=$5
FASTQ1=$6
FASTQ2=$7

## Process args:
[[ ! -d $BAM_DIR ]] && mkdir -p $BAM_DIR
BAM_OUT=$BAM_DIR/$ID.bam

## Report:
echo
date
echo "## bam1_align.sh: Script: bam1_align.sh"
echo "## bam1_align.sh: Job ID: $SLURM_JOB_ID"
echo "## bam1_align.sh: Number of nodes: $SLURM_JOB_NUM_NODES"
echo "## bam1_align.sh: Nr of tasks: $SLURM_NTASKS"
echo
echo "## bam1_align.sh: Aligning fastq files to reference genome for: $ID"
echo "## bam1_align.sh: Readgroup string: $READGROUP_STRING"
echo "## bam1_align.sh: Reference fasta: $REF"
echo "## bam1_align.sh: Fastq 1: $FASTQ1"
echo "## bam1_align.sh: Fastq 2: $FASTQ2"
echo "## bam1_align.sh: Bam (output) dir: $BAM_DIR"
echo "## bam1_align.sh: Bam (output) file: $BAM_OUT"


################################################################################
#### MAP WITH BWA ####
################################################################################
echo -e "\n## bam1_align.sh: Mapping with bwa mem..."
bwa mem $REF -aM -R "$READGROUP_STRING" -t $SLURM_NTASKS $FASTQ1 $FASTQ2 | samtools view -b -h > $BAM_OUT

## Report:
echo -e "\n## bam1_align.sh: Resulting bam file:"
ls -lh $BAM_OUT

date
echo -e "\n## bam1_align.sh: Done with script.\n"


################################################################################
################################################################################
## bwa flags:
# -t nr of threads
# -a alignments for single-end / unpaired reads are also output, as secondary alignments
# -M shorter split reads are output as secondary alignments, for Picard compatibility
# -R "@RG\tID:group1\tSM:$IND\tPL:illumina\tLB:lib1"