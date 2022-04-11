FASTQ_DIR=data/fastq/raw
CUTSITE="TGCAGG"

################################################################################
##### L1_redo #####
################################################################################
LIB=L1_redo_S2

FASTQ=$FASTQ_DIR/$LIB*R1*gz
OUTFILE=results/qc/fastq/barcodeCounts_$LIB.R1.txt
sbatch scripts/fastq_process/checkBarcodes.sh $FASTQ $CUTSITE $OUTFILE

FASTQ=$FASTQ_DIR/$LIB*R2*gz
OUTFILE=results/qc/fastq/barcodeCounts_$LIB.R2.txt
sbatch scripts/fastq_process/checkBarcodes.sh $FASTQ $CUTSITE $OUTFILE

## Check adapters:
ADAPTER="GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTACAAGATCTCGTAT"
zgrep --color=always $ADAPTER $FASTQ_DIR/$LIB*R1*gz | head -n 100

## Check Ns:
FASTQ=$FASTQ_DIR/$LIB*R1*gz
OUTFILE=results/qc/fastq/nCounts_$LIB.R1.txt
sbatch scripts/fastq_process/countNs.sh $FASTQ $OUTFILE

FASTQ=$FASTQ_DIR/$LIB*R2*gz
OUTFILE=results/qc/fastq/nCounts_$LIB.R2.txt
sbatch scripts/fastq_process/countNs.sh $FASTQ $OUTFILE


################################################################################
##### L1_redo #####
################################################################################
LIB=newLib_failedInds

FASTQ=$FASTQ_DIR/$LIB*R1*gz
OUTFILE=results/qc/fastq/barcodeCounts_$LIB.R1.txt
sbatch scripts/fastq_process/checkBarcodes.sh $FASTQ $CUTSITE $OUTFILE

FASTQ=$FASTQ_DIR/$LIB*R2*gz
OUTFILE=results/qc/fastq/barcodeCounts_$LIB.R2.txt
sbatch scripts/fastq_process/checkBarcodes.sh $FASTQ $CUTSITE $OUTFILE

## Check adapters:
ADAPTER="GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTGAGGAATCTCGTAT"
zgrep --color=always $ADAPTER $FASTQ_DIR/$LIB*R1*gz | head -n 100

## Check Ns:
FASTQ=$FASTQ_DIR/$LIB*R1*gz
OUTFILE=results/qc/fastq/nCounts_$LIB.R1.txt
sbatch scripts/fastq_process/countNs.sh $FASTQ $OUTFILE

FASTQ=$FASTQ_DIR/$LIB*R2*gz
OUTFILE=results/qc/fastq/nCounts_$LIB.R2.txt
sbatch scripts/fastq_process/countNs.sh $FASTQ $OUTFILE
