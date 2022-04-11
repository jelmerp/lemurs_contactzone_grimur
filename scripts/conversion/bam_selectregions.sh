#!/bin/bash
set -e
set -o pipefail
set -u

echo 
date
echo -e "## Script: $0 \n"

################################################################################
#### SET-UP ####
################################################################################
## Command-line args:
BAM_IN=$1
shift
BAM_OUT=$1
shift
REGION_FILE=$1
shift
BAMSTATS_DIR=$1
shift

## Optional args:
HEADER_TOP=notany
while getopts "h:" opt; do
	case $opt in
    	h)      HEADER_TOP=$OPTARG
                echo -e "## Using header top file $HEADER_TOP \n"
        ;;
        \?)     echo "## Error: Invalid option"
                exit 1
        ;;
	esac
done

## Process args:
ID=$(basename $BAM_IN)
STATSFILE=$BAMSTATS_DIR/$ID.bamfilterstats.txt

## Report:
echo "#########################################################################"
echo "## Variables:"
echo
echo "## $0: Bam input file:"
echo $BAM_IN
echo
echo "## $0: Bam output file:"
echo $BAM_OUT
echo
echo "## $0: File with genomic regions to select:"
echo $REGION_FILE
echo
echo "## $0: Bamstats dir:"
echo $BAMSTATS_DIR
echo
echo "## $0: Bamstats file:"
echo $STATSFILE
echo
echo "## $0: Header top:"
echo $HEADER_TOP
echo

[[ ! -d $BAMSTATS_DIR ]] && mkdir -p $BAMSTATS_DIR  # Make statsdir if needed


################################################################################
#### SELECT REGIONS ####
################################################################################
echo "##########################################################################"
echo "## $0: Selecting only specified regions..."
samtools view -b -L $REGION_FILE $BAM_IN > $BAM_OUT


################################################################################
#### REHEADER ####
################################################################################
if [ HEADER_TOP != "notany" ]; then
    echo "######################################################################"
    echo -e "\n## $0: Reheading bamfile..."

    ## Temporary files:
    FIXED_HEADER=tmp.$ID.header
    BAM_TMP=tmp.$ID.bam
    mv $BAM_OUT $BAM_TMP # Move subsetted bam to temp file to allow for reheadering

    ## Reheader:
    cp $HEADER_TOP $FIXED_HEADER # Get universal bam header top with all scaffold
    samtools view -H $BAM_TMP | egrep -v "@HD|@SQ" >> $FIXED_HEADER # Add sample-specific rest of bam header
    cat $FIXED_HEADER <(samtools view $BAM_TMP) | samtools view -bo $BAM_OUT - # Replace the header in the bamfile
    
    ## Check nr of chromosomes:
    NCHROM=$(samtools view -H $BAM_OUT | grep "@SQ" | wc -l)
    echo -e "\n## $0: Nr of chroms/scaffolds in new bam header: $NCHROM"

    ## Remove temp files:
    rm -f $FIXED_HEADER
    rm -f $BAM_TMP
fi


################################################################################
#### REPORT ####
################################################################################
echo "#########################################################################"
NRSEQS_IN=$(samtools view -c $BAM_IN)
NRSEQS_OUT=$(samtools view -c $BAM_OUT)
NRSEQS_REMOVED=$(($NRSEQS_IN - $NRSEQS_OUT))

echo -e "\n## $0: Nr of sequences in: $NRSEQS_IN"
echo -e "## $0: Nr of sequences removed: $NRSEQS_REMOVED"
echo -e "## $0: Nr of sequences out: $NRSEQS_OUT"
echo -e "## Nr of sequences before extracting specified scaffolds: $NRSEQS_IN" >> $STATSFILE
echo -e "## Nr of sequences after extracting specified scaffolds: $NRSEQS_OUT" >> $STATSFILE

echo -e "\n### $0: Listing output file:"
ls -lh $BAM_OUT

echo -e "\n## $0: Done with script."
date
