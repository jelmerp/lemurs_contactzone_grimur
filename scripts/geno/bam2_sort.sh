#!/bin/bash
set -e
set -o pipefail
set -u

echo -e "\n## bam2_sort.sh: Starting with script."
date
echo

################################################################################
#### SET-UP ####
################################################################################
# Help function:
Help() {
  echo "bam2_sort.sh: Sort bamfile by coordinate."
  echo
  echo "Syntax: bam2_filter.sh [-c INT | -h ] in.bam out.bam"
  echo "Options:"
  echo "c INT   Number of cores (Default: 1)"
  echo "h       Display help."
  echo "m INT	Memory in GB (Default: 4)"
  echo
}

## Defaults:
NCORES="1"
MEM="4"

## Optional args:
while getopts ':m:c:h' flag; do
  case "${flag}" in
    c)  NCORES="$OPTARG"
        echo "## bam2_sort.sh: Nr of cores: $NCORES"
        ;;
    h)  Help
        exit 0
        ;;
    m)  MEM="$OPTARG"
        echo "## bam2_sort.sh: Memory: "$MEM""
	    ;;
	\?) echo "Error: Invalid option" # https://opensource.com/article/19/12/help-bash-program
        exit 1
		;;
	:)  echo "Option -$OPTARG requires an argument." >&2
      	exit 1
      	;;
  esac
done
shift $(( OPTIND - 1 ))

## Positional args:
INFILE=$1
echo "## bam2_sort.sh: INFILE: $INFILE"
OUTFILE=$2
echo "## bam2_sort.sh: OUTFILE: $OUTFILE"
echo

## Create output dir if needed:
OUTDIR=$(dirname $OUTFILE)
[[ ! -d $OUTDIR ]] && echo "## bam2_sort.sh: Creating output dir $OUTDIR" && mkdir -p $OUTDIR
[[ ! -d tmpdir ]] && mkdir tmpdir

################################################################################
#### RUN ####
################################################################################
echo "## bam2_sort.sh: Sorting bam file:"
samtools sort -@ $NCORES -m ${MEM}G -T tmpdir $INFILE > $OUTFILE

## Report:
echo -e "\n## bam2_sort.sh: Listing output file:"
ls -lh $OUTFILE

echo -e "\n## bam2_sort.sh: Done with script bam2_sort.sh"
date