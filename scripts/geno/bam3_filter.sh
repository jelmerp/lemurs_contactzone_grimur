#!/bin/bash
set -e
set -o pipefail
set -u

echo -e "\n## bam3_filter.sh: Starting with script."
date
echo

################################################################################
#### PROCESS ARGS ####
################################################################################
## Command-line args:
# Help function:
Help() {
  echo "bam3_filter.sh: Filter bamfile by MMQ and/or properly paired reads."
  echo
  echo "Syntax: bam3_filter.sh [-c INT | -h | -f LOG | -m INT | -s STR] in.bam out.bam"
  echo "Options:"
  echo "c INT   Number of cores (Default: 4)"
  echo "h       Display help."
  echo "f LOG   Filter for properly paired reads TRUE/FALSE (Default: TRUE)"
  echo "m INT   Minimum mapping quality (Default: 30)"
  echo "s STR   Dir for bamfile stats (Default: geno/qc/bam/)"
  echo
}

## Defaults:
MINMAPQUAL="30"
FILTER_PROPPAIRS="TRUE"
BAMSTATS_DIR="geno/qc/bam/"
NCORES="4"

## Process optional args:
while getopts ':c:f:hq:s:' flag; do
  case "${flag}" in
    c)  NCORES="$OPTARG"
        echo "## bam3_filter.sh: Nr of cores: $NCORES"
        ;;
    f)  FILTER_PROPPAIRS="$OPTARG"
        echo "## bam3_filter.sh: Filter for properly paired reads: "$FILTER_PROPPAIRS""
        ;;
    h)  Help
        exit 0
        ;;
    q)  MINMAPQUAL="$OPTARG"
        echo "## bam3_filter.sh: Minimum mapping quality: "$MINMAPQUAL""
	    ;;
    s)  BAMSTATS_DIR="$OPTARG"
        echo "## bam3_filter.sh: Bamstats dir: $BAMSTATS_DIR"
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
echo "## bam3_filter.sh: INFILE: $INFILE"
OUTFILE=$2
echo "## bam3_filter.sh: OUTFILE: $OUTFILE"
echo

## Process args:
ID=$(basename $INFILE | sed 's/\..*//')
echo "## bam3_filter.sh: Bamfile ID: $ID"
OUTDIR=$(dirname $OUTFILE)
echo "## bam3_filter.sh: Outdir: $OUTDIR"
STATSFILE=$BAMSTATS_DIR/$ID.bamfilterstats.txt
echo "## bam3_filter.sh: ID: $ID"
echo


################################################################################
#### MORE SET-UP ####
################################################################################
NRSEQS_IN=$(samtools view -c $INFILE)   # Count nr of seqs in input bamfile

## Create dirs if needed:
[[ ! -d $OUTDIR ]] && echo "## bam3_filter.sh: Creating output dir $OUTDIR" && mkdir -p $OUTDIR
[[ ! -d $BAMSTATS_DIR ]] && echo "## bam3_filter.sh: Creating bamstats dir $BAMSTATS_DIR" && mkdir -p $BAMSTATS_DIR


################################################################################
#### RUN ####
################################################################################
## Sort & filter by minimum mapping quality:
echo "## bam3_filter.sh: Removing reads with MMQ smaller than $MINMAPQUAL..."
samtools view -bhu -q $MINMAPQUAL -@ $NCORES $INFILE > $OUTDIR/$ID.MQonly.bam
NRSEQS_POSTMQ=$(samtools view -c $OUTDIR/$ID.MQonly.bam)

## Filter for properly paired reads:
## Using "f 3" flag to filter both for paired reads and properly paired reads, see https://broadinstitute.github.io/picard/explain-flags.html 
if [ $FILTER_PROPPAIRS == TRUE ]
then
	echo -e "\n## bam3_filter.sh: Filtering for properly paired reads..." 
	samtools view -f 3 $OUTDIR/$ID.MQonly.bam -O bam > $OUTFILE
	NRSEQS_POSTPAIR=$(samtools view -c $OUTFILE)
else
	echo -e "\n## bam3_filter.sh: Not filtering for properly paired reads, renaming file..." 
	mv $OUTDIR/$ID.MQonly.bam $OUTFILE
fi


################################################################################
#### REPORTS STATS ####
################################################################################
echo -e "\n## bam3_filter.sh: Nr of sequences in raw bam file: $NRSEQS_IN"
echo -e "## bam3_filter.sh: Nr of sequences in MQ-filtered bam file: $NRSEQS_POSTMQ"
[[ $FILTER_PROPPAIRS == TRUE ]] && echo "## bam3_filter.sh: Nr of sequences in properly-paired-filtered bam file: $NRSEQS_POSTPAIR"

echo "$INFILE" > $STATSFILE
echo "Nr of sequences in raw bam file: $NRSEQS_IN" >> $STATSFILE 
echo "Nr of sequences in MQ-filtered bam file: $NRSEQS_POSTMQ" >> $STATSFILE
[[ $FILTER_PROPPAIRS == TRUE ]] && echo "Nr of sequences in properly-paired-filtered bam file: $NRSEQS_POSTPAIR" >> $STATSFILE

echo -e "\n## bam3_filter.sh: Listing output file:"
ls -lh $OUTFILE

## Remove temp files:
rm -f $OUTDIR/$ID.MQonly.bam

## Report:
echo -e "\n## bam3_filter.sh: Done with script bam3_filter.sh"
date


################################################################################
## Testing: https://www.biostars.org/p/17575/
#samtools view -F 14 -h $BAM | awk '$7 !~ /=/' | samtools view -Sb - > map2diffchroms.bam
#samtools view -F 1 -h $BAM -b > notflag1.bam
#samtools view -f 3 -h $BAM -b > flag3.bam

#samtools flagstat $BAM
# The -F 14 flag gives you reads that:
#  are mapped
#  have a mate that is mapped
#  but are not mapped in a proper pair