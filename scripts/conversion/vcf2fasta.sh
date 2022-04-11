#!/bin/bash

## Bash strict settings
set -euo pipefail

# SET-UP -----------------------------------------------------------------------
## Software
VCFTAB2FASTA=scripts/vcf_tab_to_fasta_alignment.pl
# export PERL5LIB=/dscrhome/rcw27/programs/vcftools/vcftools-master/src/perl/

## Get Positional args
file_id=$1
shift
indir=$1
shift
outdir=$1
shift

## Optional args
scaffold="ALL"

while getopts 's:' flag; do
    case "${flag}" in
    s) scaffold="${OPTARG}" ;;
    esac
done

## Report
echo -e "\n## Starting script vcf2fasta.sh"
date
echo

## Get VCF
if [ -e "$indir"/"$file_id".vcf.gz ]; then
    echo "## vcf2fasta.sh: Zipped VCF detected..."
    zipped_vcf=true
    infile="$indir"/"$file_id".vcf.gz
elif [ -e "$indir"/"$file_id".vcf ]; then
    echo "## vcf2fasta.sh: Unzipped VCF detected..."
    zipped_vcf=false
    infile="$indir"/"$file_id".vcf
else
    echo -e "\n## ERROR: No VCF file $indir/$file_id.vcf(.gz)" && exit 1
fi

## Make outdir if it doesn't exist
mkdir -p "$outdir"

## Report
echo "## File ID: $file_id"
echo "## Input dir: $indir"
echo "## Output dir: $outdir"
echo "## Input file: $infile"
echo "## Is input file zipped (true/false): $zipped_vcf"
echo "## Scaffold ('ALL' if entire vcf will be processed): $scaffold"


# CONVERT VCF TO FASTA ---------------------------------------------------------
if [ "$scaffold" != "ALL" ]; then

    echo -e "\n## vcf2fasta.sh: Extracting single scaffold $scaffold from vcf file..."

    file_id="$file_id"_"$scaffold"

    [[ $zipped_vcf == TRUE ]] && vcftools_option="--gzvcf"
    [[ $zipped_vcf == FALSE ]] && vcftools_option="--vcf"

    vcftools $vcftools_option "$infile" --chr "$scaffold" \
        --remove-filtered-all --recode --recode-INFO-all --stdout |
        vcf-to-tab |
        sed 's/\./N/g' >"$outdir"/"$file_id".tab

else

    echo -e "\n## Processing ALL scaffolds / entire vcf file..."
    echo -e "\n## Converting vcf to tab-delimited file using vcftools perl utility vcf-to-tab..."
    [[ $zipped_vcf == TRUE ]] && zcat "$infile" | vcf-to-tab | sed 's/\./N/g' >"$outdir"/"$file_id".tab
    [[ $zipped_vcf == FALSE ]] && cat "$infile" | vcf-to-tab | sed 's/\./N/g' >"$outdir"/"$file_id".tab

fi

echo -e "\n## Creating varpos file w/ indices of var positions - $outdir/$file_id.varpos"
cut -f 1,2 "$outdir"/"$file_id".tab | tail -n +2 >"$outdir"/"$file_id".varpos

echo -e "\n## Converting tab to fasta...: $outdir/$file_id.fasta"
perl $VCFTAB2FASTA -i "$outdir"/"$file_id".tab >"$outdir"/"$file_id".fasta

echo -e "\n## Outputting varscaffold...: $outdir/$file_id.varscaffold"
sed ':a;N;$!ba;s/\n/\t/g' "$outdir"/"$file_id".fasta |
    sed 's/\t>/\n/g' |
    sed 's/\t/ /' |
    sed 's/\t//g' |
    sed 's/>//g' >"$outdir"/"$file_id".varscaffold


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing output file:"
ls -lh "$outdir"/"$file_id".fasta
echo -e "\n## Done with script."
date
