#!/bin/bash

# SETUP ------------------------------------------------------------------------
## Bash strict settings
set -euo pipefail

## Software & scripts
PGDS=software/PGDSpider_2.1.1.3/PGDSpider2-cli.jar

## Positional args:
infile=$1
outfile=$2
spidfile=$3
mem=$4

## Report:
echo
date
echo "## Starting script vcf2newhybrids.sh"
echo "## Input file:            $infile"
echo "## Output file:           $outfile"
echo "## PGDS-Spider file:      $spidfile"
echo "## Memory:                $mem"
echo -e "---------------------\n"


# CONVERT ----------------------------------------------------------------------
## Run PGDS-Spider:
echo "## Converting with PGDS-Spider..."
java -Xmx"$mem"G -Xms"$mem"G -jar "$PGDS" \
    -inputfile "$infile" \
    -inputformat VCF \
    -outputfile "$outfile".tmp \
    -outputformat NEWHYBRIDS \
    -spid "$spidfile"

## Get sample IDs into first column instead of numbers:
echo "## Editing sample names..."
head -n 5 "$outfile".tmp > "$outfile".tmp.firstlines 
tail -n +6 "$outfile".tmp > "$outfile".tmp.cut
bcftools query -l "$infile" | paste - <( cut -d ' ' -f 2- "$outfile".tmp.cut ) > "$outfile".tmp.replaced
cat "$outfile".tmp.firstlines "$outfile".tmp.replaced > "$outfile"


# WRAP UP ----------------------------------------------------------------------
## Remove temporary files
echo -e "\n## Removing temporary files..."
rm -f "$outfile".tmp*

## Report
echo -e "\n## Listing Newhybrids input file: $outfile"
ls -lh "$outfile"

echo -e "\n## Showing first three columns of NewHybrids input file:"
cat "$outfile" | cut -f 1,2,3 -d " "

echo -e "\n ## Done with script vcf2newhybrids.sh"
date
echo
