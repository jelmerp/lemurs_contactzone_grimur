#!/bin/bash

## Bash script settings
set -euo pipefail

## Other scripts/software
PGDS=software/PGDSpider2-cli.jar # https://anaconda.org/bioconda/pgdspider

## Command-line args
infile=$1
outfile=$2
informat=$3
spidfile=$4
mem=$5

## Report
echo "## Starting script convert2nexus.sh"
date
echo "## Input file: $infile"
echo "## Output file: $outfile"
echo "## Input format: $informat"
echo "## SPID file: $spidfile"
echo "## Memory: $mem"

## Run PGDS Spider
java -Xmx${mem}G -Xms${mem}G -jar $PGDS \
    -inputfile "$infile" \
    -inputformat "$informat" \
    -outputfile "$outfile".tmp \
    -outputformat NEXUS -spid "$spidfile"

## Remove block with "Taxon Sets" from Nexus file
head -n -1 "$outfile".tmp |
    grep -v "BEGIN SETS" | grep -v "TaxSet" | grep -v "TaxPartition" |
    grep -v pop_1 >"$outfile"

rm "$outfile".tmp *log

## Report
echo -e "\n## Listing output file:"
ls -lh "$outfile"
echo -e "\n## Done with script convert2nexus.sh"
date
