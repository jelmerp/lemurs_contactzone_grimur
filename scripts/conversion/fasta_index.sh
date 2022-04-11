#!/bin/bash

set -euo pipefail

## Software
BWA=software/bwa-0.7.15/bwa

## Args 
REF=$1

## Report
echo "## Indexing reference: $REF"

## Index
$BWA index "$REF"

## Report
echo "## Done with script."
date
echo


#? bwa flags:
# -a genome type; "is" for large genomes
# $BWA index -a is $RE
