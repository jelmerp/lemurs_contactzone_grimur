#!/bin/bash

set -euo pipefail

# SET-UP -----------------------------------------------------------------------
## Scripts
SCRIPT_VCF2FASTA=scripts/conversion/vcf2fasta.sh
SCRIPT_FASTA2NEXUS=scripts/conversion/convert2nexus.sh
SPIDFILE=scripts/conversion/fasta2nexus.spid
SCRIPT_SPLITSTREE=scripts/splitstree/splitstree.sh

## Command-line args:
file_id=$1
vcf_dir=$2
fasta_dir=$3
nexus_dir=$4
outdir=$5
mem=$6 # At least 20GB for FASTA2NEXUS

## Other variables:
fasta="$fasta_dir"/"$file_id".fasta
nexus="$nexus_dir"/"$file_id".nexus
outfile="$outdir"/"$file_id".nexus

mkdir -p "$fasta_dir" "$nexus_dir" "$outdir"

## Report:
echo "## Starting script 01_splitstree_runner.sh"
date
echo
echo "## File ID:           $file_id"
echo "## Vcf dir:           $vcf_dir"
echo "## Fasta dir:         $fasta_dir"
echo "## Nexus dir:         $nexus_dir"
echo "## Output dir:        $outdir"
echo "## Memory:            $mem"
echo
echo "## Fasta file:        $fasta"
echo "## Nexus file:        $nexus"
echo "## Output file:       $outfile"
echo

# CONVERSIONS ------------------------------------------------------------------
## VCF2FASTA
echo -e "-----------------------------\n"
echo "## Running VCF2FASTA script..."
SCAFFOLD=ALL
"$SCRIPT_VCF2FASTA" "$file_id" "$vcf_dir" "$fasta_dir" "$SCAFFOLD"

## FASTA2NEXUS
echo -e "-----------------------------\n"
echo "## Running FASTA2NEXUS script..."
INFORMAT=FASTA
"$SCRIPT_FASTA2NEXUS" "$fasta" "$nexus" "$INFORMAT" "$SPIDFILE" "$mem"


# RUN SPLITSTREE ---------------------------------------------------------------
echo -e "-----------------------------\n"
echo "## Running splitstree..."
"$SCRIPT_SPLITSTREE" "$nexus" "$outfile"
