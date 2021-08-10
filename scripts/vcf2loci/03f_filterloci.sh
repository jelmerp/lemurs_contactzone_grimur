#!/bin/bash
set -euo pipefail

## Software:
# R

## Scripts:
script_filterloci=03e_filterloci.R

## Command-line args:
locusstats_intermed=$1
file_ld=$2
dir_locusfasta_intermed=$3
dir_locusfasta_final=$4

## Submit script:
Rscript "$script_filterloci" \
    "$locusstats_intermed" "$file_ld" "$dir_locusfasta_intermed" "$dir_locusfasta_final"
