#!/bin/bash

## Software/scripts
NEWHYBRIDS=software/newhybrids/newhybrids-no-gui-linux.exe

## Bash strict settings
set -euo pipefail

## Positional args:
infile=$1
outdir=$2
burnin=$3
nsweeps=$4

[[ -z $burnin ]] && burnin=10000      # 10k=default
[[ -z $nsweeps ]] && nsweeps=50000    # 50k=default

## Report:
date
echo "Starting script newhybrids.sh"
echo "Input file:                      $infile"
echo "Output dir:                      $outdir"
echo "Number of burn-in sweeps:        $burnin"
echo "Number of sweeps after burn-in:  $nsweeps"
echo

## Process args:
mkdir -p "$outdir"
cd "$outdir" || exit 1


# RUN NEWHYBRIDS ---------------------------------------------------------------
echo "## Starting newhybrids run..."
$NEWHYBRIDS -d "$infile" --burn-in "$burnin" --num-sweeps "$nsweeps" --no-gui 


# WRAP UP ----------------------------------------------------------------------
echo "## Listing files in output dir $outdir:"
ls -lh "$outdir"
echo -e "\n## Done with script newhybrids.sh"
date
echo
