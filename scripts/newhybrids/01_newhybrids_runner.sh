#!/bin/bash

# SETUP ------------------------------------------------------------------------
## Software & scripts:
SCRIPT_CONVERT=scripts/newhybrids/vcf2newhybrids.sh
SCRIPT_NEWHYBRIDS=scripts/newhybrids/newhybrids.sh

## PGDSpider file for conversion
SPIDFILE=scripts/newhybrids/vcf2newhybrids.spid

## Bash strict settings
set -euo pipefail

## Positional args
run_id=$1
vcf=$2
newhyb_in=$3
newhyb_outdir=$4
mem=$4
burnin=$5
nsweeps=$6

## Process parameters
newhyb_outdir_full="$newhyb_outdir"/"$run_id"/

## Report
date
echo "## Starting script newhybrids_runner.sh"
echo "## vcf file:                         $vcf"
echo "## Newhybrids input file:            $newhyb_in"
echo "## Newhybrids output dir:            $newhyb_outdir_full"
echo "## Memory:                           $mem"
echo "## Number of burn-in sweeps:         $burnin"
echo "## Number of sweeps after burn-in:   $nsweeps"
echo "## SPID-file:                        $SPIDFILE"
echo


# CONVERT vcf TO NEWHYBRIDS INPUT ----------------------------------------------
echo -e "\n## Calling conversion script..."
$SCRIPT_CONVERT "$vcf" "$newhyb_in" "$SPIDFILE" "$mem"


# RUN NEWHYBRIDS ---------------------------------------------------------------
echo -e "\n## Calling newhybrids script..."
$SCRIPT_NEWHYBRIDS "$newhyb_in" "$OUTDIR_NEWHYB" "$burnin" "$nsweeps"
