#!/bin/bash

set -euo pipefail

# SET-UP --------------------------------------------------------------------
## Software and scripts
SCRIPT_PREPINPUT=scripts/admixtools/admixtools_prepInput.sh
SCRIPT_ATOOLS=scripts/admixtools/admixtools.sh

## Positional args:
file_id=$1
shift
run_id=$1
shift
vcf_dir=$1
shift
plink_dir=$1
shift
atools_dir=$1
shift
indfile=$1
shift
popfile=$1
shift
vcf2plink=$1
shift
create_indfile=$1
shift
subset_indfile=$1
shift
atools_mode=$1
shift
inds_metadata=$1
shift
id_column=$1
shift
groupby=$1

## Process parameters
file_id_full=${file_id}$run_id

indir=$atools_dir/input/
outdir=$atools_dir/output/

[[ $atools_mode == "D" ]] && parfile=$indir/parfile_dmode_$file_id_full.txt
[[ $atools_mode == "F4RATIO" ]] && parfile=$indir/parfile_f4ratio_$file_id_full.txt

mkdir -p "$atools_dir"/input "$atools_dir"/output/raw "$plink_dir"

## Report:
echo
echo "## Starting script."
date
echo
echo "## File ID:                      $file_id"
echo "## Run ID:                       $run_id"
echo "## VCF dir:                      $vcf_dir"
echo "## PLINK dir:                    $plink_dir"
echo "## Create indfile (TRUE/FALSE):  $create_indfile"
echo "## Subset indfile (TRUE/FALSE):  $subset_indfile"
echo
echo "## Admixtools mode:              $atools_mode"
echo
echo "## Metadata file:                $inds_metadata"
echo "## ID column:                    $id_column"
echo "## Group-by column:              $groupby"
echo
echo "## Indfile (output):             $indfile"
echo "## Popfile (input):              $popfile"
echo "## Parfile:                      $parfile"
echo -e "--------------------\n"


# PREP INPUT -------------------------------------------------------------------
echo "## Calling script to prep input files..."
$SCRIPT_PREPINPUT "$file_id" "$vcf_dir" "$plink_dir" "$vcf2plink" \
    "$create_indfile" "$subset_indfile" "$indfile" "$popfile" "$parfile" \
    "$atools_mode" "$inds_metadata" "$id_column" "$groupby"


# RUN ADMIXTOOLS ---------------------------------------------------------------
## D-mode
if [ "$atools_mode" == "D" ]; then
	
    output=$outdir/$file_id_full.dmode.out
	
	echo -e "\n## Running admixtools in D mode..."
	echo "## Output: $output"
	
	for popfile_line in $(seq 1 "$(wc -l < "$popfile")"); do
		
        echo -e "\n#### Line nr: $popfile_line"
		head -n "$popfile_line" "$popfile" | tail -n 1
		
		$SCRIPT_ATOOLS "$file_id_full" "$popfile_line" "$parfile" "$output" "$atools_mode"
	done
	
	## Combine output into single file:
	grep -h "result" "$output".line* > "$output"
fi

## F4-ratio-mode
if [ "$atools_mode" == "F4RATIO" ]; then
	
    output=$outdir/$file_id_full.f4ratio.out
	popfile_line=ALL
	
    echo "## Running admixtools in f4ratio mode..."
	echo "## Output: $output"
	
	"$SCRIPT_ATOOLS" "$file_id_full" "$popfile_line" "$parfile" "$output" "$atools_mode"
fi


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Admixtools output: $output"
cat "$output"
echo -e "\n## Done with script."
date
echo
