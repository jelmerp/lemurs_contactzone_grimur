#!/bin/bash

## Prep input file
file_id="hz.mur3gri2c"
fasta_dir=results/vcf2loci/
gphocs_locus_dir=results/gphocs/input_prep."$file_id"
gphocs_input_dir=results/gphocs/input/

sbatch scripts/gphocs/01_gphocs_prepinput.sh \
    "$file_id" "$fasta_dir" "$gphocs_locus_dir" "$gphocs_input_dir"

## Run GPhocs
file_id=hz.mur3gri2c
dir_focal=results/gphocs/controlfiles/reps/"$file_id"/
ncores=12

for control_file in "$dir_focal"/"$file_id"*ctrl; do
    sbatch scripts/gphocs/02_gphocs_run.sh "$control_file" "$ncores"
done
