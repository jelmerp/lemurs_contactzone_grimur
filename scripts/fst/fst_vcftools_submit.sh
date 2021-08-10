#!/bin/bash

# GENERAL SETTINGS -------------------------------------------------------------

indir=results/geno/vcf/gatk
lookup=metadata/hzlookup_bysample.txt
winsize=50000
stepsize=50000

# FOR CONTACT ZONE INDS ONLY ---------------------------------------------------

## Settings
file_id=r03.all.mac1.FS6

outdir=results/fst/output/"$file_id"
popfile_dir=results/fst/input/"$file_id"
popcombs_file=$popfile_dir/popcombs.txt # file with on each line a pair of pops to compute Fst for

## Create popfiles
mkdir -p "$popfile_dir"
awk '$4 == "mgri"' $lookup | cut -f 1 >"$popfile_dir"/gri-C.txt
awk '$4 == "mmur"' $lookup | cut -f 1 >"$popfile_dir"/mur-C.txt
awk '$4 == "mgri"' $lookup | grep "sympatric" | cut -f 1 >"$popfile_dir"/gri_sym.txt
awk '$4 == "mgri"' $lookup | grep "parapatric" | cut -f 1 >"$popfile_dir"/gri_para.txt
awk '$4 == "mmur"' $lookup | grep "sympatric" | cut -f 1 >"$popfile_dir"/mur_sym.txt
awk '$4 == "mmur"' $lookup | grep "parapatric" | cut -f 1 >"$popfile_dir"/mur_para.txt

## Run VCFtools for each pop combination
while read -r pop1 pop2; do
    echo "Pop 1: $pop1 / Pop 2: $pop2"

    popfile1="$popfile_dir"/$pop1.txt
    popfile2="$popfile_dir"/$pop2.txt

    scripts/fst/fst_vcftools.sh "$indir" "$file_id" "$popfile1" "$popfile2" \
        "$winsize" "$stepsize" "$outdir" "$popfile_dir"

done <"$popcombs_file"

# FOR ALL INDS -----------------------------------------------------------------

## Settings
file_id=hzproj1.mac1.FS6

lookup_allo=metadata/hzlookup_allopatric.txt
outdir=results/fst/output/"$file_id"
popfile_dir=results/fst/input/"$file_id"
popcombs_file=$popfile_dir/popcombs.txt # file with on each line a pair of pops to compute Fst for

## Create popfiles
mkdir -p "$popfile_dir"
grep "gri-W" $lookup_allo | cut -f 1 >"$popfile_dir"/gri-W.txt
grep "mur-E" $lookup_allo | cut -f 1 >"$popfile_dir"/mur-E.txt
grep "mur-W" $lookup_allo | cut -f 1 >"$popfile_dir"/mur-W.txt
awk '$4 == "mgri"' $lookup | cut -f 1 >"$popfile_dir"/gri-C.txt
awk '$4 == "mmur"' $lookup | cut -f 1 >"$popfile_dir"/mur-C.txt

## Run VCFtools for each pop combination
while read -r pop1 pop2; do
    echo "Pop 1: $pop1 / Pop 2: $pop2"

    popfile1="$popfile_dir"/$pop1.txt
    popfile2="$popfile_dir"/$pop2.txt

    scripts/fst/fst_vcftools.sh "$indir" "$file_id" "$popfile1" "$popfile2" \
        "$winsize" "$stepsize" "$outdir" "$popfile_dir"

done <"$popcombs_file"
