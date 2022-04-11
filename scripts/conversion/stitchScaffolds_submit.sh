#!/bin/bash

## Settings
ref_dir=data/refgenome/
id_in=GCF_000165445.2_Mmur_3.0_genomic            # As available on NCBI
id_out=GCF_000165445.2_Mmur_3.0_genomic_stitched
nr_n=1000
maxlen=100000000
identifier=NW
scaf_exclude_infile=$ref_dir/scaffolds.nonAutosomal.txt
scaf_sizes_infile=$ref_dir/scaffolds_withLength.txt

## Get rid of superfluous info in scaffold names, and sort alphabetically:
sed 's/ Microcebus.*//' $ref_dir/$id_in.fasta | seqkit sort > $ref_dir/${id_in}2.fasta  
id_in2=${id_in}2

## Run `stitchScaffolds.sh` script
sbatch scripts/conversion/stitchScaffolds.sh \
    "$id_in2" "$id_out" "$nr_n" "$maxlen" "$identifier" "$ref_dir" \
    "$scaf_exclude_infile" "$scaf_sizes_infile"
