#!/bin/bash
set -euo pipefail

## Constants:
CALLABLE_COMMAND="--minDepth 3"

## software:
JAVA=/path/to/java
GATK=/path/to/gatk3.8.0
# bedtools
# samtools
# bgzip

## Command-line args:
ind=$1
set_id=$2
bam=$3
vcf_altref=$4
bed_removed_sites=$5
dir_indfasta=$6
dir_bed=$7
ref=$8

## Process args:
[[ ! -d $dir_indfasta ]] && mkdir -p "$dir_indfasta"
[[ ! -d $dir_bed ]] && mkdir -p "$dir_bed"

callable_summary="$dir_bed"/"$ind".callableloci.sumtable.txt

bed_out="$dir_bed"/"$ind".callablelocioutput.bed
bed_notcallable="$dir_bed"/"$ind".noncallable.bed
bed_callable="$dir_bed"/"$ind".callable.bed

fasta_altref="$dir_indfasta"/"$ind".altref.fasta
fasta_masked="$dir_indfasta"/"$ind".altrefmasked.fasta

## Report:
date
echo "## Starting script $0"
echo "## Sample ID:                               $ind"
echo "## VCF ID:                                  $set_id"
echo "## BAM file:                                $bam"
echo "## VCF file - for producing altref:         $vcf_altref"
echo "## BEDfile with sites removed by filtering: $bed_removed_sites"
echo "## FASTA dir:                               $dir_indfasta"
echo "## BED dir:                                 $dir_bed"
echo "## Reference genome:                        $ref"

[[ ! -e $bam.bai ]] && echo "## Indexing bam..." && samtools index "$bam"

# -----------------------------------------------------------------------------
# Step 1 -- run gatk callable-loci
## Using ref genome and bamfiles, produce bedfile for sites that are (non-)callable for a single sample
echo "## Running GATK CallableLoci..."

## Run callableloci:
$JAVA -Xmx4G -jar $GATK -T CallableLoci \
    -R "$ref" \
    -I "$bam" \
    -summary "$callable_summary" \
    "$CALLABLE_COMMAND" -o "$bed_out"

## Edit bedfile to include only non-callable loci:
grep -v "callable" "$bed_out" >"$bed_notcallable"
grep "callable" "$bed_out" >"$bed_callable"

# -----------------------------------------------------------------------------
# Step 2 -- run gatk fasta-alternate-reference-maker
## Using ref genome and vcf file, produce whole-genome fasta file for a single sample:
echo -e "## running gatk FastaAlternateReferenceMaker...\n"

## Run GATK:
$JAVA -Xmx4G -jar $GATK -T FastaAlternateReferenceMaker \
    -IUPAC "$ind" -R "$ref" -V "$vcf_altref" -o "$fasta_altref"

## Edit FASTA headers:
sed -i -e 's/:1//g' -e 's/>[0-9]* />/g' "$fasta_altref"

## Count bases:
n_acgt=$(grep -Eo "A|C|G|T" "$fasta_altref" | wc -l)
n_ambig=$(grep -Eo "M|R|W|S|Y|K" "$fasta_altref" | wc -l)
echo "## number of A/C/G/Ts in fasta_altref: $n_acgt"
echo "## number of het sites (as counted by ambig codes) in fasta_altref: $n_ambig"

# -----------------------------------------------------------------------------
# Step 3 -- run bedtools maskfasta ---------------------------------------------
## In whole-genome fasta for a given sample, mask:
## a) sites identified as non-callable by callableloci, and
## b) sites removed during vcf-filtering.

echo -e "## running bedtools maskfasta...\n"
## Masking non-callable sites:
bedtools maskfasta -fi "$fasta_altref" -bed "$bed_notcallable" -fo "$fasta_masked".intermed.fasta
## Masking removed (filtered-out) sites:
bedtools maskfasta -fi "$fasta_masked".intermed.fasta -bed "$bed_removed_sites" -fo "$fasta_masked"
## Counting Ns in the FASTA files:
ncount_fasta_altref=$(grep -Fo "N" "$fasta_altref" | wc -l)
ncount_fasta_masked_intermed=$(grep -Fo "N" "$fasta_masked".intermed.fasta | wc -l)
ncount_fasta_masked=$(grep -Fo "N" "$fasta_masked" | wc -l)
## Report:
echo "## Nr Ns in fasta_altref:           $ncount_fasta_altref"
echo "## Nr Ns in fasta_masked_intermed:  $ncount_fasta_masked_intermed"
echo "## Nr Ns in fasta_masked:           $ncount_fasta_masked"

# ------------------------------------------------------------------------------
## Remove intermediate files:
rm -f "$fasta_masked".intermed.fasta

## Report:
echo "## done with script."
date
