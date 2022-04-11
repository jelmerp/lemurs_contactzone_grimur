#!/bin/bash

set -euo pipefail

# SET-UP -----------------------------------------------------------------------
echo "## Starting script: stitchScaffolds.sh"
date
echo

## Software and scripts
PYTHON3=software/Python-3.6.3/python
JAVA=software/java_1.8.0/jre1.8.0_144/bin/java
PICARD=software/picard_2.13.2/picard.jar
SAMTOOLS=software/samtools-1.6/samtools
SCAFFOLD_STITCHER=software/scaffoldStitcher/ScaffoldStitcher.py # https://bitbucket.org/dholab/scaffoldstitcher/src

SCRIPT_LOOKUP scripts/pipeline_misc/stitchScaffolds_extract.R

## Command-line args
ID_IN=$1
shift
ID_OUT=$1
shift
NR_N=$1
shift
MAXLENGTH=$1
shift
IDENTIFIER=$1
shift
REF_DIR=$1
shift
SCAF_EXCLUDE_INFILE=$1
shift
SCAF_SIZES_INFILE=$1
shift

## Process args:
FASTA_IN=$REF_DIR/$ID_IN.fasta
FASTA_OUT=$REF_DIR/$ID_OUT.fasta

SCAF_INDEX_INFILE=$REF_DIR/$ID_OUT.scaffoldIndex.txt # Created by scaffoldStitcher
SCAF_INDEX_OUTFILE=$REF_DIR/$ID_OUT.scaffoldIndexLookup.txt # Created by stitchScaffolds_extract.R
SCAF_EXCLUDE_OUTFILE=$REF_DIR/$ID_OUT.nonAutosomalCoords.bed # Created by stitchScaffolds_extract.R
SCAFLIST_FILE=$REF_DIR/$ID_OUT.scaffoldList.txt # Created by stitchScaffolds_extract.R

## Report:
echo "## Ref dir: $REF_DIR"
echo "## ID in: $ID_IN"
echo "## ID out: $ID_OUT"
echo "## Nr Ns between scaffolds: $NR_N"
echo "## Max length of superscaffold: $MAXLENGTH"
echo "## Identifier to distinguish scaffolds that should be merged: $IDENTIFIER"
echo
echo "## Fasta in: $FASTA_IN"
echo "## Fasta out: $FASTA_OUT"
echo
echo "## Infile with scaffold sizes: $SCAF_SIZES_INFILE"
echo "## Infile with index of superscaffolds-to-scaffolds: $SCAF_INDEX_INFILE"
echo "## Infile with scaffolds to exclude: $SCAF_EXCLUDE_INFILE"
echo "## Outfile with lookup for superscaffolds-to-scaffolds $SCAF_INDEX_OUTFILE"
echo "## Outfile (bed) with regions to exclude from bam: $SCAF_EXCLUDE_OUTFILE"
echo "## Outfile with list of scaffolds: $SCAFLIST_FILE"
echo -e "-----------------\n"


# STITCH SCAFFOLDS -------------------------------------------------------------
echo "## Creating stitched ref genome..."
"$PYTHON3" "$SCAFFOLD_STITCHER" \
    -fasta "$FASTA_IN" \
    -identifier "$IDENTIFIER" \
    -nlength "$NR_N" \
    -maxlength "$MAXLENGTH" > "$FASTA_OUT"

#nlength = N spacer length between scaffolds; maxlength = max. super scaffold length

mv "$REF_DIR"/"$ID_IN"_scaffold_index.txt "$SCAF_INDEX_INFILE"


# INDEX NEW FASTA --------------------------------------------------------------
## Index new fasta with samtools, picard, and bwa:
echo -e "\n## Indexing ref genome with samtools..."
$SAMTOOLS faidx "$FASTA_OUT"

echo -e "\n## Indexing ref genome with picard..."
[[ -f $REF_DIR/$ID_OUT.dict ]] && rm "$REF_DIR"/"$ID_OUT".dict
$JAVA -Xmx4g -jar $PICARD CreateSequenceDictionary \
    R="$FASTA_OUT" O="$REF_DIR"/"$ID_OUT".dict

echo -e "\n## Indexing ref genome with bwa..."
scripts/misc/indexGenome.sh "$FASTA_OUT"


# CREATE BEDFILE AND LOOKUP TABLE ----------------------------------------------
## Create bedfile with regions to exclude, and superscaffold-to-scaffold location lookup table:
echo "## Running R script for superscaffold-to-scaffold location lookup table..."
Rscript "$SCRIPT_LOOKUP" "$SCAF_SIZES_INFILE" "$SCAF_INDEX_INFILE" \
	"$SCAF_EXCLUDE_INFILE" "$SCAF_INDEX_OUTFILE" "$SCAF_EXCLUDE_OUTFILE" \
    "$SCAFLIST_FILE" "$NR_N"


## Report:
echo -e "\nDone with script stitchScaffolds.sh"
date
echo