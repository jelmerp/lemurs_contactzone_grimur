## Settings
VCF_DIR=data/vcf/gatk/final/           # Dir with VCF files
PREP_INPUT=TRUE                        # Prep Treemix input files TRUE/FALSE
MINMIG=0                               # Min. number of migration events to run Treemix for
MAXMIG=10                              # Max. number of migration events to run Treemix for
TREEMIX_DIR=results/treemix/           # Treemix dir
POPFILE=metadata/hzlookup_bysample.txt # Metadata file with sample-to-pop mappings
POPCOLUMN=pop2                         # Column in $POPFILE with population info
ROOT=ruf                               # Outgroup population

FILE_IDS=(hzproj1.mac1.FS6 hzproj1.mac3.FS6) # VCF file IDs

## Run Treemix
for FILE_ID in "${FILE_IDS[@]}"; do
    echo "## File ID: $FILE_ID"
    scripts/treemix/01_treemix_runner.sh "$FILE_ID" "$VCF_DIR" "$PREP_INPUT" "$MINMIG" "$MAXMIG" "$ROOT" "$TREEMIX_DIR"" $POPFILE" "$POPCOLUMN"
done
