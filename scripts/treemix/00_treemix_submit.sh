## Settings
FILE_IDS=(hzproj1.mac1.FS6)            # VCF file IDs
VCF_DIR=data/geno/vcf/                 # Dir with VCF files
TREEMIX_DIR=results/treemix/           # Treemix dir
POPFILE=metadata/hzlookup_bysample.txt # Metadata file with sample-to-pop mappings
POPCOLUMN=pop                          # Column in $POPFILE with population info
PREP_INPUT=TRUE                        # Prep Treemix input files TRUE/FALSE
MINMIG=0                               # Min. number of migration events to run Treemix for
MAXMIG=10                              # Max. number of migration events to run Treemix for
ROOT=ruf                               # Outgroup population

## Run Treemix
for FILE_ID in "${FILE_IDS[@]}"; do
    echo "## File ID: $FILE_ID"
    scripts/treemix/01_treemix_runner.sh \
        "$FILE_ID" "$VCF_DIR" "$PREP_INPUT" "$MINMIG" "$MAXMIG" "$ROOT" \
        "$TREEMIX_DIR"" $POPFILE" "$POPCOLUMN"
done
