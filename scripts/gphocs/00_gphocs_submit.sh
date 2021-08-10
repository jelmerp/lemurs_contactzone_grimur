## Prep input
FILE_ID="r03.wOutgroups.hz.mur2gri2c"
FASTA_DIR=/work/jwp37/hybridzone/seqdata/fasta_full//byLocus.final*mur2gri2c*
GPHOCS_LOCUS_DIR=analyses/gphocs/input_prep."$FILE_ID"
GPHOCS_INPUT_DIR=analyses/gphocs/input/

sbatch -p yoderlab,common,scavenger -o slurm.gphocs1createLoci."$FILE_ID".txt \
    /datacommons/yoderlab/users/jelmer/scripts/genomics/gphocs/gphocs_1_createLoci.sh "$FILE_ID" "$FASTA_DIR" "$GPHOCS_LOCUS_DIR" "$GPHOCS_INPUT_DIR"

## Run GPhocs
FILE_ID=hz.mur3gri2c
DIR_FOCAL=results/gphocs/controlfiles/reps/"$FILE_ID"/
NCORES=12

for control_file in "$DIR_FOCAL"/"$FILE_ID"*ctrl; do

    echo "$control_file"

    sbatch -n $NCORES -o slurm.gphocs_run."$(basename "$control_file")" \
        scripts/gphocs/02_gphocs_run.sh "$control_file" "$NCORES"

done
