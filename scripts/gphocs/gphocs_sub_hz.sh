# mkdir -p analyses/gphocs/output gphocs_logfiles/

################################################################################
#### PREP INPUT ####
################################################################################
# mv /work/jwp37/hybridzone/seqdata/fasta_full/byLocus.final.r03.all.gphocs1.mac3.FS7.callableDP3.ov0.9.ovt0.8.ls100/ /work/jwp37/hybridzone/seqdata/fasta_full//byLocus.final.r03.all.gphocs1/
# mv /work/jwp37/hybridzone/seqdata/fasta_full/byLocus.final.r03.all.gphocs2.mac1.FS7.callableDP3.ov0.9.ovt0.8.ls100/ /work/jwp37/hybridzone/seqdata/fasta_full//byLocus.final.r03.all.gphocs2/
# mv /work/jwp37/hybridzone/seqdata/fasta_full/byLocus.final.r03.wOutgroups.hz.mur3gri2.mac1.FS7.callableDP3.ov0.9.ovt0.8.ls100/ /work/jwp37/hybridzone/seqdata/fasta_full//byLocus.final.hz.mur3gri2/
# mv /work/jwp37/hybridzone/seqdata/fasta_full/byLocus.final.r03.wOutgroups.hz.mur2gri2.mac1.FS7.callableDP3.ov0.9.ovt0.8.ls100/ /work/jwp37/hybridzone/seqdata/fasta_full//byLocus.final.hz.mur2gri2/

#FILE_ID=r03.all.gphocs1; FASTA_DIR=/work/jwp37/hybridzone/seqdata/fasta_full//byLocus.final.r03.all.gphocs1/
#FILE_ID=r03.all.gphocs2; FASTA_DIR=/work/jwp37/hybridzone/seqdata/fasta_full//byLocus.final.r03.all.gphocs2/
FILE_ID="r03.wOutgroups.hz.mur2gri2c"; FASTA_DIR=/work/jwp37/hybridzone/seqdata/fasta_full//byLocus.final*mur2gri2c*
#FILE_ID=r03.wOutgroups.hz.mur3gri2c; FASTA_DIR=/work/jwp37/hybridzone/seqdata/fasta_full//byLocus.final*mur3gri2c*

GPHOCS_LOCUS_DIR=analyses/gphocs/input_prep.$FILE_ID
GPHOCS_INPUT_DIR=analyses/gphocs/input/

sbatch -p yoderlab,common,scavenger -o slurm.gphocs1createLoci.$FILE_ID.txt \
	/datacommons/yoderlab/users/jelmer/scripts/gphocs/gphocs_1_createLoci.sh $FILE_ID $FASTA_DIR $GPHOCS_LOCUS_DIR $GPHOCS_INPUT_DIR

## If nr of loci given in first line of file is 0, because ls didn't work (too many files) -- see below.


################################################################################
#### RUN GPHOCS ####
################################################################################
#FILE_ID=r03.all.gphocs1
#FILE_ID=r03.all.gphocs2
#FILE_ID=hz.mur2gri2c
FILE_ID=hz.mur3gri2c

DIR_FOCAL=analyses/gphocs/controlfiles/reps/$FILE_ID/
NCORES=12
GPHOCS_COPY=1

for CFILE in $(ls $DIR_FOCAL/$FILE_ID*gac2mac2anc*190302*ctrl)
do
	#CFILE=analyses/gphocs/controlfiles/reps/hz.mur3gri2c//hz.mur3gri2c_noMig_190301_rep1.ctrl
	echo $CFILE
	
	sbatch -p yoderlab,common,scavenger -N 1 -n $NCORES -o gphocs_logfiles/slurm.gphocs_run.$(basename $CFILE).$(date +%Y%m%d-%H%M) --exclude=dcc-biostat-01,dcc-biostat-02,dcc-biostat-03 \
	/datacommons/yoderlab/users/jelmer/scripts/gphocs/gphocs_4_run.sh $CFILE $NCORES $GPHOCS_COPY
	
	printf "\n"
done



################################################################################
################################################################################
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/hybridzone/scripts/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/scripts/
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/scripts/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/scripts/

# rsync -avr /home/jelmer/Dropbox/sc_lemurs/hybridzone/analyses/gphocs/controlfiles/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/analyses/gphocs/controlfiles/

# rsync -avr jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/analyses/gphocs/output/* /home/jelmer/Dropbox/sc_lemurs/hybridzone/analyses/gphocs/output/


################################################################################
################################################################################
# find analyses/gphocs/output/*log* -maxdepth 1 -mmin +$((60)) -exec rm -f {} \;
# find gphocs_logfiles/* -maxdepth 1 -mmin +$((5)) -exec rm -f {} \;

## To change max nr of migration bands: Changed "#define MAX_MIGS 10" in src/path.h to max 20.

## If nr of loci given in first line of file is 0, because ls didn't work (too many files), then:
# SEQFILE=analyses/gphocs/input/msp3proj.eastwest2.gphocsInput.txt
# NLOCI=$(grep "fa" $SEQFILE | wc -l)
# cat $SEQFILE | sed "s/^0$/$NLOCI/" > $SEQFILE.tmp
# cp $SEQFILE $SEQFILE.backup
# mv $SEQFILE.tmp $SEQFILE

################################################################################
################################################################################
## Troubleshooting failed runs with mur3gri2:
## Check input:
INPUT=analyses/gphocs/input/r03.wOutgroups.hz.mur3gri2.gphocsInput.txt
cat $INPUT | egrep -oh "m[m|g][a-z][a-z][0-9][0-9][0-9]" | sort | uniq
#INPUT2=analyses/gphocs/input/r03.wOutgroups.hz.mur2gri2.gphocsInput.txt

INPUT=analyses/gphocs/input/r03.wOutgroups.hz.mur3gri2.gphocsInput.txt
INDS=( $(cat metadata/indSel/hz.mur3gri2.indsel.txt) )
for IND in ${INDS[@]}
do
	echo "Indiv: $IND"
	grep $IND $INPUT | wc -l
	printf "\n"
done