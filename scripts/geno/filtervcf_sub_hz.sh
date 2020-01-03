cd /datacommons/yoderlab/users/jelmer/hybridzone

################################################################################
##### COMMON OPTIONS #####
################################################################################
VCF_DIR_MAIN=/work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/intermed
VCF_DIR_FINAL=/work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/final
QC_DIR=analyses/qc/vcf/map2mmur/map2mmur.gatk4.paired.joint/
REF=/datacommons/yoderlab/users/jelmer/seqdata/reference/mmur/GCF_000165445.2_Mmur_3.0_genomic_stitched.fasta
DP_MEAN=15
MAC=3
MEM=12
SKIP_COMMON_STEPS="-456789tew"
SKIP_FINAL_STEPS="-123"
SKIP_IN_PIP="-C" # -C: skip common steps


################################################################################
##### ALL #####
################################################################################
## r03.all:
FILE_ID=r03.all
INPUT_NAME=$FILE_ID.S3
OUTPUT_NAME=$FILE_ID
INDFILE=notany
SELECT_INDS_BY_FILE=FALSE
FILTER_INDS_BY_MISSING=TRUE
JOBNAME=filterVCF.$OUTPUT_NAME
sbatch -p yoderlab,common,scavenger -o slurm.filterVCFpip.$FILE_ID.$OUTPUT_NAME \
scripts/rad/filtering/filterVCF_FS6_pip.sh $INPUT_NAME $OUTPUT_NAME $VCF_DIR_MAIN $VCF_DIR_FINAL $QC_DIR $REF \
	$DP_MEAN $MAC $FILTER_INDS_BY_MISSING $SELECT_INDS_BY_FILE $INDFILE $MEM $JOBNAME \
	$SKIP_COMMON_STEPS $SKIP_FINAL_STEPS $SKIP_IN_PIP
	
## Without worst inds, but keep all putative hybs:
FILE_ID=r03.all
INPUT_NAME=$FILE_ID.rawSNPs.ABHet
OUTPUT_NAME=r03.keepHybs
INDFILE=metadata/indSel/r03.keepHybs.txt
SELECT_INDS_BY_FILE=TRUE
FILTER_INDS_BY_MISSING=FALSE
JOBNAME=filterVCF.$OUTPUT_NAME
SKIP_COMMON_STEPS="-12456789tew"
SKIP_IN_PIP="" # -C: skip common steps
cat metadata/r03/sampleIDsShort_r03.txt | \
	egrep -v 'mgan015|mhyb005|mgri068|mgri073|mgri074|mgri081|mgri083|mgri089|mgri091|mgri092|mgri096|mgri101|mmur054|mmur062|mmur074|mmur075' > $INDFILE
scripts/rad/filtering/filterVCF_FS6_pip.sh $INPUT_NAME $OUTPUT_NAME $VCF_DIR_MAIN $VCF_DIR_FINAL $QC_DIR $REF \
	$DP_MEAN $MAC $FILTER_INDS_BY_MISSING $SELECT_INDS_BY_FILE $INDFILE $MEM $JOBNAME \
	$SKIP_COMMON_STEPS $SKIP_FINAL_STEPS $SKIP_IN_PIP
	
## With outgroups:
FILE_ID=r03.wOutgroups
INPUT_NAME=$FILE_ID.rawSNPs.ABHet
OUTPUT_NAME=$FILE_ID
INDFILE=notany
SELECT_INDS_BY_FILE=FALSE
FILTER_INDS_BY_MISSING=TRUE
JOBNAME=filterVCF.$OUTPUT_NAME
SKIP_COMMON_STEPS="-12456789tew"
SKIP_IN_PIP="" # -C: skip common steps
scripts/rad/filtering/filterVCF_FS6_pip.sh $INPUT_NAME $OUTPUT_NAME $VCF_DIR_MAIN $VCF_DIR_FINAL $QC_DIR $REF \
	$DP_MEAN $MAC $FILTER_INDS_BY_MISSING $SELECT_INDS_BY_FILE $INDFILE $MEM $JOBNAME \
	$SKIP_COMMON_STEPS $SKIP_FINAL_STEPS $SKIP_IN_PIP
	
	
################################################################################	
##### TSIMELAHY + ALLO (EXCLUDE MANGATSIAKA) #####
################################################################################
OUTPUT_NAME=$FILE_ID.TM
INDFILE=metadata/r01/samples/selections/IDs_griseorufus_r01.txt
cat metadata/r01/samples/sampleIDs.txt | grep "mgri" > $INDFILE
SELECT_INDS_BY_FILE=TRUE
FILTER_INDS_BY_MISSING=TRUE
JOBNAME=filterVCF.$OUTPUT_NAME
sbatch -p yoderlab,common,scavenger -o slurm.filterVCFpip.$FILE_ID.$OUTPUT_NAME \
scripts/rad/filtering/filterVCF_FS6_pip.sh $INPUT_NAME $OUTPUT_NAME $VCF_DIR_MAIN $VCF_DIR_FINAL $QC_DIR $REF \
	$DP_MEAN $MAC $FILTER_INDS_BY_MISSING $SELECT_INDS_BY_FILE $INDFILE $MEM $JOBNAME \
	$SKIP_COMMON_STEPS $SKIP_FINAL_STEPS $SKIP_IN_PIP


################################################################################	
##### MANGATSIAKA + ALLO (EXCLUDE TSIMELAHY) #####
################################################################################
OUTPUT_NAME=$FILE_ID.MG
INDFILE=metadata/r01/samples/selections/IDs_griseorufus_r01.txt
cat metadata/r01/samples/sampleIDs.txt | grep "mgri" > $INDFILE
SELECT_INDS_BY_FILE=TRUE
FILTER_INDS_BY_MISSING=TRUE
JOBNAME=filterVCF.$OUTPUT_NAME
sbatch -p yoderlab,common,scavenger -o slurm.filterVCFpip.$FILE_ID.$OUTPUT_NAME \
scripts/rad/filtering/filterVCF_FS6_pip.sh $INPUT_NAME $OUTPUT_NAME $VCF_DIR_MAIN $VCF_DIR_FINAL $QC_DIR $REF \
	$DP_MEAN $MAC $FILTER_INDS_BY_MISSING $SELECT_INDS_BY_FILE $INDFILE $MEM $JOBNAME \
	$SKIP_COMMON_STEPS $SKIP_FINAL_STEPS $SKIP_IN_PIP
	
	
################################################################################
################################################################################
#### Copy files ####
# rsync -avr jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/analyses/qc/ /home/jelmer/Dropbox/sc_lemurs/hybridzone/analyses/qc/

# rsync -avr jwp37@dcc-slogin-02.oit.duke.edu:/work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/final/*FS6* /home/jelmer/Dropbox/sc_lemurs/hybridzone/seqdata/vcf/

# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/radseq/scripts/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/radseq/scripts/
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/hybridzone/scripts/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/scripts/
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/hybridzone/metadata/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/metadata/
