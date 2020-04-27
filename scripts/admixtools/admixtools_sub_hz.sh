#mkdir -p analyses/admixtools/input/ analyses/admixtools/output/raw /work/jwp37/hybridzone/seqdata/plink

## TO DO: PUT THIS IN FILTERING PIPELINE
# VCF=/work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/final/backup/r03.wOutgroups.mac3.FS6.vcf.gz
# bcftools query -l $VCF > analyses/qc/vcf/map2mmur.gatk4.paired.joint/filtering/r03.wOutgroups.mac3.FS6_indlist.txt
# VCF=/work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/final/backup/r03.wOutgroups.mac1.FS6.vcf.gz
# bcftools query -l $VCF > analyses/qc/vcf/map2mmur.gatk4.paired.joint/filtering/r03.wOutgroups.mac1.FS6_indlist.txt

## General settings:
FILE_ID=hzproj1.mac1.FS6
VCF_DIR=/work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/final/
PLINK_DIR=/work/jwp37/hybridzone/seqdata/plink/
ATOOLS_DIR=analyses/admixtools/
VCF2PLINK=TRUE
CREATE_INDFILE=TRUE
SUBSET_INDFILE=TRUE
INDS_METADATA=/datacommons/yoderlab/users/jelmer/radseq/metadata/lookup_IDshort.txt

## D-mode:
ATOOLS_MODE="D"
#VCF2PLINK=FALSE
RUN_ID=bySupersite2
INDFILE=analyses/admixtools/input/indfile_$FILE_ID.txt
POPFILE=analyses/admixtools/input/popfile_dstat_hzproj1.$RUN_ID.txt
ID_COLUMN=ID.short
GROUPBY=supersite2
/datacommons/yoderlab/users/jelmer/scripts/genomics/admixtools/admixtools_pip.sh $FILE_ID $RUN_ID $VCF_DIR $PLINK_DIR $ATOOLS_DIR $INDFILE $POPFILE \
	$VCF2PLINK $CREATE_INDFILE $SUBSET_INDFILE $ATOOLS_MODE $INDS_METADATA $ID_COLUMN $GROUPBY

## F4-ratio-mode:
ATOOLS_MODE="F4RATIO"
RUN_ID=bySupersite2
VCF2PLINK=FALSE
INDFILE=analyses/admixtools/input/indfile_$FILE_ID.txt
POPFILE=analyses/admixtools/input/popfile_f4ratio_hzproj1.$RUN_ID.txt
/datacommons/yoderlab/users/jelmer/scripts/genomics/admixtools/admixtools_pip.sh $FILE_ID $RUN_ID $VCF_DIR $PLINK_DIR $ATOOLS_DIR $INDFILE $POPFILE \
	$VCF2PLINK $CREATE_INDFILE $SUBSET_INDFILE $ATOOLS_MODE $INDS_METADATA $ID_COLUMN $GROUPBY



################################################################################
################################################################################
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/hybridzone/scripts/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/scripts/
# rsync -avr --no-perms /home/jelmer/Dropbox/scripts/genomics/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/scripts/genomics/

# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/hybridzone/metadata/* jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/metadata/
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/radseq/metadata/* jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/radseq/metadata/

# rsync -avr /home/jelmer/Dropbox/sc_lemurs/hybridzone/analyses/admixtools/input/* jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/analyses/admixtools/input/
# rsync -avr jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/analyses/admixtools/output/* /home/jelmer/Dropbox/sc_lemurs/hybridzone/analyses/admixtools/output/


################################################################################
## Check inds in pedfile:
#cat /work/jwp37/radseq/seqdata/plink/msp3proj.all.mac3.FS6.ped | cut -f 1

## F3-mode, each species as a pop:
# INDFILE_ID=msp3
# POPFILE_ID=f3stat_$INDFILE_ID
# ATOOLS_MODE=F3
# /datacommons/yoderlab/users/jelmer/scripts/genomics/admixtools/admixtools_pip.sh $FILE_ID $VCF_DIR $PLINK_DIR $ATOOLS_DIR $INDFILE_ID $POPFILE_ID \
#	$VCF2PLINK $CREATE_INDFILE $SUBSET_INDFILE $ATOOLS_MODE $INDS_METADATA
