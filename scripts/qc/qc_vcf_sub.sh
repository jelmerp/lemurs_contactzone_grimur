VCF_DIR=/work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/final/
QC_DIR=analyses/qc/vcf/map2mmur.gatk4.paired.joint/
RUN_BCF=TRUE
RUN_BCF_BY_IND=FALSE
MEM=4

FILE_ID=r03.keepHybs.mac1.FS6
sbatch --mem ${MEM}G -p common,yoderlab,scavenger -o slurm.qcVCF.$FILE_ID \
scripts/qc/qc_vcf.sh $FILE_ID $VCF_DIR $QC_DIR $RUN_BCF $RUN_BCF_BY_IND



################################################################################
################################################################################
# rsync -r --verbose jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/radseq/analyses/qc/vcf/map2msp3.gatk4.paired.joint/ /home/jelmer/Dropbox/sc_lemurs/radseq/analyses/qc/vcf/
