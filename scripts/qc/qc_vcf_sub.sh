VCF_DIR=seqdata/vcf/map2mmur.gatk4.paired.joint/final/
QC_DIR=analyses/qc/vcf/map2mmur.gatk4.paired.joint/
RUN_BCF=TRUE
RUN_BCF_BY_IND=FALSE
MEM=4

FILE_ID=r03.keepHybs.mac1.FS6
sbatch --mem ${MEM}G -p common,yoderlab,scavenger -o slurm.qcVCF.$FILE_ID \
scripts/qc/qc_vcf.sh $FILE_ID $VCF_DIR $QC_DIR $RUN_BCF $RUN_BCF_BY_IND



################################################################################
################################################################################
# rsync -avr dcc:~/dc/proj/hybridzone/analyses/qc/vcf/map2mmur.gatk4.paired.joint/ ~/Dropbox/sc_lemurs/proj/hybridzone/analyses/qc/vcf/

VCF=seqdata/vcf/gatk/final/hzproj1.mac1.FS6.vcf.gz
vcftools --gzvcf $VCF --idepth