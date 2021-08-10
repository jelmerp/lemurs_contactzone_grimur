# JOINT GENOTYPING -------------------------------------------------------------

SCR_JGENO=scripts/geno/2_gatk/02_joint-geno_runner.sh
REF_DIR=metadata/reference/mmur/
REF_ID=GCF_000165445.2_Mmur_3.0_genomic_stitched
REF=$REF_DIR/$REF_ID.fasta
SCAFFOLD_FILE=$REF_DIR/$REF_ID.scaffoldList_NC.txt # Only do mapped (NC_) scaffolds + exclude mtDNA
GVCF_DIR=seqdata/vcf/map2mmur.gatk.ind/gvcf
VCF_DIR=seqdata/vcf/map2mmur.gatk.joint/
QC_DIR_VCF=analyses/qc/vcf/map2mmur.gatk.joint
ADD_COMMANDS="none"
MEM_JOB=36
MEM_GATK=24
NCORES=1
SKIP_GENO=FALSE
DP_MEAN=5

FILE_ID=grimurche.og
IDs=($(cut -f 1 /datacommons/yoderlab/users/jelmer/proj/iim/indsel/stacks_popmap/grimurche.og.txt))

sbatch --job-name=jgeno -o slurm.jgeno.pip.$FILE_ID \
    $SCR_JGENO $FILE_ID $SCAFFOLD_FILE $GVCF_DIR $VCF_DIR $QC_DIR_VCF \
    $REF "$ADD_COMMANDS" $MEM_JOB $MEM_GATK $NCORES $SKIP_GENO $DP_MEAN ${IDs[@]}
