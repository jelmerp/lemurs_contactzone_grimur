## Outgroup inds: cmed001 cmed003 cmed006 cmed010 ccro003 ccro021 ccro022 mzaz004

cd ~/dc/proj/hybridzone

################################################################################
#### JOINT GENOTYPING ####
################################################################################
SCR_JGENO=scripts/genomics/geno/gatk/jgeno_pip.sh
REF_DIR=seqdata/reference/mmur/
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
IDs=( $(cut -f 1 /datacommons/yoderlab/users/jelmer/proj/iim/indsel/stacks_popmap/grimurche.og.txt) )

sbatch -p yoderlab,common,scavenger --job-name=jgeno -o slurm.jgeno.pip.$FILE_ID \
	$SCR_JGENO $FILE_ID $SCAFFOLD_FILE $GVCF_DIR $VCF_DIR $QC_DIR_VCF \
	$REF "$ADD_COMMANDS" $MEM_JOB $MEM_GATK $NCORES $SKIP_GENO $DP_MEAN ${IDs[@]}
	
## Process QC:
rsync -avr dcc:dc/proj/radseq/analyses/qc/vcf/map2mmur.gatk4.joint/vcftools/*FS6* ~/l/proj/radseq/analyses/qc/vcf/map2mmur.gatk4.joint/vcftools/
rsync -avr dcc:dc/proj/radseq/analyses/qc/vcf/map2mmur.gatk4.joint/filtering/*filterstats ~/l/radseq/analyses/qc/vcf/map2mmur.paired.joint/filtering/

################################################################################
################################################################################
rsync -avr --no-perms ~/sg/ dcc:dc/scripts/genomics/
rsync -avr ~/l/proj/radseq/metadata/ dcc:dc/proj/radseq/metadata/

# ln -s ~/dc/scripts/genomics scripts/genomics
# ln -s ~/dc/seqdata/reference seqdata/reference