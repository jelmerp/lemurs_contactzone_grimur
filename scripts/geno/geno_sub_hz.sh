## Outgroup inds: cmed001 cmed003 cmed006 cmed010 ccro003 ccro021 ccro022 mzaz004

################################################################################
#### JOINT GENOTYPING ####
################################################################################
SCR_JGENO=/datacommons/yoderlab/users/jelmer/scripts/geno/gatk/jgeno_pip.sh
REF_DIR=/datacommons/yoderlab/users/jelmer/seqdata/reference/mmur/
REF_ID=GCF_000165445.2_Mmur_3.0_genomic_stitched
REF=$REF_DIR/$REF_ID.fasta
SCAFFOLD_FILE=$REF_DIR/$REF_ID.scaffoldList_NC.txt # Only do mapped (NC_) scaffolds + exclude mtDNA
GVCF_DIR=/datacommons/yoderlab/data/radseq/vcf/map2mmur.gatk.ind/gvcf
VCF_DIR=/datacommons/yoderlab/data/radseq/vcf/map2mmur.gatk.joint/
QC_DIR_VCF=/datacommons/yoderlab/users/jelmer/proj/radseq/analyses/qc/vcf/map2mmur.gatk.joint
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
rsync -avr jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/proj/radseq/analyses/qc/vcf/map2mmur.gatk4.joint/vcftools/*FS6* /home/jelmer/Dropbox/sc_lemurs/radseq/analyses/qc/vcf/map2mmur.gatk4.joint/vcftools/
rsync -avr jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/proj/radseq/analyses/qc/vcf/map2mmur.gatk4.joint/filtering/*filterstats /home/jelmer/Dropbox/sc_lemurs/radseq/analyses/qc/vcf/map2mmur.paired.joint/filtering/

	
################################################################################
################################################################################
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/scripts/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/scripts/
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/proj/radseq/metadata/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/proj/radseq/metadata/

# grep -l -Z 'Done with script igeno_pip.sh' slurm.igeno* | xargs -0 -I{} mv {} radseq/geno/logfiles/
# grep -l -Z 'Done with script gatk2_jointgeno' slurm.jgeno5a* | xargs -0 -I{} mv {} proj/sisp/geno/logfiles/
# grep -l -Z 'CANCELLED' slurm.genoInd* | xargs -0 -I{} mv {} jobcancelled/
