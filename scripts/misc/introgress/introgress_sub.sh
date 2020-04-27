cd /datacommons/yoderlab/users/jelmer/hybridzone


VCF=/work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/final/r03.all.mac1.FS6.vcf.gz
ASSIGNMENT_FILE=analyses/introgress/assignment_files/dummyAssignment_all.txt
ASSIGN_CUTOFF="0.9"
FREQ_CUTOFF="0.5"
OUTPUT_PREFIX="Mtk"
python3 scripts/introgress/vcf2introgress.py $VCF $ASSIGNMENT_FILE $ASSIGN_CUTOFF $FREQ_CUTOFF $OUTPUT_PREFIX



################################################################################
################################################################################
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/hybridzone/scripts/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/scripts/
# rsync -avr --no-perms /home/jelmer/Dropbox/scripts/genomics/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/scripts/genomics/

# rsync -avr /home/jelmer/Dropbox/sc_lemurs/hybridzone/analyses/introgress/assignment_files/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/analyses/introgress/assignment_files/

# python3 -m pip install PyVCF --use
# python3 scripts/introgress/vcf2introgress.py -h

## IndSel:
VCF_IN=/work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/final/r03.all.mac1.FS6.vcf.gz
VCF_OUT=introgress_vcf.vcf.gz
INDS="-mruf007_r01_p3d12"
sbatch -p yoderlab,common,scavenger --mem=32G -o slurm.splitVCF \
	/datacommons/yoderlab/users/jelmer/scripts/genomics/conversion/splitVCF_byIndv.sh $VCF_IN $VCF_OUT $INDS