
for MAXMISS in 75 80 85 90 95 
do
	echo "Max missing: $MAXMISS"
	VCF_IN=/work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/intermed/r03.all.rawSNPs.ABHet.vcf
	VCF_OUT=/work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/intermed/customfilter/r03.all.rawSNPs.goodInds.maxmiss$MAXMISS.mac3.vcf.gz
	INDFILE=metadata/r03.all.passedInds.txt
	bcftools query -l /work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/final/r03.all.mac1.FS6.vcf.gz > $INDFILE
	vcftools --vcf $VCF_IN --mac 3 --max-missing 0.$MAXMISS --keep $INDFILE --recode --recode-INFO-all --stdout | gzip > $VCF_OUT
done


# rsync -avr jwp37@dcc-slogin-02.oit.duke.edu:/work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/intermed/customfilter/ /home/jelmer/Dropbox/sc_lemurs/hybridzone/seqdata/vcf/


## Fastq stats:
SCRIPT_STATS_FASTQ=/datacommons/yoderlab/users/jelmer/radseq/scripts/fastq_process/fastq.stats.sh
DIR_STATS=/datacommons/yoderlab/users/jelmer/hybridzone/analyses/qc/fastq/byInd
for FASTQ in /work/jwp37/radseq/seqdata/fastq/r03/processed/*fastq.gz
do
	FASTQ_SHORT=$(basename -s .fastq.gz $FASTQ)
	echo -e "##### fastq.process_pip.sh: Step 4 -- getting stats for $FASTQ_SHORT..."
	sbatch -p yoderlab,common,scavenger -o slurm.fastqstats.$FASTQ_SHORT $SCRIPT_STATS_FASTQ $FASTQ $DIR_STATS
	printf "\n"
done

# rsync -avr jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/analyses/qc/fastq/byInd/ /home/jelmer/Dropbox/sc_lemurs/hybridzone/analyses/qc/fastq/byInd/