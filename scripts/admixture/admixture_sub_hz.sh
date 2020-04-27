module load bcftools # Bcftools happens to be on cluster module system and can be loaded

VCF_DIR=/datacommons/yoderlab/users/jelmer/proj/hybridzone/seqdata/vcf/map2mmur.gatk.joint/final # Dir with existing VCF file(s)
PLINK_DIR=/work/jwp37/hybridzone/seqdata/plink # Dir for PLINK files (to be produced)
OUTDIR=analyses/admixture/output/ # Dir for Admixture files (to be produced)
MAF=0 # Minor Allele Frequency (normally set to 0, meaning that low-frequency alleles will not be removed)
LD_MAX=1 # Max Linkage Disequilibrium (LD) (normally set to 1, meaning to no "LD-pruning" will be done)
NCORES=1 # Number of cores
INDFILE=metadata/indsel/indsel.admixture.txt # Optional: file with individuals to select (i.e. discard those not in the file

FILE_IDS=( r03.wOutgroups.mac1.FS6 r03.wOutgroups.mac3.FS6 ) # File IDs for VCF files, VCF files should be: $VCF_DIR/$FILE_ID.vcf.gz

for FILE_ID in ${FILE_IDS[@]}
do
	# FILE_ID=r03.wOutgroups.mac1.FS6 # Can also assign $FILE_ID like so to avoid loop
	echo -e "\n##### File ID: $FILE_ID"
	bcftools query -l $VCF_DIR/$FILE_ID.vcf.gz | grep -v "mruf" > $INDFILE # Send list of inds minus "mruf" (outgroup) to separate textfile

	sbatch -p yoderlab,common,scavenger --mem 8G -o slurm.admixture.pip.$FILE_ID \
	/datacommons/yoderlab/users/jelmer/scripts/genomics/admixture/admixture_pip.sh $FILE_ID $VCF_DIR $PLINK_DIR $OUTDIR $MAF $LD_MAX $NCORES $INDFILE
done


################################################################################
################################################################################
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/hybridzone/scripts/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/scripts/
# rsync -avr --no-perms /home/jelmer/Dropbox/scripts/genomics/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/scripts/genomics/
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/hybridzone/metadata/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/metadata/

# rsync -avr jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/analyses/admixture/output/* /home/jelmer/Dropbox/sc_lemurs/hybridzone/analyses/admixture/output/
