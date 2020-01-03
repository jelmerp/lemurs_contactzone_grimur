# cd /datacommons/yoderlab/users/jelmer/hybridzone

BASEDIR=/datacommons/yoderlab/users/jelmer/hybridzone
ID_TYPE=ID.short
KEEP_SET=all
FILE_INDS_KEEP=notany
FILE_LOOKUP=/datacommons/yoderlab/users/jelmer/radseq/metadata/lookup_IDshort.txt
FILE_COLS=/datacommons/yoderlab/users/jelmer/radseq/metadata/colors/colors.species.txt
VCF_DIR=/work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/final/

#FILE_IDS=( r03.wOutgroups.mac1.FS6 r03.wOutgroups.mac3.FS6 r03.wOutgroups.mac1.FS7 r03.wOutgroups.mac3.FS7 )
FILE_IDS=( r03.keepHybs.mac1.FS6 r03.keepHybs.mac3.FS6 r03.keepHybs.mac1.FS7 r03.keepHybs.mac3.FS7 )

for FILE_ID in ${FILE_IDS[@]}
do
	sbatch -p yoderlab,common,scavenger --mem=12G -o slurm.PCA.$FILE_ID \
	scripts/PCA/PCA_adegenet_wrap.sh $FILE_ID $BASEDIR $ID_TYPE $KEEP_SET $FILE_INDS_KEEP $FILE_LOOKUP $FILE_COLS $VCF_DIR
done



################################################################################
################################################################################
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/hybridzone/scripts/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/scripts/
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/radseq/scripts/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/radseq/scripts/
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/hybridzone/metadata/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/metadata/

# rsync -avr jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone//analyses/PCA/dfs/ /home/jelmer/Dropbox/sc_lemurs/hybridzone/analyses/PCA/dfs/