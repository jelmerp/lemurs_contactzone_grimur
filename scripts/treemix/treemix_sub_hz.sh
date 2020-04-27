################################################################################
##### RUN TREEMIX #####
################################################################################
## General settings:
VCF_DIR=/datacommons/yoderlab/users/jelmer/proj/hybridzone/seqdata/vcf/map2mmur.gatk.joint/final/
PREP_INPUT=TRUE
MINMIG=0
MAXMIG=10
TREEMIX_DIR=/datacommons/yoderlab/users/jelmer/proj/sisp/analyses/treemix/
INDS_METADATA=/datacommons/yoderlab/users/jelmer/radseq/metadata/lookup_IDshort.txt
#SELECT_BY_COLUMN="pop"; ROOT=mruf
SELECT_BY_COLUMN="pop2"; ROOT=ruf

## Run:
FILE_IDS=( hzproj1.mac1.FS6 hzproj1.mac3.FS6 )
for FILE_ID in ${FILE_IDS[@]}
do
	#FILE_ID=hzproj1.mac1.FS6
	echo "#### File ID: $FILE_ID"
	/datacommons/yoderlab/users/jelmer/scripts/genomics/treemix/treemix_pip.sh $FILE_ID $VCF_DIR $PREP_INPUT $MINMIG $MAXMIG $ROOT $TREEMIX_DIR $INDS_METADATA $SELECT_BY_COLUMN
done


	
################################################################################
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/hybridzone/scripts/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/scripts/
# rsync -avr --no-perms /home/jelmer/Dropbox/scripts/genomics/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/scripts/genomics/

# rsync -avr /home/jelmer/Dropbox/sc_lemurs/radseq/metadata/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/radseq/metadata/
# rsync -avr /home/jelmer/Dropbox/sc_lemurs/hybridzone/metadata/* jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/metadata/
# rsync -avr /home/jelmer/Dropbox/sc_lemurs/hybridzone/analyses/treemix/popfiles/* jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/analyses/treemix/popfiles/

# rsync -avr jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/analyses/treemix/output/* /home/jelmer/Dropbox/sc_lemurs/hybridzone/analyses/treemix/output
