cd /datacommons/yoderlab/users/jelmer

## General options:
SCR_PIP=/datacommons/yoderlab/users/jelmer/scripts/geno/stacks/stacks_pip.sh
ADD_OPS_GSTACKS=""
ADD_OPS_POPSTACKS="--min-samples-overall 0.80 --max-obs-het 0.70 --treemix --radpainter --fasta-samples --hwe --fstats --vcf"
BAMDIR=/datacommons/yoderlab/data/radseq/bam/map2mmur/final_merged/
BAMSUFFIX=".sort.MQ30.dedup.bam"
STACKSDIR=/datacommons/yoderlab/users/jelmer/proj/hybridzone/seqdata/stacks/
NCORES=8
POPMAP_DIR=/datacommons/yoderlab/users/jelmer/proj/hybridzone/metadata/indsel/stacks

## Gstacks + popstacks for all individuals:
GSTACKS_ID=hzproj
POPSTACKS_ID=all
POPMAP_ALL=$POPMAP_DIR/$GSTACKS_ID.txt
cat $POPMAP_ALL
POPMAP_SEL=$POPMAP_DIR/$GSTACKS_ID.$POPSTACKS_ID.txt
cat $POPMAP_ALL > $POPMAP_SEL
TO_SKIP=""
$SCR_PIP $GSTACKS_ID $POPSTACKS_ID $STACKSDIR $BAMDIR $BAMSUFFIX \
	$POPMAP_ALL $POPMAP_SEL "$ADD_OPS_GSTACKS" "$ADD_OPS_POPSTACKS" $NCORES $TO_SKIP
	
## popstacks for mur3gri2ruf:
GSTACKS_ID=hzproj
POPSTACKS_ID=mur3gri2ruf
POPMAP_ALL=indsel/stacks_popmap/$GSTACKS_ID.txt
POPMAP_SEL=indsel/stacks_popmap/$GSTACKS_ID.$POPSTACKS_ID.txt
TO_SKIP="-G"
$SCR_PIP $GSTACKS_ID $POPSTACKS_ID $STACKSDIR $BAMDIR $BAMSUFFIX \
	$POPMAP_ALL $POPMAP_SEL "$ADD_OPS_GSTACKS" "$ADD_OPS_POPSTACKS" $NCORES $TO_SKIP
	
## popstacks for mur3gri2c:
GSTACKS_ID=hzproj
POPSTACKS_ID=mur3gri2c
POPMAP_ALL=indsel/stacks_popmap/$GSTACKS_ID.txt
POPMAP_SEL=indsel/stacks_popmap/$GSTACKS_ID.$POPSTACKS_ID.txt
TO_SKIP="-G"
$SCR_PIP $GSTACKS_ID $POPSTACKS_ID $STACKSDIR $BAMDIR $BAMSUFFIX \
	$POPMAP_ALL $POPMAP_SEL "$ADD_OPS_GSTACKS" "$ADD_OPS_POPSTACKS" $NCORES $TO_SKIP
	
#scp jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/proj/hybridzone/seqdata/stacks/stacks_*txt /home/jelmer/Dropbox/sc_lemurs/proj/hybridzone/qc/stacks/




################################################################################
#### PROCESS VCF ####
################################################################################
VCF_DIR=$STACKSDIR/$STACKS_ID/
VCF_IN=$VCF_DIR/populations.snps.vcf.gz
VCF_OUT=$VCF_DIR/populations.snps.mac3.vcf
MAC=3
vcftools --vcf $VCF_IN.vcf.gz --recode --recode-INFO-all --mac $MAC --stdout > $VCF_OUT
### TRY OTHER FILTERS


################################################################################
#### PROCESS FASTA ####
################################################################################
#[ID].markers.tsv # Can extract missing data props from this
#[ID].sumstats.tsv # Column 8 has nr of inds per site *per pop*
#[ID].hapstats.tsv # Column 8 has nr of inds per locus *per pop*
/datacommons/yoderlab/users/jelmer/proj/murclade/geno/stacks/murclade/all/murclade.all.markers.tsv
/datacommons/yoderlab/users/jelmer/proj/murclade/geno/stacks/murclade/all/murclade.all.sumstats.tsv


################################################################################
################################################################################
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/scripts/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/scripts/
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/radseq/metadata/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/radseq/metadata/
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/proj/hybridzone/metadata/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/proj/hybridzone/metadata/

#rsync -avr --no-perms jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/proj/hybridzone/seqdata/stacks/hzproj/all/fst/hzproj.all.phistats_mur-E-mruf.tsv /home/jelmer/Dropbox/sc_lemurs/proj/iim/locusstats/stacks/

# http://catchenlab.life.illinois.edu/stacks/manual/
# http://catchenlab.life.illinois.edu/stacks/comp/populations.php