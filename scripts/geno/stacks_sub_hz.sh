cd ~/dc/proj/hybridzone

## General options:
ADD_OPS_GSTACKS=""
ADD_OPS_POPSTACKS="--min-samples-overall 0.80 --max-obs-het 0.70 --treemix --radpainter --fasta-samples --hwe --fstats --vcf"
BAMDIR=seqdata/radseq/bam/map2mmur/final_merged/
BAMSUFFIX=".sort.MQ30.dedup.bam"
STACKSDIR=seqdata/stacks/
NCORES=8

## Gstacks + popstacks for all individuals:
GSTACKS_ID=hzproj
POPSTACKS_ID=all
POPMAP_GSTACKS=metadata/indsel/stacks/$GSTACKS_ID.txt
POPMAP_POPSTACKS=metadata/indsel/stacks/$GSTACKS_ID.$POPSTACKS_ID.txt
cp $POPMAP_GSTACKS $POPMAP_POPSTACKS
TO_SKIP=""
scripts/genomics/geno/stacks/stacks_pip.sh $GSTACKS_ID $POPSTACKS_ID $STACKSDIR $BAMDIR \
	$BAMSUFFIX $POPMAP_GSTACKS $POPMAP_POPSTACKS "$ADD_OPS_GSTACKS" "$ADD_OPS_POPSTACKS" \
	$NCORES $TO_SKIP
	
## popstacks for mur3gri2ruf:
GSTACKS_ID=hzproj
POPSTACKS_ID=mur3gri2ruf
POPMAP_ALL=indsel/stacks_popmap/$GSTACKS_ID.txt
POPMAP_SEL=indsel/stacks_popmap/$GSTACKS_ID.$POPSTACKS_ID.txt
TO_SKIP="-G"
scripts/genomics/geno/stacks/stacks_pip.sh $GSTACKS_ID $POPSTACKS_ID $STACKSDIR $BAMDIR \
	$BAMSUFFIX $POPMAP_GSTACKS $POPMAP_POPSTACKS "$ADD_OPS_GSTACKS" "$ADD_OPS_POPSTACKS" \
	$NCORES $TO_SKIP
	
## popstacks for mur3gri2c:
GSTACKS_ID=hzproj
POPSTACKS_ID=mur3gri2c
POPMAP_ALL=indsel/stacks_popmap/$GSTACKS_ID.txt
POPMAP_SEL=indsel/stacks_popmap/$GSTACKS_ID.$POPSTACKS_ID.txt
TO_SKIP="-G"
scripts/genomics/geno/stacks/stacks_pip.sh $GSTACKS_ID $POPSTACKS_ID $STACKSDIR $BAMDIR \
	$BAMSUFFIX $POPMAP_GSTACKS $POPMAP_POPSTACKS "$ADD_OPS_GSTACKS" "$ADD_OPS_POPSTACKS" \
	$NCORES $TO_SKIP
	
#scp dcc:/dc/proj/hybridzone/seqdata/stacks/stacks_*txt ~/l/proj/hybridzone/qc/stacks/


################################################################################
################################################################################
rsync -avr --no-perms ~/sg dcc:dc/scripts/genomics/
rsync -avr ~/l/proj/hybridzone/metadata dcc:dc/proj/hybridzone/metadata/

rsync -avr -dcc:/datacommons/yoderlab/users/jelmer/proj/hybridzone/seqdata/stacks/hzproj/all/fst/hzproj.all.phistats_mur-E-mruf.tsv /home/jelmer/Dropbox/sc_lemurs/proj/iim/locusstats/stacks/