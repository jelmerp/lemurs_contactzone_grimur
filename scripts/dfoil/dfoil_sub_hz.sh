#mkdir -p analyses/dfoil/input analyses/dfoil/output analyses/dfoil/summaries

################################################################################
##### SET-UP #####
################################################################################
FILE_ID=r03.wOutgroups.indsel.mac3.FS6
#FILE_ID=r03.wOutgroups.indsel.mac1.FS6
VCF_DIR=/work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/final
FASTA_DIR=/work/jwp37/hybridzone/seqdata/fasta
FASTA=$FASTA_DIR/$FILE_ID.fasta
DFOIL_INDIR=analyses/dfoil/input
DFOIL_OUTDIR=analyses/dfoil/output
OUTGROUP=mruf


################################################################################
##### STEP 1 - CONVERT VCF TO FASTA #####
################################################################################
SCAFFOLD=ALL
sbatch --mem 8G -p yoderlab,common,scavenger -o slurm.vcf2fasta.$FILE_ID.txt \
	/datacommons/yoderlab/users/jelmer/scripts/conversion/vcf2fasta.sh $FILE_ID $VCF_DIR $FASTA_DIR $SCAFFOLD

	
################################################################################
##### STEP 2 - PREP FASTA FOR DFOIL #####
################################################################################
## All sequences belonging to same focal unit (species, pop, whatever level the test is run on)
## should have the exact same name in the fasta file as specified for the fasta2dfoil script.
## Also, only the 5 pops/species that are being actively tested can be present in the fasta.

## Change names:
cp $FASTA $FASTA.tmp
cat $FASTA.tmp | \
	sed 's,mgri084\|mgri085\|mgri086\|mgri087\|mgri088\|mgri090\|mgri093\|mgri094\|mgri095\|mgri097\|mgri098\|mgri103\|mgri104,mgri_hz,' | \
	sed 's,mgri001\|mgri003\|mgri004\|mgri036\|mgri048\|mgri049\|mgri075\|mgri076\|mgri077\|mgri078\|mgri079\|mgri080\|mgri082,mgri_se,' | \
	sed 's,mgri005\|mgri006\|mgri007\|mgri008\|mgri037\|mgri038\|mgri040\|mgri041\|mgri043\|mgri044\|mgri045\|mgri046\|mgri047\|mgri050\|mgri051,mgri_sw,' | \
	sed 's,mgan007\|mgan008\|mgan010\|mgan011\|mgan014\|mgan016\|mgan017\|mgan018\|mgan019\|mgan021\|mgan022\|mgan023,mmur_gan,' | \
	sed 's,mmur045\|mmur047\|mmur048\|mmur049\|mmur050\|mmur051\|mmur052\|mmur053\|mmur055\|mmur056\|mmur057\|mmur061\|mmur063\|mmur064\|mmur065\|mmur066\|mmur068\|mmur069\|mmur070,mmur_hz,' | \
	sed 's,mmur037\|mmur038\|mmur039\|mmur040\|mmur041\|mmur042\|mmur043\|mmur044,mmur_se,' | \
	sed 's,mmur001\|mmur002\|mmur004\|mmur005\|mmur006\|mmur008\|mmur009\|mmur010\|mmur012\|mmur013\|mmur014,mmur_w,' | \
	sed 's,mruf007_r01_p3d12,mruf,' > $FASTA
grep ">" $FASTA | sort | uniq

## Remove inds:
cat $FASTA | awk '/^>/ {P=index($0,"mgri_se")==0} {if(P) print} ' | awk '/^>/ {P=index($0,"mmur_se")==0} {if(P) print} ' | \
	awk '/^>/ {P=index($0,"mmur_gan")==0} {if(P) print} ' > $FASTA.HZ
grep ">" $FASTA.HZ | sort | uniq

cat $FASTA | awk '/^>/ {P=index($0,"mgri_hz")==0} {if(P) print} ' | awk '/^>/ {P=index($0,"mmur_hz")==0} {if(P) print} ' | \
	awk '/^>/ {P=index($0,"mmur_gan")==0} {if(P) print} ' > $FASTA.SE
grep ">" $FASTA.SE | sort | uniq

cat $FASTA | awk '/^>/ {P=index($0,"mgri_se")==0} {if(P) print} ' | awk '/^>/ {P=index($0,"mmur_se")==0} {if(P) print} ' | \
	awk '/^>/ {P=index($0,"mmur_w")==0} {if(P) print} ' > $FASTA.HZgan
grep ">" $FASTA.HZgan | sort | uniq

cat $FASTA | awk '/^>/ {P=index($0,"mgri_hz")==0} {if(P) print} ' | awk '/^>/ {P=index($0,"mmur_hz")==0} {if(P) print} ' | \
	awk '/^>/ {P=index($0,"mmur_w")==0} {if(P) print} ' > $FASTA.SEgan
grep ">" $FASTA.SEgan | sort | uniq

cat $FASTA | awk '/^>/ {P=index($0,"mgri_se")==0} {if(P) print} ' | awk '/^>/ {P=index($0,"mmur_se")==0} {if(P) print} ' | \
	awk '/^>/ {P=index($0,"mmur_hz")==0} {if(P) print} ' > $FASTA.gan
grep ">" $FASTA.gan | sort | uniq


################################################################################
##### STEP 3 - PREP DFOIL INPUT FROM FASTA #####
################################################################################
#module load python/2.7.11
FASTA2DFOIL=/datacommons/yoderlab/programs/dfoil/fasta2dfoil.py

POPCOMB_HZ="mgri_sw,mgri_hz,mmur_hz,mmur_w,$OUTGROUP"
INFILE_HZ=$DFOIL_INDIR/$FILE_ID.HZ.dfoil.in
$FASTA2DFOIL $FASTA.HZ --out $INFILE_HZ --names $POPCOMB_HZ

POPCOMB_SE="mgri_sw,mgri_se,mmur_se,mmur_w,$OUTGROUP"
INFILE_SE=$DFOIL_INDIR/$FILE_ID.SE.dfoil.in
$FASTA2DFOIL $FASTA.SE --out $INFILE_SE --names $POPCOMB_SE

POPCOMB_HZ_GAN="mgri_sw,mgri_hz,mmur_hz,mmur_gan,$OUTGROUP" 
INFILE_HZ_GAN=$DFOIL_INDIR/$FILE_ID.HZgan.dfoil.in
$FASTA2DFOIL $FASTA.HZgan --out $INFILE_HZ_GAN --names $POPCOMB_HZ_GAN

POPCOMB_SE_GAN="mgri_sw,mgri_se,mmur_se,mmur_gan,$OUTGROUP"
INFILE_SE_GAN=$DFOIL_INDIR/$FILE_ID.SEgan.dfoil.in
$FASTA2DFOIL $FASTA.SEgan --out $INFILE_SE_GAN --names $POPCOMB_SE_GAN

POPCOMB_GAN="mgri_sw,mgri_hz,mmur_gan,mmur_w,$OUTGROUP"
INFILE_GAN=$DFOIL_INDIR/$FILE_ID.gan.dfoil.in
$FASTA2DFOIL $FASTA.gan --out $INFILE_GAN --names $POPCOMB_GAN


################################################################################
##### STEP 4 - RUN DFOIL #####
################################################################################
module load Python/3.6.4
DFOIL=/datacommons/yoderlab/programs/dfoil/dfoil.py

## Normal mode:
for SUFFIX in HZ SE HZgan SEgan gan
do
	echo "##### Suffix: $SUFFIX"
	INFILE=$DFOIL_INDIR/$FILE_ID.$SUFFIX.dfoil.in
	OUTFILE=$DFOIL_OUTDIR/$FILE_ID.$SUFFIX.dfoil.out
	$DFOIL --infile $INFILE --out $OUTFILE --mode dfoil
	printf "\n\n\n"
done

## Alt mode:
for SUFFIX in HZ SE HZgan SEgan gan
do
	echo "##### Suffix: $SUFFIX"
	INFILE=$DFOIL_INDIR/$FILE_ID.$SUFFIX.dfoil.in
	OUTFILE=$DFOIL_OUTDIR/$FILE_ID.$SUFFIX.altMode.dfoil.out
	$DFOIL --infile $INFILE --out $OUTFILE --mode dfoilalt
	printf "\n\n\n"
done



################################################################################
################################################################################
# rsync -avr /home/jelmer/Dropbox/sc_lemurs/hybridzone/scripts/* jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/scripts/
# rsync -avr /home/jelmer/Dropbox/sc_lemurs/scripts/* jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/scripts/

# rsync -avr jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/analyses/dfoil/output/* /home/jelmer/Dropbox/sc_lemurs/hybridzone/dfoil/output/


################################################################################
################################################################################
## Get regular dstats:
# $DFOIL --infile $DFOIL_INFILE.noCfus --out $DFOIL_OUTFILE.noCfus.dstats --mode dstat

## Edit fasta:
# cat $FASTA | awk '/^>/ {P=index($0,"mrav01")==0} {if(P) print} ' | awk '/^>/ {P=index($0,"msam01")==0} {if(P) print} ' | \
#	awk '/^>/ {P=index($0,"mmac01")==0} {if(P) print} ' | awk '/^>/ {P=index($0,"mtav01")==0} {if(P) print} ' | \
#	awk '/^>/ {P=index($0,"mzaz01")==0} {if(P) print} ' > $FASTA.dfoil

# TAXA=mleh01,mmit01,mmyo01,mber01,$OUTGROUP
# scripts/dfoil/dfoil_run.sh $FILE_ID $FASTA $DFOIL_INFILE $DFOIL_OUTFILE "$TAXA"

## Dfoil-alt mode:
# $DFOIL --infile $DFOIL_INFILE --out $OUTFILE_NO_MMIT.alt --mode dfoilalt

#SEQKIT=/datacommons/yoderlab/programs/seqkit # https://github.com/shenwei356/seqkit