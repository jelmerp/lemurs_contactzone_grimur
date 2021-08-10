
## ID:
#VCF_ID=r03.all.mac1.FS6
#VCF_ID=r03.all.mac3.FS6
#VCF_ID=hzproj2.mac1.FS6
VCF_ID=hzproj1.mac1.FS6

## Set-up:
cd /home/jelmer/Dropbox/sc_lemurs/hybridzone
SCAFFOLDS=$( cat /home/jelmer/Dropbox/sc_lemurs/other/seqdata_misc/reference/mmur/scaffolds_stitched.txt )
VCF_IN=/home/jelmer/Dropbox/sc_lemurs/hybridzone/seqdata/vcf/$VCF_ID.vcf
DIR_VCF_OUT=/home/jelmer/Dropbox/sc_lemurs/hybridzone/seqdata/vcf/by_scaffold/$VCF_ID/
mkdir -p $DIR_VCF_OUT

## Run:
gunzip $VCF_IN.gz
bgzip $VCF_IN
tabix -p vcf $VCF_IN.gz

for SCAFFOLD in ${SCAFFOLDS[@]}
do
	OUTFILE=$DIR_VCF_OUT/$SCAFFOLD.vcf.gz
	echo "$SCAFFOLD $OUTFILE"
	bcftools view $VCF_IN.gz $SCAFFOLD > $OUTFILE
done
