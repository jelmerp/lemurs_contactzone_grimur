## Script to separate loci in G-PhoCS input sequence file based on whether
## they are genic or intergenic.


################################################################################
#### SET-UP ####
################################################################################
cd /home/jelmer/Dropbox/sc_lemurs/proj/hybridzone/

## Input files:
INPUT_ORG=analyses/gphocs/input/r03.wOutgroups.hz.mur3gri2c.gphocsInput.txt
GENESBED=/home/jelmer/Dropbox/sc_lemurs/seqdata/reference/mmur/GCF_000165445.2_Mmur_3.0_genes.bed
EXONBED=/home/jelmer/Dropbox/sc_lemurs/seqdata/reference/mmur/GCF_000165445.2_Mmur_3.0_exons.bed

## Output files:
LOCUSLIST=analyses/gphocs/input/genic/mur3gri2c.locuslist.txt
LOCUSLIST2=analyses/gphocs/input/genic/mur3gri2c.locuslist.onlyNC.txt
LOCUSBED_EXONIC=analyses/gphocs/input/genic/mur3gri2c.locusbed.onlyNC.exonic.txt
LOCUSBED_GENIC=analyses/gphocs/input/genic/mur3gri2c.locusbed.onlyNC.genic.txt
LOCUSBED_INTERGENIC=analyses/gphocs/input/genic/mur3gri2c.locusbed.onlyNC.intergenic.txt
INPUT_EXONIC=analyses/gphocs/input/genic/mur3gri2c.onlyNC.exonic.gphocsinput.txt
INPUT_GENIC=analyses/gphocs/input/genic/mur3gri2c.onlyNC.genic.gphocsinput.txt
INPUT_INTERGENIC=analyses/gphocs/input/genic/mur3gri2c.onlyNC.intergenic.gphocsinput.txt


################################################################################
#### CREATE BEDFILE WITH LOCI ####
################################################################################
grep "fa" $INPUT_ORG | sed -r 's/\./\t/' | sed 's/.fa//' | sed -r 's/-/\t/' | sed -r 's/:/\t/' | sed -r 's/ /\t/g' | cut -f 2,3,4,5,6,1 > $LOCUSLIST.tmp  #sed -r 's/^[0-9]+\.//'
paste $LOCUSLIST.tmp $LOCUSLIST.tmp | cut -f2,3,4,5,6,7 > $LOCUSLIST
grep -v "Super_Scaffold" $LOCUSLIST > $LOCUSLIST2 # REMOVE SUPER_SCAFFOLD LOCI! (n=121)
#cut -f 4 $LOCUSLIST2 | sort | uniq # All have 15 inds

bedtools intersect -u -a $LOCUSLIST2 -b $EXONBED > $LOCUSBED_EXONIC # Get bedfile with exonic loci
bedtools intersect -u -a $LOCUSLIST2 -b $GENESBED > $LOCUSBED_GENIC # Get bedfile with genic loci
bedtools intersect -v -a $LOCUSLIST2 -b $GENESBED > $LOCUSBED_INTERGENIC # Get bedfile with intergenic loci

wc -l $LOCUSBED_EXONIC
wc -l $LOCUSBED_GENIC
wc -l $LOCUSBED_INTERGENIC


################################################################################
#### EXTRACT LOCI FROM INPUT FILE FOR GENIC/INTERGENIC ####
################################################################################
## Exonnic:
cat $LOCUSBED_EXONIC | wc -l > $INPUT_EXONIC
echo -en '\n' >> $INPUT_EXONIC
cut -f 6 $LOCUSBED_EXONIC | while read LOCUSID
do
	grep -A 15 "^${LOCUSID}\." $INPUT_ORG >> $INPUT_EXONIC
done

## Genic:
cat $LOCUSBED_GENIC | wc -l > $INPUT_GENIC
echo -en '\n' >> $INPUT_GENIC
cut -f 6 $LOCUSBED_GENIC | while read LOCUSID
do
	grep -A 15 "^${LOCUSID}\." $INPUT_ORG >> $INPUT_GENIC
done

## Intergenic:
cat $LOCUSBED_INTERGENIC | wc -l > $INPUT_INTERGENIC
echo -en '\n' >> $INPUT_INTERGENIC
cut -f 6 $LOCUSBED_INTERGENIC | while read LOCUSID
do
	grep -A 15 "^${LOCUSID}\." $INPUT_ORG >> $INPUT_INTERGENIC
done

## Check if all loci were recovered:
#grep " 15 " $INPUT_GENIC | sed -r 's/^([0-9]+\.).*/\1/' | sort | uniq -d

wc -l $LOCUSBED_EXONIC
grep " 15 " $INPUT_EXONIC | wc -l

wc -l $LOCUSBED_GENIC
grep " 15 " $INPUT_GENIC | wc -l

wc -l $LOCUSBED_INTERGENIC
grep " 15 " $INPUT_INTERGENIC | wc -l

## Remove tmp files:
rm analyses/gphocs/input/genic/*tmp


################################################################################
#### CONVERT TO BPP INPUT ####
################################################################################
INPUT_BPP_EXONIC=analyses/bpp/input/genic/mur3gri2c.onlyNC.exonic.txt
INPUT_BPP_GENIC=analyses/bpp/input/genic/mur3gri2c.onlyNC.genic.txt
INPUT_BPP_INTERGENIC=analyses/bpp/input/genic/mur3gri2c.onlyNC.intergenic.txt

sed 's/.*fa //' $INPUT_EXONIC > $INPUT_BPP_EXONIC
sed 's/.*fa //' $INPUT_GENIC > $INPUT_BPP_GENIC
sed 's/.*fa //' $INPUT_INTERGENIC > $INPUT_BPP_INTERGENIC
