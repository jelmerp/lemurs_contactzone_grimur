## Script to separate loci in G-PhoCS input sequence file based on whether
## they are genic or intergenic.


################################################################################
#### SET-UP ####
################################################################################
cd /home/jelmer/Dropbox/sc_lemurs/proj/hybridzone/

## Input files:
INPUT_ORG=analyses/gphocs/input/r03.wOutgroups.hz.mur3gri2c.gphocsInput.txt
GENESBED=/home/jelmer/Dropbox/sc_lemurs/seqdata/reference/mmur/GCF_000165445.2_Mmur_3.0_genes.bed

## Output files:
LOCUSLIST=analyses/gphocs/input/genic/mur3gri2c.locuslist.txt
LOCUSLIST2=analyses/gphocs/input/genic/mur3gri2c.locuslist.onlyNC.txt
LOCUSBED_GENIC=analyses/gphocs/input/genic/mur3gri2c.locusbed.onlyNC.genic.txt
LOCUSBED_INTERGENIC=analyses/gphocs/input/genic/mur3gri2c.locusbed.onlyNC.intergenic.txt
INPUT_GENIC=analyses/gphocs/input/genic/mur3gri2c.onlyNC.genic.gphocsinput.txt
INPUT_INTERGENIC=analyses/gphocs/input/genic/mur3gri2c.onlyNC.intergenic.gphocsinput.txt


################################################################################
#### CREATE BEDFILE WITH LOCI ####
################################################################################
grep "fa" $INPUT_ORG | sed -r 's/\./\t/' | sed 's/.fa//' | sed -r 's/-/\t/' | sed -r 's/:/\t/' | sed -r 's/ /\t/g' | cut -f 2,3,4,5,6,1 > $LOCUSLIST.tmp  #sed -r 's/^[0-9]+\.//'
paste $LOCUSLIST.tmp $LOCUSLIST.tmp | cut -f2,3,4,5,6,7 > $LOCUSLIST
grep -v "Super_Scaffold" $LOCUSLIST > $LOCUSLIST2 # REMOVE SUPER_SCAFFOLD LOCI! (n=121)
#cut -f 4 $LOCUSLIST2 | sort | uniq # All have 15 inds
bedtools intersect -v -a $LOCUSLIST2 -b $GENESBED > $LOCUSBED_GENIC # Get bedfile with genic loci
bedtools intersect -wa -a $LOCUSLIST2 -b $GENESBED > $LOCUSBED_INTERGENIC # Get bedfile with genic loci


################################################################################
#### EXTRACT LOCI FROM INPUT FILE FOR GENIC/INTERGENIC ####
################################################################################
## Genic:
> $INPUT_GENIC
cut -f 6 $LOCUSBED_GENIC | while read LOCUSID
do
	grep -A 15 "^${LOCUSID}\." $INPUT_ORG >> $INPUT_GENIC
done

## Intergenic:
> $INPUT_INTERGENIC
cut -f 6 $LOCUSBED_INTERGENIC | while read LOCUSID
do
	grep -A 15 "^${LOCUSID}\." $INPUT_ORG >> $INPUT_INTERGENIC
done

## Check if all loci were recovered:
wc -l $LOCUSBED_GENIC
grep " 15 " $INPUT_GENIC | wc -l
#grep " 15 " $INPUT_GENIC | sed -r 's/^([0-9]+\.).*/\1/' | sort | uniq -d

wc -l $LOCUSBED_INTERGENIC
grep " 15 " $INPUT_INTERGENIC | wc -l

