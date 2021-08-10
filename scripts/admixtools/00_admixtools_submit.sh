## General settings:
FILE_ID=hzproj1.mac1.FS6
VCF_DIR=results/geno/vcf/map2mmur.gatk4.paired.joint/final/
PLINK_DIR=results/geno/conversions/plink/
ATOOLS_DIR=analyses/admixtools/
VCF2PLINK=TRUE
CREATE_INDFILE=TRUE
SUBSET_INDFILE=TRUE
INDS_METADATA=metadata/hzlookup_bysample.txt

## D-mode:
ATOOLS_MODE="D"
VCF2PLINK=TRUE
RUN_ID=bySupersite2
INDFILE=results/admixtools/input/indfile_$FILE_ID.txt
POPFILE=results/admixtools/input/popfile_dstat_hzproj1.$RUN_ID.txt
ID_COLUMN=ID.short
GROUPBY=supersite2
scripts/admixtools/01_admixtools_runner.sh $FILE_ID $RUN_ID $VCF_DIR $PLINK_DIR \
    $ATOOLS_DIR $INDFILE $POPFILE $VCF2PLINK $CREATE_INDFILE $SUBSET_INDFILE \
    $ATOOLS_MODE $INDS_METADATA $ID_COLUMN $GROUPBY

## F4-ratio-mode:
ATOOLS_MODE="F4RATIO"
RUN_ID=bySupersite2
VCF2PLINK=TRUE
INDFILE=analyses/admixtools/input/indfile_$FILE_ID.txt
POPFILE=analyses/admixtools/input/popfile_f4ratio_hzproj1.$RUN_ID.txt
scripts/admixtools/01_admixtools_runner.sh $FILE_ID $RUN_ID $VCF_DIR $PLINK_DIR \
    $ATOOLS_DIR $INDFILE $POPFILE $VCF2PLINK $CREATE_INDFILE $SUBSET_INDFILE \
    $ATOOLS_MODE $INDS_METADATA $ID_COLUMN $GROUPBY

## TO DO: PUT THIS IN FILTERING PIPELINE
# VCF=results/geno/vcf/gatk/final/r03.wOutgroups.mac3.FS6.vcf.gz
# bcftools query -l $VCF > results/geno/filtering/r03.wOutgroups.mac3.FS6_indlist.txt
# VCF=results/geno/vcf/gatk/final/r03.wOutgroups.mac1.FS6.vcf.gz
# bcftools query -l $VCF > results/geno/filtering/r03.wOutgroups.mac1.FS6_indlist.txt
