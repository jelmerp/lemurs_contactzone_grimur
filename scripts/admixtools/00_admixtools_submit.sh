## General settings:
file_id=hzproj1.mac1.FS6
vcf_dir=results/gatk
plink_dir=results/admixtools/input
ATOOLS_DIR=results/admixtools/
vcf2plink=TRUE
create_indfile=TRUE
subset_indfile=TRUE
inds_metadata=metadata/hzlookup_bysample.txt
id_column=ID.short
groupby=supersite2

# RUN IN D-MODE (D-STATS) ------------------------------------------------------
atools_mode="D"
vcf2plink=TRUE
run_id=bySupersite
indfile=results/admixtools/input/indfile_$file_id.txt
popfile=results/admixtools/input/popfile_dstat_hzproj1.$run_id.txt
scripts/admixtools/01_admixtools_runner.sh $file_id $run_id $vcf_dir $plink_dir \
    $ATOOLS_DIR $indfile $popfile $vcf2plink $create_indfile $subset_indfile \
    $atools_mode $inds_metadata $id_column $groupby


# RUN IN F4-RATIO-MODE ---------------------------------------------------------
atools_mode="F4RATIO"
run_id=bySupersite
vcf2plink=TRUE
indfile=results/admixtools/input/indfile_$file_id.txt
popfile=results/admixtools/input/popfile_f4ratio_hzproj1.$run_id.txt
scripts/admixtools/01_admixtools_runner.sh $file_id $run_id $vcf_dir $plink_dir \
    $ATOOLS_DIR $indfile $popfile $vcf2plink $create_indfile $subset_indfile \
    $atools_mode $inds_metadata $id_column $groupby
