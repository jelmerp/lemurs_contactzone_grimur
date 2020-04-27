cd /datacommons/yoderlab/data/radseq/fastq/

################################################################################
#### PREP ####
################################################################################
## Combine the "failed" and 2nd run:
R1_A=r03_20190129/original/by_run/POELSTRA_5372_190116A1/raw/lemurRadseqHybridZone_r03_S25_L002_R1_001.fastq
R1_B=r03_20190129/original/by_run//POELSTRA_5372_190125B1/raw/lemurRadseqHybridZone_r03_S25_L002_R1_001.fastq
R1_COMBINED=r03_20190129/original/murgriHZ.r03_S25_L002_R1_001.fastq
cat $R1_A $R1_B > $R1_COMBINED

R2_A=r03_20190129/original/by_run/POELSTRA_5372_190116A1/raw/lemurRadseqHybridZone_r03_S25_L002_R2_001.fastq
R2_B=r03_20190129/original/by_run/POELSTRA_5372_190125B1/raw/lemurRadseqHybridZone_r03_S25_L002_R2_001.fastq
R2_COMBINED=r03_20190129/original/murgriHZ.r03_S25_L002_R2_001.fastq
cat $R2_A $R2_B > $R2_COMBINED

################################################################################
#### RUN SCRIPT TO PROCESS FASTQS ####
################################################################################
# Q) SKIP_QC F) SKIP_FLIP M) SKIP_DEMULT D) SKIP_DEDUP T) SKIP_TRIM S) SKIP_STATS
LIBRARY_ID="murgriHZ.r03"
INDIR=r03/original/
OUTDIR=r03/processed/
BARCODES=r03/metadata/barcodes_df_r03.txt #Two-column file with barcodes and corresponding sample IDs
SKIP="" #"-QFM"
PIP_SCRIPT=scripts/genomics/geno/fastq/rad/fqprocess_pip.sh
sbatch -p yoderlab,common,scavenger -o slurm.fqprocess_pip.$LIBRARY_ID \
	$PIP_SCRIPT -i $INDIR -o $OUTDIR -l $LIBRARY_ID -b $BARCODES -s $SKIP

## Process QC stats:
grep "Input Read Pairs" slurm.dedup* > r03/processed/QC/trimstats_sumr.txt
grep "pairs of reads input" slurm.dedup* > r03/processed/QC/dedupstats_sumr.txt


################################################################################
################################################################################
rsync -av --no-perms /home/jelmer/Dropbox/scripts/genomics/geno/fastq/rad/* dcc:/datacommons/yoderlab/data/radseq/fastq/scripts/genomics/geno/fastq/rad/
rsync -av --no-perms /home/jelmer/Dropbox/scripts/genomics/qc/*fastq* dcc:/datacommons/yoderlab/data/radseq/fastq/scripts/genomics/qc/

rsync -av --no-perms dcc:/datacommons/yoderlab/data/radseq/fastq/r03/processed/QC/*sumr* /home/jelmer/Dropbox/sc_lemurs/proj/hybridzone/analyses/qc/fastq/