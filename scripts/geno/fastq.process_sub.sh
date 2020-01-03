cd /datacommons/yoderlab/users/jelmer/hybridzone

## Combine the "failed" and 2nd run:
R1_A=/work/jwp37/radseq/seqdata/fastq/r03/POELSTRA_5372_190116A1/raw/lemurRadseqHybridZone_r03_S25_L002_R1_001.fastq
R1_B=/work/jwp37/radseq/seqdata/fastq/r03/POELSTRA_5372_190125B1/raw/lemurRadseqHybridZone_r03_S25_L002_R1_001.fastq
R1_COMBINED=/work/jwp37/radseq/seqdata/fastq/r03/raw/murgriHZ.r03_S25_L002_R1_001.fastq
cat $R1_A $R1_B > $R1_COMBINED

R2_A=/work/jwp37/radseq/seqdata/fastq/r03/POELSTRA_5372_190116A1/raw/lemurRadseqHybridZone_r03_S25_L002_R2_001.fastq
R2_B=/work/jwp37/radseq/seqdata/fastq/r03/POELSTRA_5372_190125B1/raw/lemurRadseqHybridZone_r03_S25_L002_R2_001.fastq
R2_COMBINED=/work/jwp37/radseq/seqdata/fastq/r03/raw/murgriHZ.r03_S25_L002_R2_001.fastq
cat $R2_A $R2_B > $R2_COMBINED

## Run fastq.process scripts::
LIBRARY_ID="murgriHZ.r03"
BARCODE_LIST_FILE=metadata/r03/barcodes_r03.txt
BARCODE_TABLE_FILE=metadata/r03/barcodes_df_r03.txt
DIR_RAW=/work/jwp37/radseq/seqdata/fastq/r03/raw
DIR_FLIPPED=/work/jwp37/radseq/seqdata/fastq/r03/raw_flipped/
DIR_DEMULT=/work/jwp37/radseq/seqdata/fastq/r03/raw_demultiplexed/
DIR_FINAL=/work/jwp37/radseq/seqdata/fastq/r03/processed/
DIR_STATS=analyses/qc/fastq/
STEPS_TO_SKIP="-QFM" #"-QFM" # -Q: Skip QC / -F: Skip flipping / -M: Skip demultiplexing / D: Skip trimming and deduplication
PIP_SCRIPT=/datacommons/yoderlab/users/jelmer/radseq/scripts/fastq_process/fastq.process_pip.sh
sbatch -p yoderlab,common,scavenger -o slurm.fastq.process_pip.$LIBRARY_ID \
	$PIP_SCRIPT $LIBRARY_ID $DIR_RAW $DIR_FLIPPED $DIR_DEMULT $DIR_FINAL $DIR_STATS $BARCODE_LIST_FILE $BARCODE_TABLE_FILE $STEPS_TO_SKIP

## Process QC stats:
grep "Input Read Pairs" $DIR_STATS/byInd/*r03*trimStats* > $DIR_STATS/trimStats.sum.r03.txt
grep "pairs of reads input" $DIR_STATS/byInd/*r03*dedupStats* > $DIR_STATS/dedupStats.sum.r03.txt


################################################################################
################################################################################
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/hybridzone/scripts/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/scripts/
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/radseq/scripts/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/radseq/scripts/
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/hybridzone/metadata/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/metadata/

# rsync -av jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/analyses/qc/fastq/* /home/jelmer/Dropbox/sc_lemurs/hybridzone/analyses/qc/fastq/

#grep -l -Z 'dedup.trim.sh: Done with script dedup.trim.sh' slurm.dedup* | xargs -0 -I{} mv {} analyses/fastq_process/logfiles/