# cd /datacommons/yoderlab/users/jelmer/proj/hybridzone
# mkdir -p analyses/trees/splitstree/output/ /work/jwp37/hybridzone/seqdata/nexus/ /work/jwp37/hybridzone/seqdata/fasta/

VCF_DIR=/datacommons/yoderlab/users/jelmer/proj/hybridzone/seqdata/vcf/map2mmur.gatk.joint/final # Dir with VCF files
FASTA_DIR=/work/jwp37/hybridzone/seqdata/fasta # Dir for fasta files (convert to Nexus via fasta)
NEXUS_DIR=/work/jwp37/hybridzone/seqdata/nexus/ # Dir for nexus files (Splitstree takes Nexus as input)
OUTDIR=analyses/trees/splitstree/output # Splitstree output dir
MEM=20 # Memory in GB
FILE_IDS=( hzproj1.mac1.FS6 hzproj1.mac3.FS6 ) # VCF file ID, VCF should be $VCF_DIR/$FILE_ID.vcf.gz

for FILE_ID in ${FILE_IDS[@]}
do
	echo -e "\n##### File ID: $FILE_ID"
	sbatch -p yoderlab,common,scavenger --mem=${MEM}G -o slurm.splitstree.pip.$FILE_ID \
	/datacommons/yoderlab/users/jelmer/scripts/trees/splitstree_pip.sh $FILE_ID $VCF_DIR $FASTA_DIR $NEXUS_DIR $OUTDIR $MEM
done


################################################################################
################################################################################
# rsync -avr /home/jelmer/Dropbox/sc_lemurs/scripts/* jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/scripts/
# rsync -avr jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/proj/hybridzone/analyses/trees/splitstree/output/* /home/jelmer/Dropbox/sc_lemurs/proj/hybridzone/analyses/trees/splitstree/output

################################################################################
# scp /home/jelmer/Dropbox/sc_fish/software/splitstree4_unix_4_14_4.sh jelmerp@killdevil.unc.edu:/proj/cmarlab/users/jelmer/software
# wget http://ab.inf.uni-tuebingen.de/data/software/splitstree4/download/splitstree4_unix_4_14_6.sh
# ./splitstree4_unix_4_14_6.sh # run the bash script to install

#FILE_IDS=( r03.keepHybs.mac1.FS6 r03.keepHybs.mac3.FS6 r03.keepHybs.mac1.FS7 r03.keepHybs.mac3.FS7 )
#FILE_IDS=( r03.wOutgroups.mac1.FS6 r03.wOutgroups.mac3.FS6 r03.wOutgroups.mac1.FS7 r03.wOutgroups.mac3.FS7 )