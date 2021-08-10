## Settings
VCF_DIR=results/geno/vcf/gatk/final       # Dir with VCF files
FASTA_DIR=results/geno/conversions/fasta  # Dir for fasta files (convert to Nexus via fasta)
NEXUS_DIR=results/geno/conversions/nexus/ # Dir for nexus files (Splitstree takes Nexus as input)
OUTDIR=results/splitstree/output          # Splitstree output dir
MEM=20                                    # Memory in GB

FILE_IDS=(hzproj1.mac1.FS6 hzproj1.mac3.FS6) # VCF file ID, VCF should be $VCF_DIR/$FILE_ID.vcf.gz

## Run Splitstree
for file_id in "${FILE_IDS[@]}"; do
    echo -e "\n## File ID: $file_id"
    sbatch --mem=${MEM}G -o slurm.splitstree.pip."$file_id" \
        scripts/01_splitstree_runner.sh "$file_id" $VCF_DIR $FASTA_DIR $NEXUS_DIR $OUTDIR $MEM
done
