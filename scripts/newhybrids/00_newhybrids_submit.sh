## Prep input and run at once:
RUN_ID=r03.all
VCF=data/vcf/gatk/final/r03.all.mac1.FS6.vcf
NEWHYBRIDS_INPUT=results/newhybrids/input/r03.all.txt
MEM=12
sbatch --mem="$MEM"G -o slurm.newhybrids.pip."$RUN_ID" \
    scripts/newhybrids/newhybrids_runner.sh "$RUN_ID" "$VCF" "$NEWHYBRIDS_INPUT" "$MEM"

# RUN WITH PREDEFINED INDS -----------------------------------------------------

#mgri - allo: mgri075 mgri076 mgri077 mgri078 mgri079 mgri080 mgri082
#mmur - allo: mmur037 mmur038 mmur039 mmur040 mmur041 mmur042 mmur043 mmur044

## Edit input file
NEWHYBRIDS_INPUT=/datacommons/yoderlab/users/jelmer/hybridzone/seqdata/newhybrids_input/r03.all.txt
cat "$NEWHYBRIDS_INPUT" |
    sed 's/mgri075/mgri075 z1/' | sed 's/mgri076/mgri076 z1/' | sed 's/mgri077/mgri077 z1/' |
    sed 's/mgri078/mgri078 z1/' | sed 's/mgri079/mgri079 z1/' | sed 's/mgri080/mgri080 z1/' | sed 's/mgri082/mgri082 z1/' |
    sed 's/mmur037/mmur037 z0/' | sed 's/mmur038/mmur038 z0/' | sed 's/mmur039/mmur039 z0/' | sed 's/mmur040/mmur040 z0/' |
    sed 's/mmur041/mmur041 z0/' | sed 's/mmur042/mmur042 z0/' | sed 's/mmur043/mmur043 z0/' | sed 's/mmur044/mmur044 z0/' >$NEWHYBRIDS_INPUT.predef

## Run Newhybrids
RUN_ID=r03.all.predef
OUTDIR=results/newhybrids/output/"$RUN_ID"
BURNIN=10000
NSWEEPS=50000
sbatch --mem=12G -o slurm.newhybrids.run."$RUN_ID" \
    scripts/newhybrids/newhybrids_run.sh "$NEWHYBRIDS_INPUT".predef "$OUTDIR" "$BURNIN" "$NSWEEPS"

# RUN WITH POOR-QUALITY HYBRID INDS KEPT ---------------------------------------

RUN_ID=r03.keepHybs
VCF=data/vcf/gatk/final/r03.keepHybs.mac1.FS6.vcf
NEWHYBRIDS_INPUT=results/newhybrids/input/r03.keepHybs.txt
MEM=12
sbatch --mem="$MEM"G -o slurm.newhybrids.pip."$RUN_ID" \
    scripts/newhybrids/newhybrids_runner.sh "$RUN_ID" "$VCF" "$NEWHYBRIDS_INPUT" "$MEM"
