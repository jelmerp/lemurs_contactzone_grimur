# MAIN DATASET -----------------------------------------------------------------
## Prep input and run at once:
run_id=r03.all
vcf=results/gatk/r03.all.mac1.FS6.vcf
newhyb_in=results/newhybrids/input/r03.all.txt
newhyb_outdir=results/newhybrids/output
mem=12
sbatch --mem="$mem"G scripts/newhybrids/newhybrids_runner.sh "$run_id" "$vcf" "$newhyb_in" "$newhyb_outdir" "$mem"


# MAIN DATASET - WITH PREDEFINED POPS ------------------------------------------
#mgri - allo: mgri075 mgri076 mgri077 mgri078 mgri079 mgri080 mgri082
#mmur - allo: mmur037 mmur038 mmur039 mmur040 mmur041 mmur042 mmur043 mmur044

## Edit input file
newhyb_in=results/newhybrids/input/r03.all.txt
cat "$newhyb_in" |
    sed 's/mgri075/mgri075 z1/' | sed 's/mgri076/mgri076 z1/' | sed 's/mgri077/mgri077 z1/' |
    sed 's/mgri078/mgri078 z1/' | sed 's/mgri079/mgri079 z1/' | sed 's/mgri080/mgri080 z1/' | sed 's/mgri082/mgri082 z1/' |
    sed 's/mmur037/mmur037 z0/' | sed 's/mmur038/mmur038 z0/' | sed 's/mmur039/mmur039 z0/' | sed 's/mmur040/mmur040 z0/' |
    sed 's/mmur041/mmur041 z0/' | sed 's/mmur042/mmur042 z0/' | sed 's/mmur043/mmur043 z0/' | sed 's/mmur044/mmur044 z0/' >$newhyb_in.predef

## Run Newhybrids
run_id=r03.all.predef
outdir=results/newhybrids/output/"$run_id"
burnin=10000
nsweeps=50000
sbatch --mem=12G -o scripts/newhybrids/newhybrids.sh "$newhyb_in".predef "$outdir" "$burnin" "$nsweeps"


# DATASET WITH POOR-QUALITY PUTATIVE HYBRIDS -----------------------------------
vcf=results/gatk/r03.keepHybs.mac1.FS6.vcf
newhyb_in=results/newhybrids/input/r03.keepHybs.txt
run_id=r03.keepHybs
mem=12
sbatch --mem="$mem"G scripts/newhybrids/newhybrids_runner.sh "$run_id" "$vcf" "$newhyb_in" "$mem"
