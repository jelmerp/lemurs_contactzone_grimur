## Prep input and run at once:
RUN_ID=r03.all
VCF=/work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/final/r03.all.mac1.FS6.vcf
#gunzip -c $VCF.gz > $VCF
NEWHYBRIDS_INPUT=/datacommons/yoderlab/users/jelmer/hybridzone/seqdata/newhybrids_input/r03.all.txt
MEM=12
sbatch -p yoderlab,common,scavenger --mem=${MEM}G -o slurm.newhybrids.pip.$RUN_ID \
/datacommons/yoderlab/users/jelmer/scripts/genomics/newhybrids/newhybrids_pip.sh $RUN_ID $VCF $NEWHYBRIDS_INPUT $MEM


################################################################################
## Run with predefined inds:
#mgri - allo: mgri075 mgri076 mgri077 mgri078 mgri079 mgri080 mgri082
#mmur - allo: mmur037 mmur038 mmur039 mmur040 mmur041 mmur042 mmur043 mmur044
# cat $NEWHYBRIDS_INPUT | cut -f 1,2 -d " "

## Edit input file:
NEWHYBRIDS_INPUT=/datacommons/yoderlab/users/jelmer/hybridzone/seqdata/newhybrids_input/r03.all.txt
cat $NEWHYBRIDS_INPUT | sed 's/mgri075/mgri075 z1/' | sed 's/mgri076/mgri076 z1/' | sed 's/mgri077/mgri077 z1/' | \
	sed 's/mgri078/mgri078 z1/' | sed 's/mgri079/mgri079 z1/' | sed 's/mgri080/mgri080 z1/' | sed 's/mgri082/mgri082 z1/' | \
	sed 's/mmur037/mmur037 z0/' | sed 's/mmur038/mmur038 z0/' | sed 's/mmur039/mmur039 z0/' | sed 's/mmur040/mmur040 z0/' | \
	sed 's/mmur041/mmur041 z0/' | sed 's/mmur042/mmur042 z0/' | sed 's/mmur043/mmur043 z0/' | sed 's/mmur044/mmur044 z0/' > $NEWHYBRIDS_INPUT.predef
# cat $NEWHYBRIDS_INPUT.predef | cut -f 1,2,3 -d " "

## Run Newhybrids:
SCRIPT_NEWHYBRIDS=/datacommons/yoderlab/users/jelmer/scripts/genomics/newhybrids/newhybrids_run.sh
RUN_ID=r03.all.predef.quick
OUTDIR=analyses/newhybrids/output/$RUN_ID/
BURNIN=5000
NSWEEPS=10000
sbatch -p yoderlab,common,scavenger --mem=12G -o slurm.newhybrids.run.$RUN_ID \
$SCRIPT_NEWHYBRIDS $NEWHYBRIDS_INPUT.predef $OUTDIR $BURNIN $NSWEEPS

RUN_ID=r03.all.predef
OUTDIR=analyses/newhybrids/output/$RUN_ID/
BURNIN=10000
NSWEEPS=50000
sbatch -p yoderlab,common,scavenger --mem=12G -o slurm.newhybrids.run.$RUN_ID \
$SCRIPT_NEWHYBRIDS $NEWHYBRIDS_INPUT.predef $OUTDIR $BURNIN $NSWEEPS


################################################################################
## Run with bad hybrid inds kept:
RUN_ID=r03.keepHybs
VCF=/work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/final/r03.keepHybs.mac1.FS6.vcf
#gunzip -c $VCF.gz > $VCF
NEWHYBRIDS_INPUT=/datacommons/yoderlab/users/jelmer/hybridzone/seqdata/newhybrids_input/r03.keepHybs.txt
MEM=12
sbatch -p yoderlab,common,scavenger --mem=${MEM}G -o slurm.newhybrids.pip.$RUN_ID \
/datacommons/yoderlab/users/jelmer/scripts/genomics/newhybrids/newhybrids_pip.sh $RUN_ID $VCF $NEWHYBRIDS_INPUT $MEM

RUN_ID=r03.keepHybs.quick
NEWHYBRIDS_INPUT=/datacommons/yoderlab/users/jelmer/hybridzone/seqdata/newhybrids_input/r03.keepHybs.txt
OUTDIR=analyses/newhybrids/output/$RUN_ID/
BURNIN=5000
NSWEEPS=10000
sbatch -p yoderlab,common,scavenger --mem=12G -o slurm.newhybrids.run.$RUN_ID \
/datacommons/yoderlab/users/jelmer/scripts/genomics/newhybrids/newhybrids_run.sh $NEWHYBRIDS_INPUT $OUTDIR $BURNIN $NSWEEPS



################################################################################
## TO DO: Run separately for the two sites:
#VCF=/work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/final/r03.all.mac1.FS6.vcf



################################################################################
## Conversion only:
SCRIPT_CONVERT=/datacommons/yoderlab/users/jelmer/scripts/genomics/conversion/vcf2newhybrids.sh
VCF=/work/jwp37/hybridzone/seqdata/vcf/map2mmur.gatk4.paired.joint/final/r03.all.mac1.FS6.vcf
NEWHYBRIDS_INPUT=/datacommons/yoderlab/users/jelmer/hybridzone/seqdata/newhybrids_input/r03.all.txt
SPIDFILE=/datacommons/yoderlab/users/jelmer/scripts/genomics/conversion/vcf2newhybrids.spid
MEM=4
$SCRIPT_CONVERT $VCF $NEWHYBRIDS_INPUT $SPIDFILE $MEM



################################################################################
################################################################################
# rsync -avr --no-perms /home/jelmer/Dropbox/sc_lemurs/hybridzone/scripts/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/scripts/
# rsync -avr --no-perms /home/jelmer/Dropbox/scripts/genomics/ jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/scripts/genomics/

# rsync -avr --no-perms jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/analyses/newhybrids/output/ /home/jelmer/Dropbox/sc_lemurs/hybridzone/analyses/newhybrids/output/

# /datacommons/yoderlab/programs/newhybrids/newhybrids-no-gui-linux.exe --help
# /datacommons/yoderlab/programs/newhybrids/newhybrids-no-gui-linux.exe --help-full
# INFILE=/datacommons/yoderlab/programs/newhybrids/test_data/TestDat.txt
# rsync -avr --no-perms jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/hybridzone/seqdata/newhybrids_input/template_VCF_NEWHYBRIDS.spid /home/jelmer/Dropbox/sc_lemurs/hybridzone/
