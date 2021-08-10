ssh anne@potto.biol.ucl.ac.uk

/d/home2/anne/Microcebus-gfox/3introgression.r3/mcmc.txt
/d/home2/anne/Microcebus-gfox/3introgression.r4/mcmc.txt

rsync -avr anne@potto.biol.ucl.ac.uk:/d/home2/anne/Microcebus-gfox/ ~/Dropbox/sc_lemurs/proj/hybridzone/analyses/bpp/ziheng_dir/




################################################################################
#### BPP TESTRUN ####
################################################################################
rsync -avr ~/Dropbox/sc_lemurs/proj/hybridzone/analyses/bpp/input/all/runtest/ dcc:/datacommons/yoderlab/users/jelmer/proj/hybridzone/analyses/bpp/input/all/runtest/
/home/jelmer/bin/bpp-dev/src/bpp # USE THIS ONE

cd /datacommons/yoderlab/users/jelmer/proj/hybridzone/analyses/bpp/input/all/runtest/
$BPP --cfile bpp.20200401_ed.ctl > bpp_stdout_cluster.txt

less bpp_stdout_cluster.txt

## Modify input file:
tail -n +3 mur3gri2c.txt | sed 's/^[ \t]*//' > new.txt
sed -i '/^m/s/^/aap^/' mur3gri2c.txt

cd /home/jelmer/Dropbox/sc_lemurs/proj/hybridzone/analyses/bpp/input/all/rundir

