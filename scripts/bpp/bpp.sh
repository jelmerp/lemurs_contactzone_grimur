
cd /datacommons/yoderlab/users/jelmer/proj/hybridzone/analyses/bpp/input/all/runtest/

BPP=/datacommons/yoderlab/programs/bpp-4.1.4-linux-x86_64/bin/bpp
#$BPP --cfile bpp.20191105.ctl
#$BPP --cfile bpp.20200401.ctl
$BPP --cfile bpp.20200401_ed.ctl > bpp_stdout_cluster.txt

less bpp_stdout_cluster.txt

################################################################################
rsync -avr /home/jelmer/Dropbox/sc_lemurs/proj/hybridzone/analyses/bpp/input/all/runtest/ dcc:/datacommons/yoderlab/users/jelmer/proj/hybridzone/analyses/bpp/input/all/runtest/

## Modify input file:
tail -n +3 mur3gri2c.txt | sed 's/^[ \t]*//' > new.txt
sed -i '/^m/s/^/aap^/' mur3gri2c.txt

#cd /home/jelmer/Dropbox/sc_lemurs/proj/hybridzone/analyses/bpp/input/all/rundir

/home/jelmer/bin/bpp-dev/src/bpp # USE THIS ONE