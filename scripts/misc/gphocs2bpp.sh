cd /home/jelmer/Dropbox/sc_lemurs/proj/hybridzone/

INPUT_GPHOCS=analyses/gphocs/input/r03.wOutgroups.hz.mur3gri2c.gphocsInput.txt
INPUT_BPP=analyses/bpp/input/mur3gri2c.txt
sed 's/.*fa//' $INPUT_GPHOCS > $INPUT_BPP

#rsync -avr jwp37@dcc-slogin-02.oit.duke.edu:/datacommons/yoderlab/users/jelmer/proj/hybridzone/analyses/bpp/input/mur3gri2c.txt /home/jelmer/Dropbox/sc_lemurs/proj/hybridzone/analyses/bpp/input/