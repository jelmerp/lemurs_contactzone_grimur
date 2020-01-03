## Get seqs from sftp to laptop:
sftp poelstra_5372@dnaseq2.igsp.duke.edu # password: pKtcvprwG9aE
#ls -l poelstra_5372/POELSTRA_5372_190116A1 # First, failed, run. Report at http://seqweb.gcb.duke.edu/19/01/ix495rqtw1v6rld_5372_190116A1.html
#ls -l poelstra_5372/POELSTRA_5372_190125B1 # Second, redone, run. Report at http://seqweb.gcb.duke.edu/19/01/r07rxiyck9t64uv_5372_190125B1.html
get poelstra_5372/POELSTRA_5372_190116A1/* /home/jelmer/Dropbox/sc_lemurs/radseq/seqdata/fastq/POELSTRA_5372_190116A1/
get poelstra_5372/POELSTRA_5372_190125B1/* /home/jelmer/Dropbox/sc_lemurs/radseq/seqdata/fastq/POELSTRA_5372_190125B1/

## Get seqs from laptop to cluster:
rsync -avr --progress /home/jelmer/Dropbox/sc_lemurs/radseq/seqdata/fastq/POELSTRA_5372_190116A1/* jwp37@dcc-slogin-02.oit.duke.edu:/work/jwp37/radseq/seqdata/fastq/r03/POELSTRA_5372_190116A1/raw/
rsync -avr --progress --partial-dir=.rsync-partial /home/jelmer/Dropbox/sc_lemurs/radseq/seqdata/fastq/POELSTRA_5372_190125B1/* jwp37@dcc-slogin-02.oit.duke.edu:/work/jwp37/radseq/seqdata/fastq/r03/POELSTRA_5372_190125B1/raw/

## Make a copy in /datacommons::
rsync -avr --progress --partial-dir=.rsync-partial /work/jwp37/radseq/seqdata/fastq/r03/POELSTRA_5372_190116A1/raw/* /datacommons/yoderlab/data/radseq/r03_20190129_hybridzone/POELSTRA_5372_190116A1/
rsync -avr --progress --partial-dir=.rsync-partial /work/jwp37/radseq/seqdata/fastq/r03/POELSTRA_5372_190125B1/raw/* /datacommons/yoderlab/data/radseq/r03_20190129_hybridzone/POELSTRA_5372_190125B1/


## From GCB email:
# Your run performed within specifications. This is the data from the re-run we did. You can combined the data from both lane.
# See your report for the first run at: http://seqweb.gcb.duke.edu/19/01/ix495rqtw1v6rld_5372_190116A1.html
# See your report for the second run at: http://seqweb.gcb.duke.edu/19/01/r07rxiyck9t64uv_5372_190125B1.html