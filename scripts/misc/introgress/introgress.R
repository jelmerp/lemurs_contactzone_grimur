#install.packages('introgress')

# prepare.data() estimates the counts of alleles inherited from each of two parental population
# at each locus for admixed individuals

# The function est.h() uses the list output returned from prepare.data() and the raw genotype data for
# the admixed population(s) to compute maximum likelihood estimates of hybrid index.

# The list output returned by prepare.data and the data.frame output returned by est.h are used to estimate
# genomic clines, which is accomplished with the function genomic.clines().

################################################################################
##### SET-UP  #####
################################################################################
library(tidyverse)
library(introgress)
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')




