## TO DO:
# https://cran.r-project.org/web/packages/PopGenome/vignettes/Whole_genome_analyses_using_VCF_files.pdf
# Add GFF file - compute by site class
# CLR sweep test
# Missing sites correction
# SFS
# biallelic.structure Shared and fixed polymorphisms

# Slots Description Module
# 1 nucleotide.diversity Nucleotide diversity FST
# 2 haplotype.diversity Haplotype diversity FST
# 3 haplotype.counts Haplotype distribution FST
# 4 minor.allele.freqs Minor allele frequencies Detail
# 5 linkage.disequilibrium Linkage disequilibrium Linkage
# 6 biallelic.structure Shared and fixed polymorphisms Detail

## LD by window?
#> GENOME.class.slide <- linkage.stats(GENOME.class.slide)
#> get.linkage(GENOME.class.slide)[[1]]

## READ FASTAS
# https://rdrr.io/cran/PopGenome/man/read.big.fasta.html
# GENOME.class <- readData("FASTA_DIR") # Reads multiple fasta files in dir
# But does not seem to understand ambiguity codes??

################################################################################
##### SET-UP #####
################################################################################
library(tidyverse)
library(PopGenome)
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')

## Variables:
fileID_in <- 'r03.all.mac1.FS6'
scaffold <- 'NC_033660.1'
do_allpops <- FALSE
winsize <- 10
stepsize <- 5
wintype <- 1
pop1 <- 'griseorufus'
pop2 <- 'murinus'
pop3 <- NA
pops <- c(pop1, pop2, pop3) %>% .[!is.na(.)]
triplet <- c(pop1, pop2, pop3)
outgroup <- 'XXX'

## Process vars:
triplet_name <- paste0(triplet, collapse = '-')
fileID_out <- paste0(fileID_in, '_win', winsize, '_step', stepsize,
                     '_', scaffold, '_', triplet_name)

## Files:
infile_lookup <- '/home/jelmer/Dropbox/sc_lemurs/radseq/metadata/lookup_IDshort.txt'
infile_IDs <- 'metadata/r03/sampleIDsShort_r03.txt'

infile_vcf <- paste0('seqdata/vcf/', fileID_in, '.vcf.gz') # CHANGE TO VCF-DIR
infile_vcf.tabix <- paste0(infile_vcf, '.tbi')
infile_scafs <- '/home/jelmer/Dropbox/sc_lemurs/other/seqdata_misc/reference/mmur/scaffolds_withLength.txt'

outfile_allstats <- paste0('analyses/gscan/popgenome/output/', fileID_out, '.txt')

## Individuals:
lookup <- read.delim(infile_lookup, as.is = TRUE)
inds_all <- readLines(infile_IDs)

hybs_mgri <- c('mhyb001', 'mhyb003', 'mhyb004', 'mhyb005', 'mhyb006',
               'mhyb011', 'mhyb015', 'mhyb016')
inds_pop1 <- c(inds_all[grepl('mgri', inds_all)],
               inds_all[inds_all %in% hybs_mgri])
hybs_mmur <- c('mhyb002', 'mhyb007', 'mhyb008', 'mhyb009',
               'mhyb010', 'mhyb012', 'mhyb013', 'mhyb014')
inds_pop2 <- c(inds_all[grepl('mmur', inds_all)],
               inds_all[inds_all %in% hybs_mmur])


## Scaffolds:
scafs_all <- read.table(infile_scafs,
                        colClasses = c('character', 'integer'), header = TRUE)
scaf_length <- scafs_all$scaffold.length[match(scaffold, scafs_all$scaffold.name)]

## Report:
cat('#### File ID:', fileID_in, '\n')
cat('#### Scaffold:', scaffold, '\n')
cat('#### Triplet:', triplet, '\n')
cat('#### Do all pops:', do_allpops, '\n')
cat('#### Window size:', winsize / 1000, 'k \n')
cat('#### Step size:', stepsize / 1000, 'k \n')


################################################################################
##### FUNCTIONS #####
################################################################################
## Change population names from "pop1" etc to actual names, in a df with stats:
namepops <- function(df, stat) {
  for(i in 1:length(pops)) {
    tofind <- paste0('pop ', i, '|pop', i)
    colnames(df) <- gsub(tofind, pops[i], colnames(df))
  }
  colnames(df) <- gsub('/', '.', colnames(df))
  colnames(df) <- paste0(stat, '_', colnames(df))

  return(df)
}


################################################################################
##### READ AND PROCESS VCF #####
################################################################################
## Read vcf:
if(!file.exists(infile_vcf)) cat('#### VCF FILE NOT FOUND\n\n\n\n\n')
if(!file.exists(infile_vcf.tabix)) system(paste('tabix', infile_vcf))

# vcf <- readVCF(infile_vcf,
#                numcols = 10000,
#                include.unknown = TRUE,
#                tid = scaffold,
#                frompos = 1,
#                topos = scaf_length)

# VCF_split_into_scaffolds(infile_vcf, output.folder = 'seqdata/vcf/by_scaffold')

# cd /home/jelmer/Dropbox/sc_lemurs/hybridzone/seqdata/vcf
# bcftools view r03.all.mac1.FS6.vcf.gz NC_033660.1 > by_scaffold/NC_033660.1.vzf.gz
# bcftools view r03.all.mac1.FS6.vcf.gz NC_033661.1 > by_scaffold/NC_033661.1.vzf.gz

#?concatenate.regions()

vcf.dir <- paste0('seqdata/vcf/by_scaffold/')
vcf <- readData(vcf.dir, format = 'VCF',
                include.unknown = TRUE) # TRUE to include sites with Ns in calculations

cat('\n#### Summary of data:\n')
get.sum.data(vcf)
#str(vcf)
#vcf@region.data@biallelic.sites[[1]][1:10]
#vcf@region.data@sites.with.unknowns
#vcf@region.data@biallelic.matrix[[1]]

## Populations:
inds_vcf <- unlist(vcf@region.data@populations2)
inds_vcf <- inds_vcf[-grep('\\.2', inds_vcf)]

inds_pop1 <- inds_pop1[inds_pop1 %in% inds_vcf]
inds_pop2 <- inds_pop2[inds_pop2 %in% inds_vcf]
#if(exists(inds_pop3)) inds_pop3 <- inds_pop3[inds_pop3 %in% inds_vcf]

poplist <- list('griseorufus' = inds_pop1, 'murinus' = inds_pop2)
#pops.list <- list(inds_pop1, inds_pop2, inds_pop3)

inds_outgroup <- NA

## Set the populations & outgroup:
#vcf <- set.outgroup(vcf, inds_outgroup, diploid = TRUE)
vcf <- set.populations(vcf, poplist, diploid = TRUE)

cat("#### Populations:\n")
print(vcf@populations)

## Transform the data into windows:
slide <- sliding.window.transform(vcf, winsize, stepsize, type = wintype)
# type=1 - only SNPs; type=2 - nucleotides

## Get the genomic positions for each window:
cat('Getting window positions:\n')
win.pos <- sapply(slide@region.names, function(x){
  split <- strsplit(x, " ")[[1]][c(1, 3)]
  val <- mean(as.numeric(split))
  return(as.numeric(val))
})
names(win.pos) <- NULL
cat('\n #### Window positions:\n')
print(win.pos)


################################################################################
##### BASIC STATS #####
################################################################################
nwindows <- length(slide@region.names)
smr <- get.sum.data(vcf)
nsites <- sum(smr[, 1]) #smr[[1]] for single scaffold
nsites_biallelic <- sum(smr[, 2]) #smr[[2]] for single scaffold
cat('#### Nr of windows:', nwindows, '   Nr of sites:', nsites,
    '    Nr of biallelic sites:', nsites_biallelic, '\n')

vcf.pops <- substr(sapply(slide@populations, '[', 1), 1, 16)
cat('\n #### First ind for each population:\n'); print(vcf.pops)


################################################################################
##### NEUTRALITY STATS #####
################################################################################
## Overall:
vcf <- neutrality.stats(vcf)
neutrality.stats <- get.neutrality(vcf)

neutrality.stats[[1]] # For first population
neutrality.stats[[2]]

vcf@Tajima.D

slide <- neutrality.stats(slide)
slide.neut <- get.neutrality(slide)
#slide.neut[[1]]


################################################################################
##### DIFFERENTIATION #####
################################################################################
## Fst, overall:
vcf <- F_ST.stats(vcf, mode = 'nucleotide')
vcf@nuc.F_ST.pairwise
vcf@nuc.F_ST.vs.all

## Fst, by window:
slide <- F_ST.stats(slide, mode = 'nucleotide')
#slide <- F_ST.stats(slide, mode = 'nucleotide', new.populations = allpops)
fst.pr <- namepops(t(slide@nuc.F_ST.pairwise), stat = 'fst')
fst.sl <- namepops(slide@nuc.F_ST.vs.all, stat = 'fst')

fst.pr.mean <- round(apply(fst.pr, 2, mean, na.rm = TRUE), 4)
cat('mean pairwise fst:\n'); print(fst.pr.mean)

#if(do_allpops == TRUE) {
#   print(fst.pr.mean[mam.columns])
#   print(fst.pr.mean[gui.columns])
# }

fst.sl.mean <- round(apply(fst.sl, 2, mean, na.rm = TRUE), 4)
cat('mean single pop fst:\n'); print(fst.sl.mean)

## Dxy, by window:
slide <- diversity.stats.between(slide)
dxy <- namepops((slide@nuc.diversity.between / winsize), stat = 'dxy')
dxy.mean <- round(apply(dxy, 2, mean, na.rm = TRUE), 4)
cat('mean dxy:\n'); print(dxy.mean)

#if(do_allpops == TRUE) {
#  print(dxy.mean[mam.columns])
#  print(dxy.mean[gui.columns])
#}

## Fixed and shared SNPs:
vcf <- detail.stats(vcf, biallelic.structure = TRUE)
get.detail(vcf, biallelic.structure = TRUE)
# CHECK! AND SUM TO GET NR OF FIXED DIFFS...
# Can also do per-site FST using detail.stats(vcf, site.FST = TRUE)
# ?`detail.stats,GENOME-method`

################################################################################
#### DIVERSITY ####
################################################################################
## Overall:
vcf <- diversity.stats(vcf)
div <- get.diversity(vcf)
vcf@nuc.diversity.within / scaf_length

## By window:
slide <- diversity.stats(slide)
nucdiv <- namepops((slide@nuc.diversity.within / winsize) * 100, stat = 'nucdiv')
nucdiv.mean <- round(apply(nucdiv, 2, mean, na.rm = TRUE), 4)
cat('mean nucleotide diversity:\n'); print(nucdiv.mean)

slide <- diversity.stats(slide, pi = TRUE)
pi <- namepops((slide@Pi / winsize) * 100, stat = 'pi')
pi.mean <- round(apply(pi, 2, mean, na.rm = TRUE), 4)
cat('mean nucleotide diversity pi:\n'); print(pi.mean)


################################################################################
#### SFS ####
################################################################################
g_con <- concatenate.regions(vcf)
g_con <- detail.stats(g_con, site.spectrum = TRUE, site.FST = TRUE)
sfs_results <- get.detail(g_con)

MAFs <- g_con@region.stats@minor.allele.freqs[[1]]

sfs_pop1 <- data.frame(table(MAFs['pop 1', ])) %>%
  transmute(MAF = as.numeric(as.character(Var1)),
            freq = as.integer(as.character(Freq)))
#head(sfs_pop1)

sfs_pop2 <- data.frame(table(MAFs['pop 2', ])) %>%
  transmute(MAF = as.numeric(as.character(Var1)),
            freq = as.integer(as.character(Freq)))

sfs <- merge(sfs_pop1, sfs_pop2, by = 'MAF',
             all.x = TRUE, all.y = TRUE)

# ggplot(sfs_pop1, aes(x = MAF, y = freq, width = 0.02)) +
#   #geom_bar(stat = 'identity', position = 'identity') +
#   geom_density(stat = 'identity', position = 'identity') +
#   scale_x_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   theme_bw()


################################################################################
#### INTROGRESSION STATS: D, df, BDF ####
################################################################################

## Overall:
vcf <- introgression.stats(vcf, do.D = TRUE)
vcf@D
vcf@f
vcf <- introgression.stats(vcf, do.BDF = TRUE)
vcf@BDF

## By window:
slide <- introgression.stats(slide, do.BDF = TRUE)
BDF <- slide@BDF; head(BDF)

slide <- introgression.stats(slide, do.D = TRUE)
D <- slide@D; head(D)
f <- slide@f; head(f)
f[f > 1] <- 1
f[f < -1] <- -1

mean.bdf <- mean(BDF, na.rm = TRUE)
mean.D <- mean(D, na.rm = TRUE)
mean.f <- mean(f, na.rm = TRUE)
cat('Mean BDF:', mean.bdf, '\n')
cat('Mean D:', mean.D, '\n')
cat('Mean Fd:', mean.f, '\n')

#BDF_bayes <- slide@BDF_bayes # Doesnt exist
# B-dbf = Pr(D|M1) / Pr(D|M2) # M1 = 2 & 3 introgression; M2 = 1 & 3 introgression

## Bayescan
if(perform.bayescan == TRUE) {
  bayescan.input <- getBayes(vcf, snps = TRUE)
  bayescan.output <- BayeScanR(bayescan.input,
                               nb.pilot = 10, pilot.runtime = 2500,
                               main.runtime = 100000, discard = 50000)
}


################################################################################
#### WRITE OUTPUT #####
################################################################################
allstats <- cbind(pi, fst.sl, fst.pr, dxy)

# if(do_allpops == FALSE) {
#   introgstats <- cbind(BDF, f, D)
#   colnames(introgstats) <- c('BDF', 'fd', 'D')
#   allstats <- cbind(introgstats, allstats)
# }

write.table(allstats, outfile_allstats,
            quote = FALSE, row.names = FALSE, sep = '\t')
