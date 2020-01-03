## TO DO:
# Add GFF file - compute by site class
# CLR sweep test?
# Missing sites correction?

# Haplotypes? haplotype.diversity
# LD? linkage.disequilibrium

## LD by window?
#> GENOME.class.slide <- linkage.stats(GENOME.class.slide)
#> get.linkage(GENOME.class.slide)[[1]]

# Switch populations: F_ST.stats(.., new.populations = allpops)

## !! READ FASTAS
# https://rdrr.io/cran/PopGenome/man/read.big.fasta.html
# GENOME.class <- readData("FASTA_DIR") # Reads multiple fasta files in dir
# But does not seem to accept ambiguity codes??

## !! PER-SCAFFOLD - COMPARE SEX CHROM


################################################################################
##### SET-UP #####
################################################################################
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')
library(tidyverse)
library(PopGenome)
source('/home/jelmer/Dropbox/sc_lemurs/scripts/gscan/popgenome_fun.R')

## Variables:
fileID_in <- 'r03.all.mac1.FS6'
indir_vcf <- 'seqdata/vcf/by_scaffold/'
scaffold <- 'all'
do_bayescan <- FALSE
do_introgstats <- FALSE

wintype <- 1 # 1 = by SNP; 2 = by coordinates (regardless of nr of SNPs)
winsize <- 10
stepsize <- 5

## Populations:
pop1 <- 'gri'
pop2 <- 'mur'
pop3 <- NA
pop4 <- NA
pops <- c(pop1, pop2, pop3, pop4) %>% .[!is.na(.)]
outgroup <- 'none'
npop <- length(pops)

if(do_introgstats == TRUE) {
 triplet <- c(pop1, pop2, pop3)
 triplet_name <- paste0(triplet, collapse = '-')
}

## Input files:
infile_lookup <- '/home/jelmer/Dropbox/sc_lemurs/radseq/metadata/lookup_IDshort.txt'
infile_IDs <- 'metadata/r03/sampleIDsShort_r03.txt'
infile_scafs <- '/home/jelmer/Dropbox/sc_lemurs/other/seqdata_misc/reference/mmur/scaffolds_withLength.txt'
if(scaffold != 'all') {
  infile_vcf <- paste0('seqdata/vcf/', fileID_in, '.vcf.gz') # CHANGE TO VCF-DIR
  infile_vcf.tabix <- paste0(infile_vcf, '.tbi')
}

## Output files:
outdir <- 'analyses/gscan/popgenome/output/'
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
outfile_base1 <- paste0(outdir, '/', fileID_in, '_win', winsize, '_step', stepsize)
outfile_base2 <- paste0(outdir, '/', fileID_in)

outfile_between_con <- paste0(outfile_base2, '_meanBtwPopStats')
outfile_within_con <- paste0(outfile_base2, '_meanWthPopStats')
outfile_sitepat_smr <- paste0(outfile_base2, '_sitepatSmr.txt')
outfile_sitestats <- paste0(outfile_base2, '_sitestats.txt')
outfile_scafstats <- paste0(outfile_base2, '_scafstats.txt')
outfile_winstats <- paste0(outfile_base1, '_winstats.txt')

outfile_sfs <- paste0(outfile_base2, '_sfs.txt')
outfile_introgstats <- paste0(outfile_base2, '_introgstats.txt')

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
scaf_all <- read.table(infile_scafs, header = TRUE,
                       colClasses = c('character', 'integer'))
scaf_len <- scaf_all$scaffold.length[match(scaffold, scaf_all$scaffold.name)]

## Report:
cat('#### File ID:', fileID_in, '\n')
cat('#### Scaffold:', scaffold, '\n')

cat('#### Populations:', pops, '\n')
cat('#### Outgroup:', outgroup, '\n')
cat('#### Nr of populations:', npop, '\n')
if(do_introgstats == TRUE) cat('#### Triplet:', triplet, '\n')

cat('#### Window type:', wintype, '\n')
cat('#### Window size:', winsize, '\n')
cat('#### Step size:', stepsize, '\n')


################################################################################
##### READ VCF #####
################################################################################
# VCF_split_into_scaffolds(infile_vcf, output.folder = 'seqdata/vcf/by_scaffold') # slow
# bcftools view r03.all.mac1.FS6.vcf.gz NC_033660.1 > by_scaffold/NC_033660.1.vzf.gz

## Read vcf:
if(scaffold != 'all') {
  cat('#### Loading single-scaffold from VCF...\n')
  cat('#### Scaffold:', scaffold, '\n')
  cat('#### VCF file:', infile_vcf, '\n')
  if(!file.exists(infile_vcf)) cat('#### VCF FILE NOT FOUND\n\n\n\n\n')
  if(!file.exists(infile_vcf.tabix)) system(paste('tabix', infile_vcf))
  g_scaf <- readVCF(infile_vcf,
                    numcols = 10000,
                    include.unknown = TRUE,
                    tid = scaffold,
                    frompos = 1,
                    topos = scaf_len)
} else {
  cat('#### Loading multiple VCF files from dir...\n')
  cat('#### Dir with VCF files:', indir_vcf, '\n')
  g_scaf <- readData(indir_vcf,
                     format = 'VCF',
                     include.unknown = TRUE) # include.unknown=TRUE to include sites with Ns in calculations
}


################################################################################
##### PREP VCF AND BASIC STATS #####
################################################################################
## Summarize and report:
g_summary <- get.sum.data(g_scaf) %>%
  as.data.frame() %>%
  rename(nsite = n.sites,
         nvar = n.biallelic.sites,
         TiTv = trans.transv.ratio) %>%
  rownames_to_column('scaffold') %>%
  select(scaffold, nsite, nvar, TiTv)

cat('\n#### Summary of loaded sequence data:\n')
head(g_summary, n = 3)
cat('\n#### Number of scaffolds:', nrow(g_summary), '\n')
cat('#### Total nr of biallelic sites:', sum(g_summary$nvar), '\n')

## Set populations:
inds_g <- unlist(g_scaf@region.data@populations2)
inds_g <- inds_g[-grep('\\.2', inds_g)]
inds_pop1 <- inds_pop1[inds_pop1 %in% inds_g]
inds_pop2 <- inds_pop2[inds_pop2 %in% inds_g]
if(exists("inds_pop3")) inds_pop3 <- inds_pop3[inds_pop3 %in% inds_g]

poplist <- list(inds_pop1, inds_pop2) # TO DO: MODIFY TO ALLOW FOR ARBITRARY NR OF POPS
names(poplist) <- pops
g_scaf <- set.populations(g_scaf, poplist, diploid = TRUE)

cat("\n#### Populations and nr of inds:\n")
inds_per_pop <- unlist(lapply(g_scaf@populations, length)) / 2
print(inds_per_pop)

## Set outgroup:
if(outgroup != 'none') {
  cat("\n#### Setting outgroup...\n")
  inds_outgroup <- NA
  g_scaf <- set.outgroup(g_scaf, inds_outgroup, diploid = TRUE)
}

## Scaffolds/regions:
scaffolds_g <- gsub('.vcf.gz', '', g_scaf@region.names)


################################################################################
##### SLIDING WINDOWS AND CONCATENATED DATA #####
################################################################################
## Transform the data into windows:
cat("\n#### Creating sliding windows...\n")
g_win <- sliding.window.transform(g_scaf,
                                  winsize,
                                  stepsize,
                                  type = wintype, # type=1 - only SNPs; type=2 - across all nucleotides (i.e. using VCF coords)
                                  whole.data = TRUE) # if FALSE, scan scaffolds separately

## Get the genomic positions for each window:
cat('#### Getting window positions...\n')
win_start <- sapply(g_win@region.names,
                    function(x) as.numeric(strsplit(x[1], " ")[[1]][1]))
names(win_start) <- NULL
win_end <- sapply(g_win@region.names,
                   function(x) as.numeric(strsplit(x[1], " ")[[1]][3]))
names(win_end) <- NULL

nwin <- length(g_win@region.names)
cat('\n #### Number of windows:', nwin, '\n')

## Concatenate data from all scaffolds/regions:
cat('\n#### Concatenating regions...\n')
g_con <- concatenate.regions(g_scaf)


################################################################################
##### BETWEEN-POP: DXY #####
################################################################################
cat('\n\n#### Calculating Dxy...\n')

## Concat:
g_con <- diversity.stats.between(g_con)
dxy_con <- g_con@nuc.diversity.between

dxy_con <- t(g_con@nuc.diversity.between) %>%
  as.data.frame() %>%
  rownames_to_column('pair') %>%
  rename(value = V1) %>%
  mutate(stat = 'dxy')
cat("\n#### Mean pairwise Dxy:\n")
print(dxy_con)

## By window:
g_win <- diversity.stats.between(g_win)
dxy_win <- g_win@nuc.diversity.between %>%
  namepops(., stat = 'dxy')

## By scaffold:
#g_scaf <- diversity.stats.between(g_scaf)


################################################################################
##### BETWEEN-POP: FST (OVERALL AND PER-WINDOW) #####
################################################################################
cat('\n\n#### Calculating FST...\n')
# get.F_ST(g_con) # doesnt work

## Concatenated:
g_con <- F_ST.stats(g_con, mode = 'nucleotide')

fst_con_pair <- g_con@nuc.F_ST.pairwise %>%
  as.data.frame() %>%
  rownames_to_column('pair') %>%
  rename(value = V1) %>%
  mutate(stat = 'FST')
cat("#### Mean pairwise FST:\n")
print(fst_con_pair)

fst_con_single <- g_con@nuc.F_ST.vs.all %>%
  namepops2(., stat = 'FST')
cat("#### Mean single-pop FST:\n")
print(fst_con_single)

## By window:
g_win <- F_ST.stats(g_win, mode = 'nucleotide')
fst_win_pair <- t(g_win@nuc.F_ST.pairwise) %>% namepops(., stat = 'fst')
fst_win_single <- g_win@nuc.F_ST.vs.all %>% namepops(., stat = 'fst')

## By scaffold:
# g_scaf <- F_ST.stats(g_scaf, mode = 'nucleotide')


################################################################################
##### BETWEEN-POP: PER-SITE FST #####
################################################################################
cat('\n\n#### Calculating per-site FST...\n')

g_scaf <- detail.stats(g_scaf, site.FST = TRUE)
fst_site_list <- g_scaf@region.stats@site.FST

fst_site <- as.data.frame(fst_site_list[[1]])
colnames(fst_site) <- 'Fst'
fst_site$scaffold <- scaffolds_g[1]
for(i in 2:length(fst_site_list)) {
  if(!is.null(fst_site_list[[i]])) {
    newrows <- as.data.frame(fst_site_list[[i]])
    colnames(newrows) <- 'Fst'
    newrows$scaffold <- scaffolds_g[i]
    fst_site <- rbind(fst_site, newrows)
  }
}

fst_site <- fst_site %>%
  rownames_to_column('site') %>%
  select(scaffold, site, contains('Fst')) %>%
  mutate(Fst = round(Fst, 4)) %>%
  mutate(Fst = ifelse(Fst < 0, 0, Fst)) # If FST < 0 --> make 0
#head(fst_site)
## TO DO: ACCOMODATE MULTIPLE POP-COMPS


################################################################################
##### BETWEEN-POP: SITE PATTERNS (FIXED, SHARED, ETC) #####
################################################################################
cat('\n\n#### Calculating site patterns...\n')
# head(t(g_scaf@region.data@biallelic.substitutions[[1]]))

# 0	population is polymorphic, the remaining individuals are polymorphic
# 1	population is polymorphic, the remaining individuals are monomorphic
# 2	population is monomorphic, the remaining individuals are polymorphic
# 3	population is monomorphic, the remaining individuals are monomorphic with the same value
# 4	population is monomorphic, the remaining individuals are monomorphic with different values

g_scaf <- detail.stats(g_scaf, biallelic.structure = TRUE)
sitepat_list <- get.detail(g_scaf, biallelic.structure = TRUE)

sitepat <- as.data.frame(t(sitepat_list[[1]]))
sitepat$scaffold <- scaffolds_g[1]
for(i in 2:length(sitepat_list)) {
  if(!is.null(sitepat_list[[i]])) {
    newrows <- as.data.frame(t(sitepat_list[[i]]))
    newrows$scaffold <- scaffolds_g[i]
    sitepat <- rbind(sitepat, newrows)
  }
}
colnames(sitepat)[1:npop] <- paste0('sipat_', pops)
sitepat$site <- rownames(sitepat)
rownames(sitepat) <- NULL
sitepat <- sitepat %>% select(scaffold, site, contains('sipat'))

## Site pattern summary:
sitepat_smr <- sitepat %>%
  select(contains('sipat')) %>%
  gather("species", "sipat", 1:2) %>%
  group_by(species, sipat) %>%
  summarize(count = n()) %>%
  spread(key = species, value = count)
cat('\n#### Site pattern summary:\n')
print(sitepat_smr)

if(npop == 2) {
  stat <- paste0('sipat', sitepat_smr$sipat)
  value <- sitepat_smr[, 2]
  pair <- c('pop1/pop2')
  sitepat_smr2 <- cbind(pair, value, stat)
  colnames(sitepat_smr2)[2] <- 'value'
}

cat('\n#### Nr of fixed differences:', sum(sitepat$sipat_gri == 4), '\n')
cat('#### Total nr of sites:', nrow(sitepat), '\n')
cat('#### Nr of shared polymorphisms:', sum(sitepat$sipat_mur == 0), '\n')
cat('#### Nr of polymorphisms unique to mur:', sum(sitepat$sipat_mur == 1), '\n')
cat('#### Nr of polymorphisms unique to gri:', sum(sitepat$sipat_gri == 1), '\n')


################################################################################
#### WITHIN-POP: DIVERSITY ####
################################################################################
cat('\n\n#### Calculating diversity...\n')

## Concatenated:
g_con <- diversity.stats(g_con, pi = TRUE)
div_con <- get.diversity(g_con)

div_smr <- do.call(rbind, lapply(1:npop, function(x) div_con[[x]]))
div_smr <- as.data.frame(t(data.frame(div_smr, row.names = NULL)))
colnames(div_smr) <- pops
div_smr <- as.data.frame(apply(div_smr, 2, round, 5))
div_smr$stat <- c('nucdiv', 'hapdiv', 'pi', 'hapFst_all', 'Fst_all')
rownames(div_smr) <- NULL

cat('\n#### Summary of diversity stats:\n')
print(div_smr)

## By window:
g_win <- diversity.stats(g_win, pi = TRUE)
div_win <- get.diversity(g_win)
#pi <- namepops((g_win@Pi / winsize), stat = 'pi')
pi_win <- get_stat_df(div_win, 'Pi', 'pi')
#head(pi_win)

## By scaffold:
#g_scaf <- diversity.stats(g_scaf, pi = TRUE)


################################################################################
##### WITHIN-POP: NEUTRALITY STATS #####
################################################################################
cat('\n\n#### Calculating neutrality stats...\n')

## Concatenated:
g_con <- neutrality.stats(g_con)
neutral_con <- get.neutrality(g_con)

neutral_smr <- do.call(rbind, lapply(1:npop, function(x) neutral_con[[x]]))
neutral_smr <- as.data.frame(t(data.frame(neutral_smr, row.names = NULL)))
colnames(neutral_smr) <- pops
neutral_smr <- as.data.frame(apply(neutral_smr, 2, round, 5))
neutral_smr$stat <- c('TajD', 'nseg', 'rozasR2', 'FuLiF', 'FuLiD',
                      'FuFs', 'FayWhuH', 'ZengE', 'StrobeckS')
rownames(neutral_smr) <- NULL

cat('\n#### Summary of neutrality stats:\n')
print(neutral_smr)

## By window:
g_win <- neutrality.stats(g_win)
neutral_win <- get.neutrality(g_win)
taj_win <- get_stat_df(neutral_win, 'Tajima.D', 'tajD')
#head(taj_win)

## Per scaffold:
#g_scaf <- neutrality.stats(g_scaf)


################################################################################
#### MAFs and SFS ####
################################################################################
cat('\n\n#### Calculating MAFs and SFS...\n')

## Concatenated:
g_con <- detail.stats(g_con, site.spectrum = TRUE) #, site.FST = TRUE)
MAF_con <- g_con@region.stats@minor.allele.freqs[[1]]

MAF <- as.data.frame(t(MAF_con)) %>%
  namepops(., stat = 'MAF')
MAF <- as.data.frame(apply(MAF, 2, round, 5))
#head(MAF)

sfs_pop1 <- data.frame(table(MAF_con['pop 1', ])) %>%
  transmute(MAF = as.numeric(as.character(Var1)),
            freq = as.integer(as.character(Freq)))
colnames(sfs_pop1)[2] <- paste0('freq.', pop1)

sfs_pop2 <- data.frame(table(MAF_con['pop 2', ])) %>%
  transmute(MAF = as.numeric(as.character(Var1)),
            freq = as.integer(as.character(Freq)))
colnames(sfs_pop2)[2] <- paste0('freq.', pop2)

sfs <- merge(sfs_pop1, sfs_pop2, by = 'MAF',
             all.x = TRUE, all.y = TRUE)
#head(sfs)


################################################################################
#### INTROGRESSION STATS: D, df, BDF ####
################################################################################
if(do_introgstats == TRUE) {
  cat('\n\n#### Calculating MAFs and SFS...\n')

  ## Concatenated:
  g_con <- introgression.stats(g_con, do.D = TRUE)
  g_con@D
  g_con@f
  g_con <- introgression.stats(g_con, do.BDF = TRUE)
  g_con@BDF

  ## By window:
  g_win <- introgression.stats(g_win, do.BDF = TRUE)
  BDF <- g_win@BDF; head(BDF)

  g_win <- introgression.stats(g_win, do.D = TRUE)
  D <- g_win@D; head(D)
  f <- g_win@f; head(f)
  f[f > 1] <- 1
  f[f < -1] <- -1

  # BDF_bayes <- g_win@BDF_bayes # Doesnt exist
  # B-dbf = Pr(D|M1) / Pr(D|M2) # M1 = 2 & 3 introgression; M2 = 1 & 3 introgression
}


################################################################################
#### BAYESCAN #####
################################################################################
if(do_bayescan == TRUE) {
  cat('\n\n#### Running Bayescan...\n')
  bayescan.input <- getBayes(g_con, snps = TRUE)
  bayescan.output <- BayeScanR(bayescan.input,
                               nb.pilot = 10, pilot.runtime = 2500,
                               main.runtime = 100000, discard = 50000)
}


################################################################################
#### CREATE FINAL DFs #####
################################################################################
## Overall single-pop stats (averaged across all data):
within_con <- rbind(div_smr, fst_con_single, neutral_smr)
within_con

## Overall between-pop stats (averaged across all data):
between_con <- rbind(fst_con_pair, dxy_con, sitepat_smr2) %>%
  select(stat, pair, value) %>%
  mutate(value = round(value, 4))
between_con$pair <- gsub('pop1', pop1, between_con$pair)
between_con$pair <- gsub('pop2', pop2, between_con$pair)
between_con

## Per-site statistics:
sitepat_and_MAF <- cbind(sitepat, MAF)
sitestats <- merge(fst_site, sitepat_and_MAF, by = c('scaffold', 'site'))
head(sitestats)

## Per scaffold stats:
scafstats <- g_summary ## ADD ADDITIONAL STATS
head(scafstats)

## Windowstats:
scaffold <- sitestats$scaffold[win_start]
scaf_win_end <- sitestats$scaffold[win_end]
win_remove <- which(scaffold != scaf_win_end) # Windows that cross scaffold boundaries will be removed
start <- sitestats$site[win_start]
end <- sitestats$site[win_end]

winstats <- cbind(scaffold, start, end,
                  pi_win, fst_win_pair, fst_win_single, taj_win)
if(do_introgstats == TRUE) winstats <- cbind(winstats, introg_win)
winstats <- winstats %>% slice(-win_remove)
head(winstats)
dim(winstats)


################################################################################
#### WRITE OUTPUT #####
################################################################################
if(exists("within_con")) write_table(within_con, outfile_within_con)
if(exists("between_con")) write_table(between_con, outfile_between_con)
if(exists("sitepat_smr")) write_table(sitepat_smr, outfile_sitepat_smr)

if(exists("scafstats")) write_table(scafstats, outfile_scafstats)
if(exists("sitestats")) write_table(sitestats, outfile_sitestats)
if(exists("winstats")) write_table(winstats, outfile_winstats)

if(exists("sfs")) write_table(sfs, outfile_sfs)
if(exists("introgstats")) write_table(introgstats, outfile_introgstats)