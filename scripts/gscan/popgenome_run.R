## TO DO:
# Add GFF file - compute by site class
# CLR sweep test?
# Missing sites correction?
# Haplotypes? haplotype.diversity - for fasta?
# LD? linkage.disequilibrium
# Switch populations? F_ST.stats(.., new.populations = allpops)

## LD by window?
#> GENOME.class.slide <- linkage.stats(GENOME.class.slide)
#> get.linkage(GENOME.class.slide)[[1]]

## !! READ FASTAS
# https://rdrr.io/cran/PopGenome/man/read.big.fasta.html
# GENOME.class <- readData("FASTA_DIR") # Reads multiple fasta files in dir
# But does not seem to accept ambiguity codes??

################################################################################
##### SET-UP - INPUT #####
################################################################################
library(tidyverse)
library(PopGenome)

## Variables:
#fileID_in <- 'r03.all.mac1.FS6'
fileID_in <- 'hzproj1.mac1.FS6'
infile_triplets <- 'hybridzone/analyses/gscan/popgenome/input/triplets2.txt'
infile_pairs <- 'hybridzone/analyses/gscan/popgenome/input/pairs1.txt'

scaffold <- 'all'
popIDcolumn <- 'supersite2'
wintype <- 1 # 1 = by SNP; 2 = by coordinates (regardless of nr of SNPs)
winsize <- 25
stepsize <- 5
indir_vcf <- paste0('hybridzone/seqdata/vcf/by_scaffold/', fileID_in, '/')

do_pop_pairs <- TRUE
do_introgstats <- TRUE

do_betweenpop <- TRUE
do_withinpop <- TRUE

do_dxy <- TRUE
do_fst_win <- TRUE
do_fst_site <- TRUE
do_sitepat <- TRUE
do_diversity <- TRUE
do_neutral <- TRUE
do_sfs <- TRUE
do_summaries <- FALSE
do_bayescan <- FALSE

## Populations:
pop1 <- 'griC'
pop2 <- 'griW'
pop3 <- 'murC'
pop4 <- 'murE'
pop5 <- 'murW'
outgroup <- 'ruf'

if(do_pop_pairs == TRUE) {
  pop_pairs_df <- read.delim(infile_pairs, sep = ' ', header = FALSE, as.is = TRUE)
  pop_pairs_list <- split(pop_pairs_df, seq(nrow(pop_pairs_df)))
  npop_pairs <- length(pop_pairs_list)
  cat('#### Nr of pop_pairs:', npop_pairs, '\n')
}

if(do_introgstats == TRUE) {
  triplets_df <- read.delim(infile_triplets, sep = ' ', header = FALSE, as.is = TRUE)
  triplets_list <- split(triplets_df, seq(nrow(triplets_df)))
  ntriplets <- length(triplets_list)
  cat('#### Nr of triplets:', ntriplets, '\n')
}

## Basedir & scripts:
setwd('/home/jelmer/Dropbox/sc_lemurs/')
source('scripts/gscan/popgenome_fun.R')

## Input files:
infile_lookup <- 'radseq/metadata/lookup_IDshort.txt'
infile_scafs <- 'other/seqdata_misc/reference/mmur/scaffolds_withLength2.txt'
if(scaffold != 'all') {
  infile_vcf <- paste0('hybridzone/seqdata/vcf/', fileID_in, '.vcf.gz') # CHANGE TO VCF-DIR
  infile_vcf.tabix <- paste0(infile_vcf, '.tbi')
}

## Output:
outdir <- 'hybridzone/analyses/gscan/popgenome/output/'


################################################################################
##### SET-UP - PROCESS & REPORT #####
################################################################################
## Output files:
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
outfile_base1 <- paste0(outdir, '/', fileID_in, '_win', winsize, '_step', stepsize)
outfile_base2 <- paste0(outdir, '/', fileID_in)

outfile_between_con <- paste0(outfile_base2, '_meanBtwPopStats.txt')
outfile_within_con <- paste0(outfile_base2, '_meanWthPopStats.txt')
outfile_sitepat_smr <- paste0(outfile_base2, '_sitepatSmr.txt')
outfile_sitestats <- paste0(outfile_base2, '_sitestats.txt')
outfile_scafstats <- paste0(outfile_base2, '_scafstats.txt')
outfile_winstats <- paste0(outfile_base1, '_winstats.txt')

outfile_sfs <- paste0(outfile_base2, '_sfs.txt')
outfile_introgstats <- paste0(outfile_base2, '_introgstats.txt')

## Individuals:
lookup <- read.delim(infile_lookup, as.is = TRUE)

## Process pops:
pops <- c(pop1, pop2, pop3, pop4, pop5) %>% .[!is.na(.)]
npop <- length(pops)

## Scaffolds:
scaf_all <- read.table(infile_scafs, header = TRUE,
                       colClasses = c('character', 'integer', 'integer'))

## Report:
cat('#### File ID:', fileID_in, '\n')
cat('#### Scaffold:', scaffold, '\n')
cat('#### Populations:', pops, '\n')
cat('#### Outgroup:', outgroup, '\n')
cat('#### Nr of populations:', npop, '\n')
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
  scaf_len <- scaf_all$length[match(scaffold, scaf_all$scaffold)]

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
## Individuals in VFC:
inds_g <- unique(unlist(g_scaf@region.data@populations2))
inds_g <- inds_g[-grep('\\.2', inds_g)]

## Summarize by region and report:
g_summary <- get.sum.data(g_scaf) %>%
  as.data.frame() %>%
  rename(nsite = n.sites,
         nvar = n.biallelic.sites,
         TiTv = trans.transv.ratio) %>%
  rownames_to_column('scaffold') %>%
  mutate(scaffold = gsub('.vcf.gz', '', scaffold)) %>%
  select(scaffold, nsite, nvar, TiTv)

cat('\n#### Summary of loaded sequence data:\n')
head(g_summary, n = 3)
cat('\n#### Number of scaffolds:', nrow(g_summary), '\n')
cat('#### Total nr of biallelic sites:', sum(g_summary$nvar), '\n')
cat('#### Nr of individuals:', length(inds_g), '\n')

## Vector of scaffolds present in data:
scaffolds_g <- g_summary$scaffold


################################################################################
##### SET POPULATIONS #####
################################################################################
cat("\n#### Setting populations...\n")

## Get populations:
lookup_sel <- filter(lookup, ID.short %in% inds_g) %>%
  select(ID.short, popIDcolumn) %>%
  rename(ID = ID.short, pop = popIDcolumn) %>%
  mutate(pop = gsub('-', '', pop))
#head(lookup_sel)

## All pops -- poplist_main:
make_poplist <- function(pops) {
  if(class(pops) == 'data.frame') pops <- as.character(unlist(pops))
  poplist <- list()
  for(i in 1:length(pops)) {
    fpop <- pops[i]
    fpopname <- paste0('inds_pop', i)
    fpopinds <- lookup_sel$ID[lookup_sel$pop %in% fpop]
    assign(fpopname, fpopinds)
    cat('pop:', fpop, '\ninds:', fpopinds, '\n\n')
    poplist[[i]] <- fpopinds
  }
  names(poplist) <- pops
  return(poplist)
}
poplist_main <- make_poplist(pops)

## Pop triplets:
if(do_introgstats == TRUE) {
  cat('#### Assigning pop triplets...\n')
  #for(triplet in triplets_list) poplist <- make_poplist(triplet)
  triplets_poplists <- lapply(triplets_list, make_poplist)
}

## Pop pairs:
if(npop > 2 & do_pop_pairs == TRUE) {
  cat('#### Assigning pop pairs...\n')
  pop_pairs_poplists <- lapply(pop_pairs_list, make_poplist)
}

## Set populations:
g_scaf <- set.populations(g_scaf, poplist_main, diploid = TRUE)

## Report:
cat("\n#### Populations and nr of inds:\n")
inds_per_pop <- unlist(lapply(g_scaf@populations, length)) / 2
print(inds_per_pop)

## Set outgroup:
if(outgroup != 'none') {
  cat("\n#### Setting outgroup...\n")
  inds_outgroup <- lookup_sel$ID[lookup_sel$pop %in% outgroup]
  g_scaf <- set.outgroup(g_scaf, inds_outgroup, diploid = TRUE)
  cat("\n#### Outgroup:\n")
  print(g_scaf@outgroup)
}


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
if(do_betweenpop == TRUE & do_dxy == TRUE) {
  cat('\n\n#### Calculating Dxy...\n')

  ## Concat:
  g_con <- diversity.stats.between(g_con)
  dxy_con <- g_con@nuc.diversity.between

  dxy_con <- t(g_con@nuc.diversity.between) %>%
    as.data.frame() %>%
    rownames_to_column('pop') %>%
    rename(value = V1) %>%
    mutate(stat = 'dxy',
           pop = namepops_vec(pop)) %>%
    mutate(pop = gsub('/', '.', pop))
  cat("\n#### Mean pairwise Dxy:\n")
  print(dxy_con)

  ## By window:
  g_win <- diversity.stats.between(g_win)
  dxy_win <- g_win@nuc.diversity.between %>%
    namepops(., stat = 'dxy')
  head(dxy_win, n = 2)

  ## By scaffold:
  #g_scaf <- diversity.stats.between(g_scaf)
}


################################################################################
##### BETWEEN-POP: FST (OVERALL AND PER-WINDOW) #####
################################################################################
if(do_betweenpop == TRUE & do_fst_win == TRUE) {
  cat('\n\n#### Calculating FST...\n')
  # get.F_ST(g_con) # doesnt work

  ## Concatenated:
  g_con <- F_ST.stats(g_con, mode = 'nucleotide')

  fst_con_pair <- g_con@nuc.F_ST.pairwise %>%
    as.data.frame() %>%
    rownames_to_column('pop') %>%
    rename(value = V1) %>%
    mutate(stat = 'fst',
           pop = namepops_vec(pop)) %>%
    mutate(pop = gsub('/', '.', pop))
  cat("#### Mean pairwise FST:\n")
  print(fst_con_pair)

  fst_con_single <- g_con@nuc.F_ST.vs.all %>%
    namepops2(., stat = 'fst')
  cat("#### Mean single-pop FST:\n")
  print(fst_con_single)

  ## By window:
  g_win <- F_ST.stats(g_win, mode = 'nucleotide')

  fst_win_pair <- t(g_win@nuc.F_ST.pairwise) %>%
    namepops(., stat = 'fst')
  head(fst_win_pair, n = 2)

  #fst_win_single <- g_win@nuc.F_ST.vs.all %>%
  #  namepops(., stat = 'fstSingle')
  #head(fst_win_single, n = 2)

  ## By scaffold:
  g_scaf <- F_ST.stats(g_scaf, mode = 'nucleotide')
  fst_scaf <- t(g_scaf@nuc.F_ST.pairwise) %>%
    namepops(., stat = 'fst') %>%
    data.frame(scaffold = gsub('.vcf.gz', '', g_scaf@region.names), .)
  head(fst_scaf, n = 4)
}


################################################################################
##### BETWEEN-POP: PER-SITE FST #####
################################################################################

## Doesn't seem to work for multi-pop - only gives single FST value per site
# cat('\n\n#### Calculating per-site FST...\n')
# g_scaf <- detail.stats(g_scaf, site.FST = TRUE)
# fst_site_list <- g_scaf@region.stats@site.FST
#
# fst_site <- as.data.frame(fst_site_list[[1]])
# colnames(fst_site) <- 'fst'
# fst_site$scaffold <- scaffolds_g[1]
#
# for(i in 2:length(fst_site_list)) {
#   if(!is.null(fst_site_list[[i]])) {
#     newrows <- as.data.frame(fst_site_list[[i]])
#     colnames(newrows) <- 'fst'
#     newrows$scaffold <- scaffolds_g[i]
#     fst_site <- rbind(fst_site, newrows)
#   }
# }
#
# fst_site <- fst_site %>%
#   mutate(site = as.integer(names(unlist(fst_site_list)))) %>%
#   select(scaffold, site, contains('fst')) %>%
#   mutate(fst = round(fst, 4)) %>%
#   mutate(fst = ifelse(fst < 0, 0, fst)) # If FST < 0 --> make 0


if(do_betweenpop == TRUE & do_fst_site == TRUE) {
  cat('\n\n#### Calculating per-site FST...\n')

  sitefst_list <- list()
  for(i in 1:length(pop_pairs_poplists)) {
    pop_pair_list <- pop_pairs_poplists[[i]]
    pop_pair_name <- paste0(names(pop_pair_list), collapse = '.')
    cat("\n#### Pop pair:", pop_pair_name, '...\n')

    g_scaf_pair <- detail.stats(g_scaf, site.FST = TRUE,
                                new.populations = pop_pair_list)

    fst_site_list <- g_scaf_pair@region.stats@site.FST
    fst_site_pair <- as.data.frame(fst_site_list[[1]])
    colnames(fst_site_pair) <- 'fst'
    fst_site_pair$scaffold <- scaffolds_g[1]

    for(j in 2:length(fst_site_list)) {
      if(!is.null(fst_site_list[[j]])) {
        newrows <- as.data.frame(fst_site_list[[j]])
        colnames(newrows) <- 'fst'
        newrows$scaffold <- scaffolds_g[j]
        fst_site_pair <- rbind(fst_site_pair, newrows)
      }
    }

    fst_site_pair <- fst_site_pair %>%
      mutate(site = as.integer(names(unlist(fst_site_list)))) %>%
      select(scaffold, site, contains('fst')) %>%
      mutate(fst = round(fst, 4)) %>%
      mutate(fst = ifelse(fst < 0, 0, fst)) # If FST < 0 --> make 0

    colnames(fst_site_pair)[3] <- paste0('fst_', pop_pair_name)
    #head(fst_site_pair)

    sitefst_list[[i]] <- fst_site_pair
  }
  fst_site <- do.call(cbind, sitefst_list)
  fst_site <- fst_site[c(1:2, grep('fst', colnames(fst_site)))]
  head(fst_site)
}


################################################################################
##### BETWEEN-POP: SITE PATTERNS (FIXED, SHARED, ETC) #####
################################################################################
if(do_betweenpop == TRUE & do_sitepat == TRUE) {
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
  rownames(sitepat) <- NULL
  sitepat <- select(sitepat, contains('sipat'))
  head(sitepat, n = 3)

  ## Site pattern summary:
  sitepat_smr <- sitepat %>%
    select(contains('sipat')) %>%
    gather("species", "sipat", 1:npop) %>%
    group_by(species, sipat) %>%
    summarize(count = n()) %>%
    spread(key = species, value = count)
  cat('\n#### Site pattern summary:\n')
  print(sitepat_smr)

  if(npop == 2) { ## TO DO: DO FOR EACH POPULATION PAIR
    stat <- paste0('sipat', sitepat_smr$sipat)
    value <- sitepat_smr[, 2]
    pair <- c('pop1/pop2')
    sitepat_smr2 <- cbind(pair, value, stat)
    colnames(sitepat_smr2)[2] <- 'value'

    cat('\n#### Total nr of sites:', nrow(sitepat), '\n')
    cat('#### Nr of fixed differences:', sum(sitepat[, 1] == 4), '\n')
    cat('#### Nr of shared polymorphisms:', sum(sitepat[, 1] == 0), '\n')
    cat('#### Nr of polymorphisms unique to pop1:', sum(sitepat[, 1] == 1), '\n')
    cat('#### Nr of polymorphisms unique to pop2:', sum(sitepat[, 2] == 1), '\n')
  }
}


################################################################################
#### WITHIN-POP: DIVERSITY ####
################################################################################
if(do_withinpop == TRUE & do_diversity == TRUE) {
  cat('\n\n#### Calculating diversity...\n')

  ## Concatenated:
  g_con <- diversity.stats(g_con, pi = TRUE)
  div_con <- get.diversity(g_con)

  div_smr <- do.call(rbind, lapply(1:npop, function(x) div_con[[x]]))
  div_smr <- as.data.frame(t(data.frame(div_smr, row.names = NULL)))
  colnames(div_smr) <- pops
  div_smr <- as.data.frame(apply(div_smr, 2, round, 5))
  div_smr$stat <- c('nucdiv', 'hapdiv', 'pi', 'hapFst_all', 'fst_all')
  rownames(div_smr) <- NULL

  cat('\n\n#### Summary of diversity stats (concatenated):\n')
  print(div_smr)

  ## By window:
  g_win <- diversity.stats(g_win, pi = TRUE)
  div_win <- get.diversity(g_win)

  pi_win <- get_stat_df(div_win, 'Pi', 'pi')
  nucdiv_win <- get_stat_df(div_win, 'nuc.diversity.within', 'nucdiv')
  pi_win <- cbind(pi_win, nucdiv_win)
  head(pi_win, n = 2)

  ## By scaffold:
  #g_scaf <- diversity.stats(g_scaf, pi = TRUE)
}


################################################################################
##### WITHIN-POP: NEUTRALITY STATS #####
################################################################################
if(do_withinpop == TRUE & do_neutral == TRUE) {
  cat('\n\n#### Calculating neutrality stats...\n')

  ## Concatenated:
  g_con <- neutrality.stats(g_con)
  neutral_con <- get.neutrality(g_con)

  neutral_smr <- do.call(rbind, lapply(1:npop, function(x) neutral_con[[x]]))
  colnames(neutral_smr) # Names of statistics
  neutral_smr <- as.data.frame(t(data.frame(neutral_smr, row.names = NULL)))
  colnames(neutral_smr) <- pops
  neutral_smr <- as.data.frame(apply(neutral_smr, 2, round, 5))
  neutral_smr$stat <- c('tajD', 'nseg', 'rozasR2', 'fuLiF', 'fuLiD',
                        'fuFs', 'fayWuH', 'zengE', 'strobeckS')
  rownames(neutral_smr) <- NULL

  cat('\n\n#### Summary of neutrality stats (concatenated):\n')
  print(neutral_smr)

  ## By window:
  g_win <- neutrality.stats(g_win)
  neutral_win <- get.neutrality(g_win)

  taj_win <- get_stat_df(neutral_win, 'Tajima.D', 'tajD')
  fuLiF_win <- get_stat_df(neutral_win, 'Fu.Li.F', 'fuLiF')
  nseg_win <- get_stat_df(neutral_win, 'n.segregating.sites', 'nseg')
  neutral_win_sel <- cbind(taj_win, fuLiF_win, nseg_win)
  head(neutral_win_sel, n = 2)

  ## Per scaffold:
  #g_scaf <- neutrality.stats(g_scaf)
}

################################################################################
#### WITHIN-POP: MAFs and SFS ####
################################################################################
if(do_withinpop == TRUE & do_sfs == TRUE) {
  cat('\n\n#### Calculating MAFs and SFS...\n')

  ## Overall MAF:
  cat("\n#### Setting overal MAFs (with no-pops list)...\n")
  nopops <- list(inds_g)
  g_no_pops <- detail.stats(g_con, site.spectrum = TRUE,
                            new.populations = nopops)

  MAF_all <- g_no_pops@region.stats@minor.allele.freqs[[1]]
  MAF_all <- as.data.frame(t(MAF_all))
  colnames(MAF_all) <- 'MAF_all'
  MAF_all$MAF_all <- round(MAF_all$MAF_all, 5)
  #head(MAF_all)
  #dim(MAF_all)

  ## Concatenated:
  g_con <- detail.stats(g_con, site.spectrum = TRUE)
  MAF_con <- g_con@region.stats@minor.allele.freqs[[1]]

  MAF <- as.data.frame(t(MAF_con)) %>%
    namepops(., stat = 'MAF')
  MAF <- as.data.frame(apply(MAF, 2, round, 5))
  #head(MAF, n = 2)

  sfs_all <- list()
  for(i in 1:npop) {
    popname <- pops[i]
    popname_in_MAF <- paste('pop', i)

    sfs_pop <- data.frame(table(MAF_con[popname_in_MAF, ])) %>%
      transmute(MAF = as.numeric(as.character(Var1)),
                freq = as.integer(as.character(Freq)))
    colnames(sfs_pop)[2] <- paste0('freq.', popname)

    sfs_all[[i]] <- sfs_pop
  }
  sfs <- Reduce(function(x, y) merge(x, y, all = TRUE),
                sfs_all, accumulate = FALSE)
  sfs[is.na(sfs)] <- 0
  cat('\n\n')
  head(sfs, n = 2)
}


################################################################################
#### INTROGRESSION STATS: D, df, BDF ####
################################################################################
if(do_introgstats == TRUE) {
  cat('\n\n#### Calculating introgression stats...\n')

  ## Df, Fd for triplets:
  triplet_statlist <- list()
  for(i in 1:length(triplets_poplists)) {
    triplet_list <- triplets_poplists[[i]]
    triplet_name <- paste0(names(triplet_list), collapse = '.')
    cat("#### Triplet:", triplet_name, '\n')

    ## Concatenated:
    # g_con_triplet <- set.populations(g_con, new.populations = triplet_list)
    # g_con_triplet <- introgression.stats(g_con_triplet, do.D = TRUE, do.df = TRUE)
    # print(g_con@D) # g_con@D.z, g_con@D.pval
    # print(g_con@f)
    # print(g_con@df) # g_con@df.z, g_con@df.pval
    # ## TO DO: MAKE OVERVIEW DF WITH STATS AND PVALS

    ## By window:
    g_win_triplet <- set.populations(g_win, new.populations = triplet_list)

    g_win_triplet <- introgression.stats(g_win_triplet, do.D = TRUE, do.df = TRUE)
    introg_win <- data.frame(cbind(g_win_triplet@D,
                                   g_win_triplet@f,
                                   g_win_triplet@df))
    colnames(introg_win) <- c('D', 'fd', 'df')
    introg_win <- introg_win %>%
      mutate(fd = ifelse(fd < -1, -1, ifelse(fd > 1, 1, fd)))
    colnames(introg_win) <- paste0(colnames(introg_win), '_', triplet_name)

    triplet_statlist[[i]] <- introg_win
  }
  triplet_stats <- do.call(cbind, triplet_statlist)

  ## RND for pop-pairs:
  pop_pair_statlist <- list()
  if(do_pop_pairs == TRUE) {
    for(i in 1:length(pop_pairs_poplists)) {
      pop_pair_list <- pop_pairs_poplists[[i]]
      pop_pair_name <- paste0(names(pop_pair_list), collapse = '.')
      cat("#### Pop pair:", pop_pair_name, "i:", i, '\n')

      # g_con_pair <- set.populations(g_con, new.populations = pop_pair_list)
      # g_con_pair <- introgression.stats(g_con_pair, do.RNDmin = TRUE,
      #                                   do.D = FALSE, do.df = FALSE)
      # print(g_con_pair@RNDmin)

      g_win_rnd <- set.populations(g_win, new.populations = pop_pair_list)
      g_win_rnd <- introgression.stats(g_win_rnd, do.RNDmin = TRUE,
                                       do.D = FALSE, do.df = FALSE)

      rnd_win_pair <- data.frame(g_win_rnd@RNDmin)
      colnames(rnd_win_pair) <- paste0('RNDmin_', pop_pair_name)

      pop_pair_statlist[[i]] <- rnd_win_pair
    }
    RND_win <- do.call(cbind, pop_pair_statlist)
    #head(RND_win)
  }
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
if(do_withinpop == TRUE) {
  within_con <- rbind(div_smr, fst_con_single, neutral_smr)
  print(within_con)
}

## Overall between-pop stats (averaged across all data):
if(do_betweenpop == TRUE & npop == 2) {
  between_con <- rbind(fst_con_pair, dxy_con, sitepat_smr2)
}
if(do_betweenpop == TRUE & npop > 2) {
  between_con <- rbind(fst_con_pair, dxy_con)
}

between_con <- between_con %>%
  select(stat, pop, value) %>%
  mutate(value = round(value, 4))
print(between_con)

## Per-site stats:
if(do_sfs == TRUE & do_sitepat == TRUE) {
  #sitepat_and_MAF <- cbind(sitepat, MAF, MAF_all)
  sitestats <- cbind(fst_site, sitepat) %>%
    rename(start = site) %>%
    mutate(start = as.integer(start)) %>%
    mutate(end = start)
  head(sitestats)

  MAF_merged <- cbind(MAF, MAF_all) ## DIFFERENT NR OF ROWS THAN OTHER SITESTATS
}

## Per-scaffold stats:
scafstats <- g_summary ## add stats?
head(scafstats)

## Per-window stats:
scaffold_vec <- fst_site$scaffold[win_start]
scaf_win_end <- fst_site$scaffold[win_end]
win_remove <- which(scaffold_vec != scaf_win_end) # Windows that cross scaffold boundaries will be removed
cat('### Number of windows to be removed:', length(win_remove), '\n')
start_vec <- fst_site$site[win_start]
end_vec <- fst_site$site[win_end]
winstats <- cbind(scaffold_vec, start_vec, end_vec,
                  pi_win, neutral_win_sel,
                  fst_win_pair, dxy_win) #fst_win_single
if(do_introgstats == TRUE) winstats <- cbind(winstats, triplet_stats, RND_win)
winstats <- winstats %>%
  slice(-win_remove) %>%
  rename(scaffold = scaffold_vec,
         start = start_vec,
         end = end_vec)
head(winstats)
dim(winstats)


################################################################################
#### SANITY-CHECK RESULTS #####
################################################################################
if(do_summaries == TRUE) {
  sitesmr <- sitestats %>%
    group_by(scaffold) %>%
    summarize(max_site = max(start),
              fst_site = round(mean(fst), 3)) %>%
    merge(., g_summary, by = 'scaffold') %>%
    mutate(gap = nsite - max_site) %>%
    select(scaffold, nsite, nvar, fst_site, gap)
  sitesmr

  all(sitesmr$gap == 0)
  if(all(sitesmr$gap == 0)) sitesmr <- select(sitesmr, -gap)

  winsmr <- winstats %>%
    group_by(scaffold) %>%
    summarize(nwin = n(),
              fst_win = round(mean(fst_gri.mur), 3))

  smr <- merge(sitesmr, winsmr, by = 'scaffold') %>%
    merge(., fst_scaf, by = 'scaffold') %>%
    merge(., scaf_all[, c('scaffold', 'length')], by = 'scaffold') %>%
    mutate(gap2 = length - nsite)

  head(smr)
  plot(smr$fst_site, smr$fst_win)
  plot(smr$fst_scaf, smr$fst_win)
  plot(smr$nsite, smr$length)
  plot(smr$length, smr$nvar)
}


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
