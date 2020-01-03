## TO DO: INCLUDE PERCENT DUPLICATES + INCLUDE TRIMSTATS
## aa_dedupStats.sum.r99_201807.txt
## aa_trimStats.sum.r99_201807.txt

################################################################################
##### SET-UP  #####
################################################################################
library(tidyverse)
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')
source('/home/jelmer/Dropbox/sc_lemurs/radseq/scripts/qc/qc_fun.R')

## IDs:
fileID_finalVcf <- 'r03.all.mac1.FS6'
fileID_rescueVcf <- 'r03.keepHybs.mac1.FS6'
cat("##### Final VCF ID:", fileID_finalVcf, '\n')

## Input files:
infile_lookup <- '/home/jelmer/Dropbox/sc_lemurs/radseq/metadata/lookup_IDshort.txt'
IDs.long <- readLines('metadata/r03/sampleIDs_r03.txt')
IDs.short <- readLines('metadata/r03/sampleIDsShort_r03.txt')

dir_bamstats <- 'analyses/qc/bam/'
dir_fastqstats <- 'analyses/qc/fastq/byInd/'
dir_vcf.joint <- 'analyses/qc/vcf/map2mmur.gatk4.paired.joint/vcftools/'

infile_prop.pairs <- 'analyses/qc/bam/aa_r03_propPairedReads.txt'
infile_bam.dp <- 'analyses/qc/bam/aa_r03_meanDepth.txt'
infile_idepth.ind <- 'analyses/qc/vcf/map2mmur.gatk4.paired.ind/vcftools/aa_r03_idepth.txt'

infile_idepth.joint <- paste0(dir_vcf.joint, fileID_finalVcf, '.idepth')
infile_imiss.joint <- paste0(dir_vcf.joint, fileID_finalVcf, '.imiss')
infile_idepth.rescue <- paste0(dir_vcf.joint, fileID_rescueVcf, '.idepth')
infile_imiss.rescue <- paste0(dir_vcf.joint, fileID_rescueVcf, '.imiss')

## Output files:
outfile_stats.df <- 'analyses/qc/r03_qc_summary.txt'
outfile_stats.rescue <- 'analyses/qc/r03_rescue_qc_summary.txt'

## Read metadata:
lookup <- read.delim(infile_lookup, header = TRUE, as.is = TRUE)

## Bamstats column names:
bamstats.colnames.initial <- readLines('analyses/qc/bam/eautils_colnames.txt')
bamstats.colnames.final <- gsub(' ', '.', bamstats.colnames.initial)
bamstats.colnames.final <- gsub('%', 'pct.', bamstats.colnames.final)


################################################################################
##### FASTQ STATS #####
################################################################################
fastqstats.R1 <- read.fastqstats.all(IDs.long, read = 'R1')
fastqstats.R2 <- read.fastqstats.all(IDs.long, read = 'R2')
fastqstats <- rbind(fastqstats.R1, fastqstats.R2) %>%
  group_by(ID) %>%
  summarise(fastq.reads = sum(fastq.reads),
            fastq.qual = round(mean(fastq.qual), 4),
            fastq.pct.dup = round(mean(fastq.pct.dup), 2),
            fastq.meanlen = round(mean(fastq.meanlen), 2))
head(fastqstats)


################################################################################
##### BAM: NR OF PROP PAIRED READS, DEPTH, MQ  #####
################################################################################
## Samtools flagstats:
flagstats <- read.flagstats.all(IDs.long)
head(flagstats)

## Mean depth:
bam.dp <- read.delim(infile_bam.dp, header = FALSE, as.is = TRUE, sep = " ",
                      col.names = c('ID', 'bam.dp')) %>%
  filter(grepl('_', ID)) %>%
  mutate(bam.dp = round(bam.dp, 2))
head(bam.dp)

## Mapping quality, etc:
bamstats.b <- read.bamstats.all(IDs.long) %>%
  select(ID, reads, pct.ambiguous, mapq.mean) %>%
  rename(bam.eau.reads = reads, bam.eau.pct.ambig = pct.ambiguous, bam.eau.mapq = mapq.mean)
head(bamstats.b)

## Merge:
bamstats.a <- merge(flagstats, bam.dp, by = 'ID', all.x = TRUE)
head(bamstats.a)
bamstats <- merge(bamstats.a, bamstats.b, by = 'ID', all.x = TRUE)
head(bamstats)


################################################################################
##### VCF  #####
################################################################################
## idepth for by-ind-VCFs:
idepth.ind <- read.delim(infile_idepth.ind, header = FALSE, as.is = TRUE,
                         colClasses = c('character', 'integer', 'numeric'),
                         col.names = c('ID', 'ivcf.nsite', 'ivcf.dp')) %>%
  mutate(ivcf.dp = round(ivcf.dp, 2))

## idepth for joint-VCFs:
idepth.joint <- read.delim(infile_idepth.joint, header = TRUE, as.is = TRUE,
                           colClasses = c('character', 'integer', 'numeric'),
                           col.names = c('ID', 'jvcf.nsite', 'jvcf.dp')) %>%
  mutate(jvcf.dp = round(jvcf.dp, 2))

vcfstats <- merge(idepth.ind, idepth.joint, by = 'ID', all.x = TRUE)
inds.pass <- vcfstats$ID[!is.na(vcfstats$jvcf.dp)]
cat("##### Number of individuals in final joint VCF:", length(inds.pass), '\n')

## imiss for joint-VCFs:
imiss.joint <- read.delim(infile_imiss.joint, header = TRUE, as.is = TRUE,
                          colClasses = c('character', 'integer', 'integer', 'integer', 'numeric'),
                          col.names = c('ID', 'jvcf.nsite.all', 'jvcf.nfiltered',
                                        'jvcf.nmiss', 'jvcf.fmiss')) %>%
  select(ID, jvcf.nsite.all, jvcf.fmiss) %>%
  mutate(jvcf.fmiss = round(jvcf.fmiss, 3))

nsite.vcfJoint.all <- imiss.joint$jvcf.nsite.all[1]
imiss.joint$jvcf.nsite.all <- NULL
cat("##### Number of SNPs in final joint VCF:", nsite.vcfJoint.all, '\n')

vcfstats <- merge(vcfstats, imiss.joint, by = 'ID', all.x = TRUE)
head(vcfstats)


################################################################################
##### VCF STATS FOR "RESCUED" INDIVS #####
################################################################################
idepth.rescue <- read.delim(infile_idepth.rescue, header = TRUE, as.is = TRUE,
                            colClasses = c('character', 'integer', 'numeric'),
                            col.names = c('ID', 'jvcf.nsite', 'jvcf.dp'))

inds.rescue <- idepth.rescue$ID[! idepth.rescue$ID %in% inds.pass]
cat("##### Number of rescued individuals in rescue-VCF:", length(inds.rescue), '\n')

imiss.rescue <- read.delim(infile_imiss.rescue, header = TRUE, as.is = TRUE,
                           colClasses = c('character', 'integer', 'integer', 'integer', 'numeric'),
                           col.names = c('ID', 'jvcf.nsite.all', 'jvcf.nfiltered',
                                         'jvcf.nmiss', 'jvcf.fmiss'))

stats.rescue <- merge(idepth.rescue, imiss.rescue) %>%
  mutate(jvcf.fmiss = round(jvcf.fmiss, 3),
         Sample_ID = lookup$Sample_ID[lookup$ID.short %in% idepth.rescue$ID]) %>%
  select(ID, Sample_ID, jvcf.nsite, jvcf.dp, jvcf.nsite.all, jvcf.fmiss) %>%
  filter(ID %in% inds.rescue)


################################################################################
##### MERGE DFs #####
################################################################################
## Merge fastq and bam stats:
stats.df.tmp <- merge(fastqstats, bamstats, by = 'ID')

## GET RID OF BOTH LIBRARIES FOR TECH-DUP INDS, FOR STATS:
stats.df.tmp$ID <- substr(stats.df.tmp$ID, 1, 7)
to_exclude <- stats.df.tmp$ID[which(duplicated(stats.df.tmp$ID))]
stats.df.tmp <- stats.df.tmp %>% filter(! ID %in% to_exclude)

## Merge with VCF stats:
stats.df <- merge(stats.df.tmp, vcfstats, by = 'ID', all.x = TRUE) %>%
  mutate(
    Sample_ID = lookup$Sample_ID[lookup$ID.short %in% stats.df.tmp$ID],
    pass = ifelse(is.na(jvcf.nsite), 0, 1),
    rescue = ifelse(ID %in% inds.rescue, 1, 0),
    species.short = substr(ID, 1, 4),
    pct.mapped = round(bam.flag.mapped / bam.flag.total, 4),
    pct.proppaired = round(bam.flag.proppaired / bam.flag.paired, 4),
    fastq.score = (fastq.reads * fastq.qual) / 1000000,
    ivcf.score = round((ivcf.nsite * ivcf.dp) / 1000, 1),
    jvcf.score = jvcf.dp * (1 - jvcf.fmiss),
    ivcf.score.by.fastq.score = ivcf.score / fastq.score,
    jvcf.score.by.fastq.score = jvcf.score / fastq.score
    )
head(stats.df)
#sum(stats.df$pass)

fstats <- stats.df %>% filter(pass == 1)


################################################################################
##### PLOT #####
################################################################################
xlabs <- c(expression(italic("griseorufus")), 'hybrid', expression(italic("murinus")))
xtitle <- 'Previous msat assignment'

ggbox('fastq.reads', xlabs = xlabs, xtitle = xtitle,
      ytitle = 'Number of reads in fastq', ptitle = 'FASTQ: number of reads')
ggbox('pct.mapped', xlabs = xlabs, xtitle = xtitle,
      ytitle = 'Percentage of successfully mapped reads', ptitle = 'BAM: mapped reads')
ggbox('pct.proppaired', xlabs = xlabs, xtitle = xtitle,
      ytitle = 'Percentage of properly paired reads', ptitle = 'BAM: properly paired reads')
ggbox('bam.eau.mapq', xlabs = xlabs, xtitle = xtitle,
      ytitle = 'Mean mapping quality', ptitle = 'BAM: mapping quality')
ggbox('bam.dp', xlabs = xlabs, xtitle = xtitle,
      ytitle = 'Mean depth in BAM', ptitle = 'BAM: depth')
ggbox('jvcf.dp', xlabs = xlabs, xtitle = xtitle,
      ytitle = 'Mean depth in joint-VCF', ptitle = 'Joint-VCF: depth')
ggbox('jvcf.fmiss', xlabs = xlabs, xtitle = xtitle,
      ytitle = 'Percentage missing data in joint-VCF', ptitle = 'Joint-VCF: missing data')
ggbox('ivcf.nsite', xlabs = xlabs, xtitle = xtitle,
      ytitle = 'Number of sites in individual-VCF', ptitle = 'Individual-VCF: number of sites')
ggbox('ivcf.dp', xlabs = xlabs, xtitle = xtitle,
      ytitle = 'Mean depth in individual-VCF', ptitle = 'Individual-VCF: depth')
ggbox('ivcf.score', xlabs = xlabs, xtitle = xtitle,
      ytitle = 'Score (depth + nr sites) for individual-VCF', ptitle = 'Individual-VCF: score')
ggbox('jvcf.score', xlabs = xlabs, xtitle = xtitle,
      ytitle = 'Score (depth + % missing) for joint-VCF', ptitle = 'Joint-VCF: score')
ggbox('fastq.score', xlabs = xlabs, xtitle = xtitle,
      ytitle = 'Score (qual + nr sites) for FASTQ', ptitle = 'FASTQ: score')
ggbox('ivcf.score.by.fastq.score', xlabs = xlabs, xtitle = xtitle,
      ytitle = 'Individual-VCF score over FASTQ score\n(~ yield per equivalent read)',
      ptitle = 'Indidividual-VCF: "yield"')
ggbox('jvcf.score.by.fastq.score', xlabs = xlabs, xtitle = xtitle,
      ytitle = 'Joint-VCF score over FASTQ score\n(~ yield per equivalent read)',
      ptitle = 'Joint-VCF: "yield"')


################################################################################
##### SUMMARIZE #####
################################################################################
## Somehow, bam.flag.paired == fastq.reads

fstats %>%
  group_by(species.short) %>%
  summarize_all(funs(mean)) %>%
  t()

## Linear model?
#readnrs.lm <- lm(reads.passed ~ genus + library + barcode, data = inds.df)
#(anova.res <- anova(readnrs.lm))

################################################################################
##### WRITE FILES #####
################################################################################
# write.table(stats.df, outfile_stats.df, sep = '\t', quote = FALSE, row.names = FALSE)
# write.table(stats.rescue, outfile_stats.rescue, sep = '\t', quote = FALSE, row.names = FALSE)
