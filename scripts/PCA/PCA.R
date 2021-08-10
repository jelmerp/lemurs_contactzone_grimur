# SET-UP -----------------------------------------------------------------------

## Script with plotting functions
source("scripts/PCA/PCA_R_fun.R")

## Packages
library(ggrepel)

## Input files and dirs
infile_inds <- here("metadata/hzlookup_bysample.txt")
infile_qc <- here('results/qc/qc_byind_r03.txt')
indir_vcf <- here("results/geno/vcf/gatk")

# Output files
outdir_figs <- here('results/figs_ms/')
outdir_snpgds <- here('results/PCA/gdsFiles/')

## Read metadata
qc <- read.delim(infile_qc, as.is = TRUE) %>%
  filter(repd == FALSE) %>%
  mutate(filter_pass = factor(filter_pass, levels = c('FS6', 'FS7', 'rescue', 'fail'))) %>%
  select(ID, filter_pass, bam_dp, vcf_dp, vcf_miss_pct, vcf_nSNP)
inds_df <- read.delim(infile_inds, as.is = TRUE) %>%
  select(ID, sp, site, pop, poptype, put_hyb)
  left_join(., qc, by = "ID")

## Colors
col_mur <- '#FF00FF'
col_gri <- '#D69B12'
col_hyb <- 'black'

## Create output dirs if necessary
if(!dir.exists(outdir_snpgds)) dir.create(outdir_snpgds, recursive = TRUE)
if(!dir.exists(outdir_figs)) dir.create(outdir_figs, recursive = TRUE)

## Plot prep
species_order <- c('mgri', 'mmur', 'mhyb')
col_by_labs <- c(expression(italic("griseorufus")),
                 expression(italic("murinus")), 'hybrid')
shape_by_labs <- c('parapatric', 'sympatric')


# FS7 --------------------------------------------------------------------------

fileID <- "r03.all.mac1.FS7"
infile_vcf <- here(indir_vcf, paste0(fileID, '.vcf.gz'))
outfile_snpgds <- here(outdir_snpgds, paste0(fileID, '.gds'))

snps <- snps_get(infile_vcf, outfile_snpgds)
inds_vcf <- read.gdsn(index.gdsn(snps, "sample.id")) # Get a list of individuals in the VCF

## Run PCA and process results
pca_raw <- snpgdsPCA(snps, autosome.only = FALSE, num.thread = 1)
pca <- process_pca(pca_raw, my_lookup = inds_df) %>%
  mutate(sp_msat = ifelse(put_hyb == FALSE, sp, "mhyb"),
         sp_msat = factor(sp_msat, levels = species_order)) %>%
  arrange(sp_msat)
eig_all <- get_eig(pca_raw) # Get eigenvalues

## Plot
p <- pcplot(pca,
             eigenvals = eig_all$eig,
             col_by = 'sp_msat',
             col_by_name = 'microsats:',
             col_by_labs = col_by_labs,
             cols = c(col_gri, col_mur, col_hyb),
             shape_by = 'poptype',
             shape_by_name = 'site type:',
             shape_by_labs = shape_by_labs,
             shapes = c(1, 0))

p + geom_label_repel(aes(label = filter_pass))


# KEEPHYBS ---------------------------------------------------------------------

fileID <- "r03.keepHybs.mac1.FS6"
infile_vcf <- here(indir_vcf, paste0(fileID, '.vcf.gz'))
outfile_snpgds <- here(outdir_snpgds, paste0(fileID, '.gds'))

snps <- snps_get(infile_vcf, outfile_snpgds)
inds_vcf <- read.gdsn(index.gdsn(snps, "sample.id")) # Get a list of individuals in the VCF

## Run PCA and process results
pca_raw <- snpgdsPCA(snps, autosome.only = FALSE, num.thread = 1)
pca <- process_pca(pca_raw, my_lookup = inds_df) %>%
  mutate(sp_msat = ifelse(put_hyb == FALSE, sp, "mhyb"),
         sp_msat = factor(sp_msat, levels = species_order)) %>%
  arrange(sp_msat)
eig_all <- get_eig(pca_raw) # Get eigenvalues

## Plot
p <- pcplot(pca,
            eigenvals = eig_all$eig,
            col_by = 'sp_msat',
            col_by_name = 'microsats:',
            col_by_labs = col_by_labs,
            cols = c(col_gri, col_mur, col_hyb),
            shape_by = 'poptype',
            shape_by_name = 'site type:',
            shape_by_labs = shape_by_labs,
            shapes = c(1, 0))

p + geom_label_repel(aes(label = ID))


# COMBINE PLOTS AND SAVE FIGURE ------------------------------------------------

## Arrange  panels
p <- ggarrange(p_mgri, p_mmur,
               ncol = 1, nrow = 2, heights = c(1, 1)) +
  draw_plot_label(label = c('A', 'B'), size = 25,
                  x = c(0, 0), y = c(1, 0.5))

## Save figure
figfile_eps <- here(outdir_figs, paste0(fileID, '_singlespecies.eps'))
ggsave(figfile_eps, p, width = 6, height = 10)