# SET-UP -----------------------------------------------------------------------

## Script with plotting functions
source("scripts/PCA/PCA_R_fun.R")

## Variables
fileID <- "r03.all.mac3.FS6"
ID.type <- "ID.short"

## Input files and dirs
infile_inds <- here("metadata/hzlookup_bysample.txt")

indir_vcf <- here("data/vcf/gatk")
infile_vcf <- here(indir_vcf, paste0(fileID, '.vcf.gz'))

# Output files
outdir_figs <- here('results/figs_ms/')
outdir_snpgds <- here('results/PCA/gdsFiles/')
outfile_snpgds <- here(outdir_snpgds, paste0(fileID, '.gds'))

## Read metadata
inds_df <- read.delim(infile_inds, as.is = TRUE)

## Colors
col_mur <- '#FF00FF'
col_gri <- '#D69B12'
col_hyb <- 'black'

## Create output dirs if necessary
if(!dir.exists(outdir_snpgds)) dir.create(outdir_snpgds, recursive = TRUE)
if(!dir.exists(outdir_figs)) dir.create(outdir_figs, recursive = TRUE)

## Read VCF
snps <- snps_get(infile_vcf, outfile_snpgds)

## Get a list of individuals in the VCF
inds_vcf <- read.gdsn(index.gdsn(snps, "sample.id"))


# PCA FOR MGRI ONLY ------------------------------------------------------------

subset_id <- 'mgri_cz'

## Select samples
hybs_mgri <- c('mhyb001', 'mhyb003', 'mhyb004', 'mhyb005', 'mhyb006',
               'mhyb011', 'mhyb015', 'mhyb016')
inds_mgri <- c(inds_vcf[grepl('mgri', inds_vcf)],
               inds_vcf[inds_vcf %in% hybs_mgri])

## Run PCA
pca_raw <- snpgdsPCA(snps, autosome.only = FALSE, num.thread = 1,
                     sample.id = inds_mgri)
pca <- process_pca(pca_raw, my_lookup = inds_df)
eig <- get_eig(pca_raw)

## Plot labels
labs_species <- c(expression(italic("griseorufus")), 'hybrid')
title_gri <- expression(paste(italic('griseorufus'), ' only'))

## Plot:
p_mgri <- pcplot(pca,
                 eigenvals = eig$eig,
                 col_by = 'site',
                 col_by_name = 'site:',
                 shape_by = 'species.short',
                 shape_by_name = 'microsats:',
                 shape_by_labs = labs_species,
                 shapes = c(1, 0),
                 plot_title = title_gri) +
  theme(plot.margin = margin(0.2, 0.2, 0.6, 0.2, "cm"))


# PCA FOR MURINUS ONLY ---------------------------------------------------------

subset_id <- 'mmur_cz'

## Select samples
hybs_mmur <- c('mhyb002', 'mhyb007', 'mhyb008', 'mhyb009',
               'mhyb010', 'mhyb012', 'mhyb013', 'mhyb014')
inds_mmur <- c(inds_vcf[grepl('mmur', inds_vcf)],
               inds_vcf[inds_vcf %in% hybs_mmur])

## Run PCA
pca_raw <- snpgdsPCA(snps, autosome.only = FALSE, num.thread = 1,
                     sample.id = inds_mmur)
pca <- process_pca(pca_raw, my_lookup = inds_df) %>%
  mutate(sp = factor(sp, levels = c('mmur', 'mhyb')))
eig <- get_eig(pca_raw)

## Plot labels
labs_species <- c(expression('hybrid', italic("murinus")))
title_mur <- expression(paste(italic('murinus'), ' only'))

## Plot:
p_mmur <- pcplot(pca,
                 eigenvals = eig$eig,
                 col_by = 'site',
                 col_by_name = 'site:',
                 shape_by = 'species.short',
                 shape_by_name = 'microsats:',
                 shape_by_labs = labs_species,
                 shapes = c(1, 0),
                 plot.title = title_mur) +
  theme(plot.margin = margin(0.6, 0.2, 0.2, 0.2, "cm"))


# COMBINE PLOTS AND SAVE FIGURE ------------------------------------------------

## Arrange  panels
p <- ggarrange(p_mgri, p_mmur,
               ncol = 1, nrow = 2, heights = c(1, 1)) +
  draw_plot_label(label = c('A', 'B'), size = 25,
                  x = c(0, 0), y = c(1, 0.5))

## Save figure
figfile_eps <- here(outdir_figs, paste0(fileID, '_singlespecies.eps'))
ggsave(figfile_eps, p, width = 6, height = 10)
