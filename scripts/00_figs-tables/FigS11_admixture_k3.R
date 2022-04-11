# SET-UP -----------------------------------------------------------------------

## Packages
library(here)

## Source script with functions
source(here("scripts/admixture/admixture_plot_fun.R"))

## Define input dir and files
infile_inds <- here("metadata/hzlookup_bysample.txt")
infile_colors <- here("metadata/colors.txt")

admixture_outdir <- here("results/admixture/output")

## Define output files
figdir <- here("results/figs_ms")

## Read metadata
color_df <- read.delim(infile_colors, as.is = TRUE) %>%
  select(-species)
inds_df <- read.delim(infile_inds, as.is = TRUE) %>%
  merge(., color_df, by.x = "sp", by.y = "pop") %>%
  rename(labcol = col) %>%
  mutate(sp_msat = ifelse(put_hyb == FALSE, sp, "hybrid?"),
         sp_msat = gsub("^m", "", sp_msat),
         poptype = ifelse(poptype == "parapatric" & sp == "mmur", "zparapatric", poptype),
         poptype = ifelse(!grepl("patric", poptype), NA, poptype))

## Colors
col_mur <- "#FF00FF"
col_gri <- "#D69B12"

## SetID
setID <- "r03.all.mac3.FS6"


# PLOT FOR K=3 -----------------------------------------------------------------

k2 <- Qdf(setID, inds_df, K = 2, sort_by = "sp") %>%
  ggax_v(.,
         group_column = c("sp_msat", "poptype"),
         barcols = c(col_gri, col_mur),
         grouplab_bgcol = my_grouplab_bgcol,
         grouplab_labeller = grouplab_labeller,
         ylab = "RADseq: K=2")

k3 <- Qdf(setID, inds_df, K = 3, sort_by = "sp") %>%
  ggax_v(.,
         group_column = c("sp_msat", "poptype"),
         barcols = c(col_mur, col_gri, "grey40"),
         strip_show = FALSE,
         ylab = "RADseq: K=3")

mt <- Qdf(setID, inds_df, K = 2, sort_by = "sp") %>%
  ggax_v(.,
         group_column = c("sp_msat", "poptype"),
         barcols = c(col_gri, col_mur),
         strip_show = FALSE,
         ylab = "mtDNA")


# COMBINE PLOTS AND SAVE FIG ---------------------------------------------------

p <- ggarrange(k2, k3, mt,
               ncol = 1, nrow = 3, heights = c(1, 0.8, 0.3))

figfile <- paste0(figdir, "/", setID, "_K3_combined.eps")
ggsave(p, filename = figfile, width = 8, height = 6)
