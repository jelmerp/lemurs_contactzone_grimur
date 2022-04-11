# SET-UP -----------------------------------------------------------------------

## Packages
library(here)

## Source script with functions
source(here("scripts/admixture/admixture_plot_fun.R"))

## Define input dir and files
infile_inds <- here("metadata/hzlookup_bysample.txt")
infile_qc <- here('results/qc/qc_byind_r03.txt')
admixture_outdir <- here("results/admixture/output")

## Define output files
figdir <- here("results/manuscript")

## Read metadata
qc <- read.delim(infile_qc, as.is = TRUE) %>%
  mutate(filter_pass = factor(filter_pass, levels = c('FS6', 'FS7', 'rescue', 'fail'))) %>%
  select(ID, filter_pass)
inds_df <- read.delim(infile_inds, as.is = TRUE) %>%
  select(ID, sp, site, pop, poptype, sp_mtDNA, sp_msat = sp_msat_overall) %>%
  mutate(sp_msat = gsub("hybrid", "hybrid?", sp_msat),
         sp_msat = gsub("^m", "", sp_msat),
         poptype = ifelse(poptype == "parapatric" & sp == "mmur", "zparapatric", poptype),
         poptype = ifelse(!grepl("patric", poptype), NA, poptype)) %>%
  left_join(., qc, by = "ID") %>%
  filter(filter_pass %in% c('FS6', 'FS7', 'rescue'))

# FS6_QC <- filter(inds_df, filter_pass == "FS6") %>% pull(ID)
# FS6_admix <- unique(Qdf("r03.all.mac1.FS6", inds_df, K = 2)$ID)
#
# FS7_QC <- c(FS6_QC, filter(inds_df, filter_pass == "FS7") %>% pull(ID))
# FS7_admix <- unique(Qdf("r03.all.mac1.FS7", inds_df, K = 2)$ID)
#
# rescue_QC <- c(FS6_QC, filter(inds_df, filter_pass == "rescue") %>% pull(ID))
# rescue_admix <- unique(Qdf("r03.keepHybs.mac1.FS6", inds_df, K = 2)$ID)

## Plot settings
col_mur <- "#FF00FF"
col_gri <- "#D69B12"

my_grouplab_bgcol <- c(col_gri, col_gri, "grey50", col_mur, col_mur)
poptype_labs <- c(parapatric = "parapatric", sympatric = "sympatric",
                  zparapatric = "parapatric")
grouplab_labeller <- labeller(poptype = poptype_labs)

setID_FS6 <- "r03.all.mac1.FS6"
setID_FS7 <- "r03.all.mac1.FS7"
setID_rescue <- "r03.keepHybs.mac1.FS6"


# CREATE PLOTS -----------------------------------------------------------------

mt_df <- rbind(
  mutate(inds_df,
         cluster = "V1",
         proportion = ifelse(sp_mtDNA == "mgri", 1, 0)),
  mutate(inds_df,
         cluster = "V2",
         proportion = ifelse(sp_mtDNA == "mgri", 0, 1))
) %>%
  arrange(ID)


## mtDNA
mt <- ggax_v(mt_df,
             group_column = c("sp_msat", "poptype"),
             barcols = c(col_gri, col_mur),
             grouplab_bgcol = my_grouplab_bgcol, grouplab_size = 10,
             grouplab_labeller = grouplab_labeller,
             ylab = "mt")

## FS6
fs6 <- Qdf(setID_FS6, inds_df, K = 2,
           sort_by = "sp", keep_all_in_lookup = TRUE) %>%
  mutate(ID = as.character(ID)) %>%
  arrange(ID) %>%
  ggax_v(.,
         group_column = c("sp_msat", "poptype"),
         barcols = c(col_mur, col_gri),
         strip_show = FALSE,
         ylab = "FS6")

## FS7
fs7 <- Qdf(setID_FS7, inds_df, K = 2,
           sort_by = "sp", keep_all_in_lookup = TRUE) %>%
  mutate(ID = as.character(ID)) %>%
  arrange(ID) %>%
  ggax_v(.,
         group_column = c("sp_msat", "poptype"),
         barcols = c(col_gri, col_mur),
         strip_show = FALSE,
         ylab = "FS7")

## KEEP-HYBRIDS
rescue <- Qdf(setID_rescue, inds_df, K = 2,
              sort_by = "sp", keep_all_in_lookup = TRUE) %>%
  mutate(ID = as.character(ID)) %>%
  arrange(ID) %>%
  ggax_v(.,
         group_column = c("sp_msat", "poptype"),
         barcols = c(col_gri, col_mur),
         strip_show = FALSE,
         ylab = "rescue")

# Qdf(setID, inds_df, K = 2, sort_by = "sp") %>% filter(proportion < 0.99, proportion > 0.9)
#> mhyb003 mgri104


# COMBINE PLOTS ----------------------------------------------------------------

p <- ggarrange(mt, fs6, fs7, rescue,
               ncol = 1, nrow = 4, heights = c(1, 1, 1, 1))

figfile_eps <- here(figdir, "FigS13_admixture_allfilters.eps")
ggsave(p, filename = figfile_eps, width = 9, height = 5)

figfile_png <- here(figdir, "FigS13_admixture_allfilters.png")
ggsave(p, filename = figfile_png, width = 10, height = 5.6)
