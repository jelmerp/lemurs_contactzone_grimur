# SET-UP - GENERAL -------------------------------------------------------------

## Load packages
library(tidyverse)
library(ggpubr)
library(cowplot)
library(here)

## Define metadata input files
infile_inds <- here("metadata/hzlookup_bysample.txt")

## Define output files
figdir <- here("results/figs_ms")
figfile_eps <- here(figdir, "Fig2_structure.eps")

## Read metadata
inds_df <- read.delim(infile_inds, header = TRUE, as.is = TRUE) %>%
  mutate(sp_msat = gsub("hybrid", "hybrid?", sp_msat),
         sp_msat = gsub("^m", "", sp_msat),
         poptype = ifelse(poptype == "parapatric" & sp == "mmur", "zparapatric", poptype),
         poptype = ifelse(!grepl("patric", poptype), NA, poptype))

## Set colors
col_mur <- "#FF00FF"
col_gri <- "#D69B12"
col.hyb <- "grey50"


# SET-UP - ADMIXTURE -----------------------------------------------------------

## Plotting functions
source("scripts/admixture/admixture_plot_fun.R")

## Input dir
admixture_res_dir <- "results/admixture/output/"

## Dataset ID
setID <- "r03.all.mac3.FS6"

## Plot prep
grouplab_bgcol <- c(col_gri, col_gri, col_hyb, col_mur, col_mur)
poptype_labs <- c(parapatric = "para", sympatric = "sym", zparapatric = "para")
grouplab_labeller <- labeller(poptype = poptype_labs)


# ADMIXTURE PLOTS --------------------------------------------------------------

## CV-plot
p_cv <- CVplot(setID) +
  labs(y = "CV error") +
  theme(axis.title.x = element_text(margin = unit(c(0.1, 0, 0, 0), "cm")),
        axis.title.y = element_text(margin = unit(c(0, 0.1, 0, 0), "cm")),
        plot.title = element_blank(),
        plot.margin = margin(1, 0.2, 0.5, 0.2, "cm"))

## Q-plots
p_k2 <- Qdf(setID, inds_df, K = 2, sort_by = "sp") %>%
  ggax_v(.,
         group_column = c("sp_msat", "poptype"),
         barcols = c(col_gri, col_mur),
         grouplab_bgcol = grouplab_bgcol,
         grouplab_labeller = grouplab_labeller,
         ylab = "RADseq")

p_mt <- Qdf(setID, inds_df, K = 2, sort_by = "sp") %>%
    ggax_v(.,
           group_column = c("sp_msat", "poptype"),
           barcols = c(col_gri, col_mur),
           strip_show = FALSE,
           ylab = "mtDNA")

p_k <- ggarrange(p_k2, p_mt, ncol = 1, nrow = 2, heights = c(1, 0.25))
p_admix <- ggarrange(p_cv, p_k, nrow = 2, heights = c(0.4, 0.6))


# PCA - SET-UP -----------------------------------------------------------------

## Plotting functions
source("scripts/PCA/PCA_R_fun.R")

## Settings
fileID <- "r03.all.mac1.FS6"

## Input VCF
infile_vcf <- here(indir_vcf, paste0(fileID, ".vcf.gz"))

## Run PCA
outfile_snpgds <- here(outdir_snpgds, paste0(fileID, ".gds"))
snps <- snps_get(infile_vcf, outfile_snpgds)
pca_raw <- snpgdsPCA(snps, autosome.only = FALSE, num.thread = 1)
pca <- process_pca(pca_raw, my_lookup = inds_df) %>%
  arrange(sp_msat)
eig_df <- get_eig(pca_raw) # Get eigenvalues


# PCA - PLOT -------------------------------------------------------------------

col_by_labs <- c(expression(italic("griseorufus")),
                 expression(italic("murinus")), "hybrid?")
shape_by_labs <- c("parapatric", "sympatric")

p_pca <- pcplot(pca.df,
                eigenvals = eig_df$eig,
                col_by = "sp_msat",
                col_by_name = "microsats:",
                col_by_labs = col_by_labs,
                cols = c(col_gri, col_mur, col_hyb),
                shape_by = "poptype",
                shape_by_name = "site type:",
                shape_by_labs = shape_by_labs,
                shapes = c(1, 0),
                legpos = "top") +
  theme(legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = margin(0.6, 0.2, 0, 1, "cm"))


# MAPS - SET-UP ---------------------------------------------------------------

## Packages
library(ggmap)
library(ggsn) #devtools::install_github("oswaldosantos/ggsn")

## Input files
infile_coords_Mtk <- here("metadata/coordinates/Mangatsiaka_coords_dec.txt")
infile_coords_Tsm <- here("metadata/coordinates/Tsimelahy_coords_dec.txt")

## Prep sampling locations
coords_Mtk <- read.delim(infile_coords_Mtk, header = TRUE, as.is = TRUE)
coords_Tsm <- read.delim(infile_coords_Tsm, header = TRUE, as.is = TRUE)

## Set plotting labels
sp.labs <- c(expression(italic("griseorufus")), expression(italic("   murinus")))
hyb.labs <- c("hybrid", "non-hybrid")
sp.cols <- c("#D69B12", "#FF00FF")
msat.name <- "\nmicrosat\nassignment:"
radseq.name <- "\nRADseq\nassignment:"

## Mangatsiaka lat and lon:
lat.center.Mtk <- -24.96275
lon.center.Mtk <- 46.55692
lon.min.Mtk <- 46.55267
lon.max.Mtk <- 46.56099
lat.min.Mtk <- -24.96914
lat.max.Mtk <- -24.96000
lon.diff.Mtk <- 46.56099 - 46.55267
lon.min.Mtk.scalebar <- lon.min.Mtk + (0.05 * (lon.max.Mtk - lon.min.Mtk))
lat.max.Mtk.scalebar <- lat.max.Mtk - (0.04 * (lat.max.Mtk - lat.min.Mtk))
lon.max.Mtk.scalebar <- lon.max.Mtk
lat.min.Mtk.scalebar <- lat.min.Mtk

## Tsimelahy lat and lon:
lat.center.Tsm <- -24.95425
lon.center.Tsm <- 46.61117
lon.min.Tsm <- 46.61004
lon.max.Tsm <- 46.62222
lat.min.Tsm <- -24.95884
lat.max.Tsm <- -24.94711
lon.min.Tsm.scalebar <- lon.min.Tsm + (0.05 * (lon.max.Tsm - lon.min.Tsm))
lat.max.Tsm.scalebar <- lat.max.Tsm - (0.04 * (lat.max.Tsm - lat.min.Tsm))
lon.max.Tsm.scalebar <- lon.max.Tsm
lat.min.Tsm.scalebar <- lat.min.Tsm

ggmap::register_google(key = "AIzaSyB0hxxyLujH5apnsRW21X5MlqO8_LnbjGc")


# PLOT MAPS ----------------------------------------------------------------

## Mangatsiaka
map_Mtk <- get_googlemap(center = c(lon = lon.center.Mtk, lat = lat.center.Mtk),
                         zoom = 15, # 3=continent, 10=city, 21=building
                         scale = 2,
                         maptype = "satellite", # roadmap, terrain, hybrid, satellite
                         color = "color")

p_Mtk <- ggmap(map_Mtk, extent = "device") +
  geom_point(data = coords_Mtk,
             aes(x = lon, y = lat, shape = msat.hybrid, fill = species.cor),
             colour = "black",
             size = 4,
             stroke = 1) +
  scale_x_continuous(limits = c(lon.min.Mtk, lon.max.Mtk), expand = c(0, 0)) +
  scale_y_continuous(limits = c(lat.min.Mtk, lat.max.Mtk), expand = c(0, 0)) +
  scale_fill_manual(values = sp.cols) +
  scale_shape_manual(values = c(25, 21)) +
  labs(title = "sympatric site 1: Mangatsiaka") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "plain"),
        legend.position = "none",
        panel.border = element_rect(colour = "grey20", fill = NA, size = 2),
        plot.margin = unit(c(1, 0.2, 0.2, 0.2), "cm")) +
  ggsn::scalebar(x.min = lon.min.Mtk.scalebar,
                 x.max = lon.max.Mtk.scalebar,
                 y.min = lat.min.Mtk.scalebar,
                 y.max = lat.max.Mtk.scalebar,
                 dist = 150,
                 dist_unit = "m",
                 location = "topleft",
                 transform = TRUE,
                 model = "WGS84",
                 st.size = 5,
                 border.size = 0.8,
                 box.fill = c("gray60", "white"),
                 st.color = "white")

## Tsimelahy
map_Tsm <- get_googlemap(center = c(lon = lon.center.Tsm, lat = lat.center.Tsm),
                         zoom = 15,
                         scale = 2,
                         maptype = "satellite",
                         color = "color")

p_Tsm <- ggmap(map_Tsm, extent = "device") +
  geom_point(data = coords_Tsm,
             aes(x = lon, y = lat, shape = msat.hybrid, fill = species.cor),
             colour = "black",
             size = 4,
             stroke = 1) +
  scale_x_continuous(limits = c(lon.min.Tsm, lon.max.Tsm),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(lat.min.Tsm, lat.max.Tsm),
                     expand = c(0, 0)) +
  scale_shape_manual(values = c(25, 21),
                     labels = hyb.labs,
                     name = msat.name) +
  scale_fill_manual(values = sp.cols,
                    labels = sp.labs,
                    name = radseq.name) +
  guides(fill = guide_legend(override.aes = list(shape = 21),
                             label.hjust = -1)) +
  labs(title = "sympatric site 2: Tsimelahy") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "plain"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold"),
        legend.key.size = unit(0.6, "cm"),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-5, 0, -5, -5),
        panel.border = element_rect(colour = "grey20", fill = NA, size = 2),
        plot.margin = unit(c(1, 0, 0.2, 0.2), "cm")) +
  scalebar(x.min = lon.min.Tsm.scalebar,
           x.max = lon.max.Tsm.scalebar,
           y.min = lat.min.Tsm.scalebar,
           y.max = lat.max.Tsm.scalebar,
           dist = 150,
           dist_unit = "m",
           location = "topleft",
           transform = TRUE,
           model = "WGS84",
           st.size = 5,
           border.size = 0.8,
           box.fill = c("gray60", "white"),
           st.color = "white")

p_map <- ggarrange(p_Mtk, p_Tsm,
                   ncol = 2, nrow = 1, widths = c(1, 1.5), heights = c(1, 1))


# COMBINE PLOTS ----------------------------------------------------------------

p_top <- ggarrange(p_admix, p_pca, ncol = 2, widths = c(0.5, 0.5))
p <- ggarrange(p_top, p_map, ncol = 1, nrow = 2, heights = c(0.5, 0.5))
p <- p +
  draw_plot_label(label = c("A", "B", "C"),
                  size = 24,
                  x = c(0, 0.52, 0), y = c(1, 1, 0.48))

## Save plot
ggsave(p, filename = figfile_eps, width = 11, height = 11)
