# SET-UP -----------------------------------------------------------------------
source("scripts/treemix/treemix_plotting_funcs.R") # Script that comes with Treemix
source("scripts/treemix/treemix_plotFun.R") # My script
library(tidyverse)
library(ggpubr)
library(cowplot)
library(png)
library(grid)

## Metadata:
infile_lookup <- "metadata/lookup_IDshort.txt"
lookup <- read.delim(infile_lookup, as.is = TRUE)

## Settings
file.ID <- "hzproj1.mac1.FS6"
root <-"None"
dir.output <- "output/bySupersite2/"
poporder <- unique(lookup$supersite2)
dir.fig <- "figures/supersite2/"
nmig.vector <- 0:8
file.open <- TRUE

## Dirs
if(!dir.exists(dir.fig)) dir.create(dir.fig)
if(!dir.exists("popfiles")) dir.create("popfiles")
if(!dir.exists("LRT")) dir.create("LRT")

## Final output figure
figfile <- paste0("figures/", file.ID, "combined.png")

## Process
poporder.file <- paste0("popfiles/", file.ID, "_poporder.txt")
writeLines(poporder, poporder.file)


# PROPORTION OF VARIANCE EXPLAINED ---------------------------------------------
cat("##### Df with proportion of variance explained: \n")
(propVar.df <- get.propVar.df(file.ID, dir.output, root, nmig.max = 8))
propVar.plot <- plot.propVar(propVar.df, file.ID, dir.fig, file.open = TRUE)


# LIKELIHOOD PLOT AND LRT ------------------------------------------------------
cat("##### Df with likelihoods: \n")
llh.df <- get.llh.df(file.ID, root, nmig.vector, dir.output)
llh.plot <- plot.llh(llh.df, file.ID, dir.fig, file.open = TRUE)


# TREE PLOT FOR EACH VALUE OF M ------------------------------------------------
cat("Plotting trees: \n")
sapply(nmig.vector, plot.tree,
       file.ID, root, dir.output, dir.fig,
       png.background = "white", file.open = file.open, filetype = "png")


# PLOT RESIDUALS ---------------------------------------------------------------
## Positive residuals: candidates for admixture
sapply(nmig.vector, plot.residuals,
       file.ID, poporder.file, dir.output, dir.figs, file.open = file.open)


# COMBINE PLOTS ----------------------------------------------------------------
treemix.plot <- rasterGrob(readPNG(treemix.plot.file), interpolate = TRUE)
treemix.plot <- ggplot(data.frame()) +
  geom_blank() +
  annotation_custom(treemix.plot, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

plots1 <- ggarrange(propVar.plot, llh.plot, ncol = 1, nrow = 2, heights = c(1, 1.2))
plots <- ggarrange(plots1, treemix.plot, ncol = 2, nrow = 1, widths = c(1, 1.5))
plots <- plots + draw_plot_label(label = c("A", "B", "C"), size = 24,
                                 x = c(0, 0, 0.4), y = c(1, 0.55, 0.9))

ggexport(plots, filename = figfile, width = 900, height = 650)
