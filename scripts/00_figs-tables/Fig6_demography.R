# SET-UP --------------------------------------------------------------------

## Packages
library(here)

## Other scripts
source(here('scripts/gphocs-bpp/gphocs_plot_fun.R'))    # Plotting functions for GPHOCS results
source(here('scripts/gphocs-bpp/bpp_plot_fun.R'))       # Plotting functions for BPP results

## Input files
infile_logs <- here('results/gphocs/output/gphocs-bpp_mergedlogs.txt')
infile_colors <- here('metadata/demography_colors.txt')

## Output files
fig_prefix <- here('results/figs_ms/', 'Fig5_demography')

## Read metadata - population colors
pop_colors <- read.delim(infile_colors, as.is = TRUE)

## Define populations
poplist <- list(murW = 'A', B = 'A',
                murC = 'B', murE = 'B',
                griC = 'C', griW = 'C',
                A = 'root', C = 'root')
childpops <- names(poplist)
parentpops <- unique(as.character(poplist))
allpops <- unique(c(childpops, parentpops))
currentpops <- setdiff(allpops, parentpops)
poporder <- c('murW', 'murC', 'murE', 'B', 'A', 'griC', 'griW', 'C', 'root')
migorder <- c('griC_2_murC', 'murC_2_griC', 'C_2_B', 'B_2_C', 'C_2_A', 'A_2_C')

## Define order of run IDs
runIDs <- c('noMig', 'g2m2anc', 'bpp_6mig')

## Load and prep GPHOCS and BPP output files
logs <- as.data.frame(fread(infile_logs, stringsAsFactors = FALSE)) %>%
  mutate(pop = factor(pop, levels = allpops),
         runID = factor(runID, levels = runIDs),
         migpattern = factor(migpattern, levels = migorder))

## Plotting labels for populations
poplabs_m <- c('gri-C > mur-C', 'mur-C > gri-C',
               'C > B', 'B > C', 'C > A', 'A > C')
poplabs_th <- c('mur-W', 'mur-C', 'mur-E', 'A', 'B', 'gri-C', 'gri-W', 'C', 'root')


# DEMOGRAPHY OVERVIEW PLOTS ----------------------------------------------------

d1 <- dplotwrap(logs, runID_f = 'noMig',
                poplist = poplist, popcols = pop_colors, poporder = poporder,
                y.max = 1000, extra.time.root = 1000, popnames.size = 5,
                xticks.by = 50, pop.spacing = 150,
                rm.y.ann = FALSE, xlab = '', plot.title = 'GPhoCS: isolation')

d2 <- dplotwrap(logs, runID_f = 'g2m2anc',
                poplist = poplist, popcols = pop_colors, poporder = poporder,
                y.max = 1000, extra.time.root = 1000, popnames.size = 5,
                xticks.by = 50, pop.spacing = 150,
                rm.y.ann = TRUE, xlab = NULL, plot.title = 'GPhoCS: migration') +
  theme(axis.title.x = element_text(margin = margin(1, 0, 0, 0, 'cm')))

d3 <- dplotwrap(logs, runID_f = 'bpp_6mig',
                poplist = poplist, popcols = pop_colors, poporder = poporder,
                y.max = 1000, extra.time.root = 1000, popnames.size = 5,
                xticks.by = 50, pop.spacing = 150,
                rm.y.ann = TRUE, xlab = '', plot.title = 'BPP (migration)')


# MIGRATION PLOTS --------------------------------------------------------------

m_bpp <- vplot(
  data = filter(logs, var == 'phi', !is.na(migpattern)),
  xvar = 'migpattern', fillvar = 'cn', colvar = 'cn', yvar = 'val',
  xlab = "", ylab = 'phi (BPP)', legpos = 'top', shade_by = 2,
  yticks.by = 0.02, rotate.x.ann = TRUE, linecols = 'black'
) +
  scale_x_discrete(labels = poplabs_m) +
  geom_vline(xintercept = 2.5, color = "grey30", size = 1) +
  geom_vline(xintercept = 4.5, color = "grey30", size = 1) +
  theme(plot.margin  = margin(0.2, 0.2, 0, 0.5, 'cm'))

m_gphocs <- vplot(
  data = filter(logs, var == '2Nm'),
  xvar = 'migpattern', fillvar = 'cn', colvar = 'cn', yvar = 'val',
  xlab = "", ylab = '2Nm (G-PhoCS)', legpos = 'top', shade_by = 2,
  yticks.by = 0.02, rotate.x.ann = TRUE, linecols = 'black'
) +
  scale_x_discrete(labels = poplabs_m) +
  geom_vline(xintercept = 2.5, color = "grey30", size = 1) +
  geom_vline(xintercept = 4.5, color = "grey30", size = 1) +
  theme(plot.margin  = margin(0.2, 0.2, 0, 0.5, 'cm'))


# COMBINE, FINALIZE AND SAVE ---------------------------------------------------

## Set colors
popname_cols <- c(rep('#FF00FF', 3), rep('#D69B12', 2))

## Combine plots
p <- (d1 + d2 + d3) / (m_bpp + m_gphocs) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 24, face = 'bold'))

## Add labels
p <- ggarrange(p) +
  draw_plot_label(label = c('mur-W', 'mur-E', 'mur-C', 'gri-C', 'gri-W'),
                  colour = popname_cols, size = 16, fontface = 'plain',
                  y = c(rep(0.575, 5)),
                  x = c(0.09, 0.155, 0.21, 0.27, 0.32)) +
  draw_plot_label(label = c('mur-W', 'mur-E', 'mur-C', 'gri-C', 'gri-W'),
                  colour = popname_cols, size = 16, fontface = 'plain',
                  y = c(rep(0.575, 5)),
                  x = c(0.09, 0.155, 0.21, 0.27, 0.32) + 0.302) +
  draw_plot_label(label = c('mur-W', 'mur-E', 'mur-C', 'gri-C', 'gri-W'),
                  colour = popname_cols, size = 16, fontface = 'plain',
                  y = c(rep(0.575, 5)),
                  x = c(0.09, 0.155, 0.21, 0.27, 0.32) + 0.606)

## Save final plot
ggsave(paste0(fig_prefix, '.png'), p, width = 12, height = 14)
ggsave(paste0(fig_prefix, '.eps'), p, width = 12, height = 14,
       device = cairo_ps, fallback_resolution = 150)
