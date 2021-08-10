#### LOAD PACKAGES -------------------------------------------------------------
if (!'pacman' %in% rownames(installed.packages())) install.packages('pacman')
library(pacman)
packages <- c('gridExtra', 'grid', 'RColorBrewer', 'scales',
              'ggpubr', 'cowplot', 'ggforce', 'patchwork',
              'forcats', 'here', 'tidyverse')
p_load(char = packages, install = TRUE)


#### GET K-VALUES FOR WHICH OUTPUT FILES ARE PRESENT ---------------------------
get_Ks <- function(
  setID,
  filedir = here('results/admixture/output')
  ) {

  cat("## Filedir:", filedir, '\n')

  Q_files <- list.files(filedir, pattern = paste0(setID, '.[0-9+].Q'))
  K <- as.integer(gsub(paste0(setID, '.*\\.([0-9]+)\\.Q$'), '\\1', Q_files))

  cat("## K's found:\n")
  print(K)

  return(K)
}


#### K CROSS-VALIDATION PLOT --------------------------------------------------

CVplot <- function(
  setID,
  title = '',
  filedir = here('results/admixture/output/'),
  figsave = FALSE,
  figdir = 'results/admixture/figures/'
  ) {

  K <- get_Ks(setID, filedir = filedir)
  CV_files <- list.files(filedir, pattern = paste0(setID, '.[0-9+].admixtureOutLog.txt'))

  get_CV <- function(K) {
    K_filename <- list.files(filedir, full.names = TRUE,
                             pattern = paste0(setID, '.', K, '\\.admixtureOutLog.txt'))
    K_file <- readLines(K_filename)
    CV <- gsub('.*: (.*)', '\\1', K_file[grep('CV', K_file)])
    return(CV)
  }
  CV <- sapply(sapply(K, get_CV), '[', 1)
  CV_df <- arrange(data.frame(K, CV), K)
  CV_df$CV <- as.numeric(as.character(CV_df$CV))

  p <- ggplot(CV_df, aes(x = K, y = CV, group = 1)) +
    geom_point() +
    geom_line(color = "grey40") +
    scale_x_continuous(breaks = 1:nrow(CV_df)) +
    labs(x = 'K (number of clusters)',
         y = 'Cross-validation error',
         title = title) +
    theme_bw() +
    theme(
      legend.text = element_text(size = 16),
      legend.spacing.x = unit(0, 'cm'),
      legend.spacing.y = unit(0.1, 'cm'),
      legend.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"),
      legend.title = element_text(size = 16, face = 'bold'),
      legend.background = element_rect(fill = "grey90", colour = "grey30"),
      legend.key = element_rect(fill = "grey90"),
      legend.key.size = unit(1, "cm"),
      legend.text.align = 0,
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.ticks.length = unit(.25, "cm"),
      axis.title.x = element_text(size = 16, margin = unit(c(4, 0, 0, 0), 'mm')),
      axis.title.y = element_text(size = 14, margin = unit(c(0, 4, 0, 0), 'mm')),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 20),
      plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm")
      )

  if (figsave == TRUE) {
    if (!dir.exists(figdir)) dir.create(figdir, recursive = TRUE)
    figfile <- paste0(figdir, '/', setID, '.CVplot.png')
    cat('## Saving to file:', figfile, '\n')
    ggsave(figfile, p, width = 5, height = 5)
  }

  return(p)
}


#### GET Q-DATAFRAME -----------------------------------------------------------

Qdf <- function(
  setID, lookup,
  keep_all_in_lookup = FALSE,
  K = 2,
  sort_by = 'ID', ID_column = 'ID', toShortID = FALSE,
  Qdf_dir = 'results/admixture/output/'
  ) {

  #setID <- setID_all; K=3; Qdf_dir = 'results/admixture/output/'

  cat("## K:", K, '\n')

  Qfile <- list.files(Qdf_dir, full.names = TRUE,
                      pattern = paste0(setID, ".", K, ".Q"))
  cat("## Qfile:", Qfile, '\n')

  indlist_file <- list.files(Qdf_dir, full.names = TRUE,
                             pattern = paste0(setID, ".indivs.txt"))
  indlist <- readLines(indlist_file)
  if (toShortID == TRUE) indlist <- substr(indlist, 1, 7)

  if (!all(indlist %in% lookup[, ID_column])) {
    cat("## NOT ALL INDS IN ADMIXTURE OUTPUT ARE FOUND IN LOOKUP - MISSING INDS:\n")
    print(indlist[! indlist %in% lookup[, ID_column]])
  }

  Q <- read.table(Qfile, as.is = TRUE) %>%
    cbind(indlist) %>%
    merge(lookup, ., by.x = ID_column, by.y = 'indlist', all.x = keep_all_in_lookup) %>%
    arrange(!!(rlang::sym(sort_by)))

  cat("## Dimensions of Q df:", dim(Q), '\n')

  Q <- Q[dim(Q)[1]:1,]
  Q[, ID_column] <- factor(Q[, ID_column], levels = Q[, ID_column][1:nrow(Q)])
  Q <- gather(Q, key = 'cluster', value = 'proportion', matches("^V[0-9]+")) %>%
    arrange(!!(rlang::sym(sort_by)))

  return(Q)
}


# VERTICAL BARPLOT -------------------------------------------------------------

ggax_v <- function(
  Q,
  barcols = NULL,
  ID_column = 'ID', prop_column = 'proportion', group_column = 'site',
  indlabs_show = FALSE, indlab_cols = NULL, indlab_column = NULL, indlab_size = 10,
  grouplab_cols = 'grey10', grouplab_size = 14,
  grouplab_labeller = NULL, grouplab_bgcol = 'white', grouplab_angle = 0,
  strip_show = TRUE,
  mar = c(0.1, 0.1, 0.1, 0.1),
  ylab = NULL, title = ''
  ) {

  ## Prep:
  stackedbar_column <- 'cluster'

  ## Plot:
  p <- ggplot(Q,
              aes_(x = as.name(ID_column),
                   y = as.name(prop_column),
                   fill = as.name(stackedbar_column))) +
    geom_bar(position = 'stack', stat = 'identity') +
    guides(fill = "none") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(mar[1], mar[2], mar[3], mar[4], 'cm'))

  ## Indiv sample labels
  if (indlabs_show == FALSE) {
    p <- p + theme(axis.text.x = element_blank())
  } else {
    p <- p + theme(
      axis.text.x = element_text(size = indlab_size, face = 'bold',
                                 angle = 90, hjust = 1)
    )
  }

  ## Colors for bars:
  if (!is.null(barcols)) p <- p + scale_fill_manual(values = barcols)

  ## y-label:
  if (!is.null(ylab)) {
    p <- p + theme(axis.title.y = element_text(angle = 90, size = 15))
    p <- p + labs(y = ylab)
  }

  ## Ind-labels - text and color:
  # if (indlabs_show == TRUE) {
  #   if (!is.null(indlab_column)) {
  #     indlabs <- Q %>% filter(cluster == 'V1') %>% pull(!!(sym(indlab_column)))
  #     cat('## Number of indlabs:', length(indlabs), '\n')
  #
  #     p <- p + scale_x_discrete(expand = c(0, 0), labels = indlabs)
  #   }
  #   if (!is.null(indlab_cols)) {
  #     indlab_cols <- as.character(Q[, indlab_cols])
  #
  #     p <- p + theme(
  #       axis.text.x = element_text(size = indlab_size, face = 'bold',
  #                                  angle = 90, hjust = 1,
  #                                  color = as.character(indlab_cols))
  #       )
  #   } else {
  #     p <- p + theme(
  #       axis.text.x = element_text(size = indlab_size, face = 'bold',
  #                                  angle = 90, hjust = 1)
  #       )
  #   }
  # }

  ## Facet grid:
  if (!is.null(grouplab_labeller)) {
    p <- p +
      facet_grid(as.formula(paste0("~", paste0(group_column, collapse = "+"))),
                 scales = "free_x", space = "free_x",
                 labeller = grouplab_labeller)
  } else {
    p <- p +
      facet_grid(as.formula(paste0("~", paste0(group_column, collapse = "+"))),
                 scales = "free_x", space = "free_x")
  }

  ## Top strips - basic text and colours:
  if (strip_show == TRUE) {
    p <- p + theme(
      strip.text = element_text(size = grouplab_size,
                                colour = grouplab_cols,
                                angle = grouplab_angle,
                                face = 'bold')
    )
  } else {
    p <- p + theme(strip.text = element_blank())
  }

  ## Top strips - background colour:
  if (length(grouplab_bgcol) == 1) {
    if (strip_show == TRUE) {
      p <- p + theme(strip.background = element_rect(fill = grouplab_bgcol))
      }
  } else {
    if (strip_show == TRUE) {
      p <- p + theme(strip.background = element_rect(fill = 'gray90'))
      g <- ggplot_gtable(ggplot_build(p))
      strip_both <- which(grepl('strip-', g$layout$name))
      k <- 1

      for (i in strip_both) {
        j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- grouplab_bgcol[k]
        k <- k + 1
      }

      p <- g
    }
  }
  if (strip_show == FALSE) p <- p + theme(strip.background = element_blank())

  grid.draw(p)
  return(p)
}


# HORIZONTAL BARPLOT -----------------------------------------------------------

ggax_h <- function(
  Q,
  barcols = NULL, labcols = NULL,
  indlab_column = NULL, indlab_size = 20, indlab_firstOnly = TRUE,
  ID_column = 'ID', prop_column = 'proportion',
  title = ''
  ) {

  if (!is.null(labcols)) labcols <- as.character(Q[, labcols])

  p <- ggplot(Q, aes_(x = as.name(ID_column), y = as.name(prop_column))) +
    geom_bar(stat = 'identity', aes(fill = cluster)) +
    coord_flip() +
    guides(fill = "none") +
    labs(title = title) +
    theme(title = element_text(size = 20, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line = element_blank()) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0))

  ## Text and color for (ind-)labels:
  if (!is.null(indlab_column)) {
    my.indlabs <- Q %>% filter(cluster == 'V1') %>% pull(!!(sym(indlab_column)))
    cat('Number of indlabs:', length(my.indlabs), '\n')

    if (indlab_firstOnly == TRUE) {
      indlabs.org <- my.indlabs
      my.indlabs <- rep("", length(indlabs.org))
      last.indices <- length(indlabs.org) - match(unique(indlabs.org), rev(indlabs.org)) + 1
      my.indlabs[last.indices] <- indlabs.org[last.indices]
    }

    cat('Changed indlabs:', my.indlabs, '\n')
    p <- p + scale_x_discrete(expand = c(0, 0), labels = my.indlabs)
  }

  if (is.null(labcols)) {
    p <- p + theme(
      axis.text.y = element_text(size = indlab_size)
      )
  } else {
    p <- p + theme(
      axis.text.y = element_text(size = indlab_size, color = as.character(labcols))
      )
  }

  ## colors for bars:
  if (!is.null(barcols)) p <- p + scale_fill_manual(values = barcols)

  ## Save plot:
  print(p)
}


# GET INTERMEDIATE COLOUR ------------------------------------------------------

get_midcol <- function(col1, col2) {
  col <- rgb(red = (col2rgb(col1)[1] + col2rgb(col2)[1]) / 2,
             green = (col2rgb(col1)[2] + col2rgb(col2)[2]) / 2,
             blue = (col2rgb(col1)[3] + col2rgb(col2)[3]) / 2,
             maxColorValue = 255)
  return(col)
}
