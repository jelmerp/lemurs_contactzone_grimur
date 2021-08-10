#### PACKAGES ------------------------------------------------------------------
library(pacman)
packages <- c('gdsfmt', 'SNPRelate', 'vcfR',
              'ggpubr', 'cowplot', 'ggforce', 'patchwork',
              'here', 'forcats', 'tidyverse')
p_load(char = packages, install = TRUE)


# LOAD VCF ---------------------------------------------------------------------

snps_get <- function(infile_vcf, outfile_snpgds,
                     force_write = FALSE, ld_treshold = 1) {

  outdir_snpgds <- dirname(outfile_snpgds)
  if(!dir.exists(outdir_snpgds)) dir.create(outdir_snpgds, recursive = TRUE)

  if((file.exists(outfile_snpgds) & force_write == TRUE) || !file.exists(outfile_snpgds)) {
    # Convert VCF to snpgds format, which is what SNPrelate uses within R
    snpgdsVCF2GDS(infile_vcf, outfile_snpgds, method = "biallelic.only")
  } else {
    message("## snpgds file exists, not overwriting...")
  }

  snps <- snpgdsOpen(outfile_snpgds, allow.duplicate = TRUE) # Conversion above was to a file, which is now loaded

  cat("## Summary of VCF:\n")
  snpgdsSummary(outfile_snpgds) # Show summary of data in VCF/snpgds

  ## LD pruning:
  if(ld_treshold < 1) {
    snpset <- snpgdsLDpruning(snps, ld.threshold = ld_treshold)
    snpset.id <- unlist(unname(snpset))
    pca <- snpgdsPCA(snpgds.file, snp.id = snpset.id, num.thread = 2)
  }

  return(snps)
}


# SELECT INDIVIDUALS -----------------------------------------------------------

snps_indsel <- function(snps, inds_rm) {

  inds_all_raw <- read.gdsn(index.gdsn(snps, "sample.id"))
  message("## Number of individuals in VCF: ", length(inds_all_raw))
  print(inds_all_raw)

  if(!is.null(inds_rm)) {
    message("## Removing the following individuals:")
    print(inds_rm)

    inds_all <- inds_all_raw[! inds_all_raw %in% inds_rm]

    message("## Number of individuals left: ", length(inds_all))
    stopifnot(length(inds_all) == (length(inds_all_raw) - length(inds_rm)))
  } else {
    inds_all <- inds_all_raw
  }

  message("## List of (remaining) individuals:")
  print(inds_all)
  return(inds_all)
}


# PROCESS PCA RESULTS ----------------------------------------------------------

process_pca <- function(pca_results,
                        my_lookup,
                        pca_type = 'snpgdsPCA') {

  ## Reorganize df:
  if(pca_type == 'glPCA') {
    pca <- data.frame(ID = rownames(pca_results$scores),
                      PC1 = pca_results$scores[, 'PC1'],
                      PC2 = pca_results$scores[, 'PC2'],
                      PC3 = pca_results$scores[, 'PC3'],
                      PC4 = pca_results$scores[, 'PC4'],
                      stringsAsFactors = FALSE, row.names = NULL)
    eigenvals <- pca_results$eig # Eigenfactor df
  }

  if(pca_type == 'snpgdsPCA') {
    pca <- data.frame(ID = pca_results$sample.id,
                      PC1 = pca_results$eigenvect[, 1],
                      PC2 = pca_results$eigenvect[, 2],
                      PC3 = pca_results$eigenvect[, 3],
                      PC4 = pca_results$eigenvect[, 4],
                      stringsAsFactors = FALSE)
    eigenvals <- pca_results$eigenval # Eigenfactor df
  }

  ## Merge with metadata df:
  if(any(! pca$ID %in% my_lookup$ID)) {
    warning("## Not all PCA IDs found in lookup!")
    not_found <- pca$ID[which(! pca$ID %in% my_lookup$ID)]
    message("## The following IDs were not found:")
    print(not_found)
  }
  pca <- merge(pca, my_lookup, by = 'ID')

  message("## Number of rows in pca df: ", nrow(pca))
  return(pca)
}


# GET EIGENVALUES --------------------------------------------------------------

get_eig <- function(pca_results) {
  eig <- data.frame(PC = 1:length(pca_results$eigenval),
                    eig = pca_results$eigenval)
  eig <- eig[complete.cases(eig), ]
}


# PLOT THE PCA -----------------------------------------------------------------

pcplot <- function(pca_df,
                   eigenvals = NULL,
                   pca_id = NULL,
                   pcX = 1,
                   pcY = 2,
                   col_by = NULL,
                   col_by_name = NULL,
                   col_by_labs = NULL,
                   cols = NULL,
                   draw_boxes = FALSE,
                   boxlabsize = 10,
                   shape_by = NULL,
                   shape_by_name = NULL,
                   shape_by_labs = NULL,
                   shapes = NULL,
                   my_shape = 21, # for single shape
                   dotsize = 3,
                   strokesize = 2.5,
                   alpha = 1,
                   legpos = 'right',
                   legstacking = 'vertical',
                   axis_title_size = 16,
                   xmin_buffer = 0.05, xmax_buffer = 0.05,
                   ymin_buffer = 0.05, ymax_buffer = 0.05,
                   plot_title = NULL,
                   plot_save = FALSE,
                   my_outdir_figs = 'tmp') {

  if(!dir.exists(my_outdir_figs)) dir.create(my_outdir_figs, recursive = TRUE)
  if(is.null(pca_id)) pca_id <- 'PCA'

  xcolumn <- paste0('PC', pcX)
  ycolumn <- paste0('PC', pcY)

  if(!is.null(eigenvals)) {
    var.x <- round(eigenvals[pcX] / sum(eigenvals, na.rm = TRUE) * 100, 1)
    var.y <- round(eigenvals[pcY] / sum(eigenvals, na.rm = TRUE) * 100, 1)
    lab.x <- paste0('PC', pcX, ' (', var.x, '%)')
    lab.y <- paste0('PC', pcY, ' (', var.y, '%)')
  } else {
    lab.x <- xcolumn
    lab.y <- ycolumn
  }

  ## Limits:
  xmin_val <- min(pca_df[, xcolumn])
  xmax_val <- max(pca_df[, xcolumn])
  xrange <- xmax_val - xmin_val
  xmin <- xmin_val - (xmin_buffer * xrange)
  xmax <- xmax_val + (xmax_buffer * xrange)

  ymin_val <- min(pca_df[, ycolumn])
  ymax_val <- max(pca_df[, ycolumn])
  yrange <- ymax_val - ymin_val
  ymin <- ymin_val - (ymin_buffer * yrange)
  ymax <- ymax_val + (ymax_buffer * yrange)

  ## By shape and colour:
  if(!is.null(shape_by) & !is.null(col_by)) {
    p <- ggplot(
      data = pca_df,
      aes_string(x = xcolumn, y = ycolumn, color = col_by, shape = shape_by)
    ) +
    geom_point(size = dotsize, stroke = strokesize, alpha = alpha)

    if(!is.null(shapes)) {
      p <- p + scale_shape_manual(
        values = shapes, name = shape_by_name, labels = shape_by_labs
        )
    } else {
      p <- p + scale_shape_discrete(
        name = shape_by_name, labels = shape_by_labs
        )
    }

    if(!is.null(cols)) {
      p <- p + scale_color_manual(
        values = cols, name = col_by_name, labels = col_by_labs
        )
    } else {
      p <- p + scale_color_discrete(
        name = col_by_name
        )
    }
  }

  ## By colour only:
  if(is.null(shape_by) & !is.null(col_by)) {
    cat("## Plotting by colour only\n")
    cat("## Colours:", cols, "\n")

    p <- ggplot(
      data = pca_df,
      aes_string(x = xcolumn, y = ycolumn, color = col_by)
    ) +
      geom_point(size = dotsize, stroke = strokesize, shape = my_shape,
                 alpha = alpha, fill = 'black')

    if(!is.null(col_by_name))
      p <- p + scale_color_manual(
        values = cols, name = col_by_name, labels = col_by_labs
        )

    if(is.null(col_by_name))
      p <- p + scale_color_manual(
        values = cols, labels = col_by_labs
      )
  }

  ## By shape only:
  if(!is.null(shape_by) & is.null(col_by)) {
    p <- ggplot(
      data = pca_df,
      aes_string(x = xcolumn, y = ycolumn, shape = shape_by)
    ) +
      geom_point(size = dotsize, stroke = strokesize, alpha = alpha)

    if(!is.null(shapes)) {
      p <- p + scale_shape_manual(
        values = shapes, name = shape_by_name, labels = shape_by_labs)
    } else {
      p <- p + scale_shape_discrete(
        name = shape_by_name, labels = shape_by_labs)
    }
  }

  if(draw_boxes == TRUE) {
    p <- p + geom_mark_ellipse(
      aes_string(label = col_by),
      show.legend = FALSE,
      label.fontsize = boxlabsize,
      con.type = 'straight',
      tol = 0.001,
      label.buffer = unit(0, 'mm')
    ) +
      guides(color = FALSE)
  }

  ## Plot formatting:
  p <- p +
    xlim(xmin, xmax) +
    ylim(ymin, ymax) +
    labs(x = lab.x, y = lab.y) +
    theme_bw() +
    theme(legend.position = legpos,
          legend.text = element_text(size = 12),
          legend.spacing.x = unit(0, 'cm'),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"),
          legend.title = element_text(size = 14, face = 'bold'),
          legend.background = element_rect(fill = "grey90", colour = "grey30"),
          legend.key = element_rect(fill = "grey90"),
          legend.key.size = unit(1, "cm"), #unit(1.3, "cm"),
          legend.text.align = 0,
          legend.box = legstacking,
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.length = unit(.25, "cm"),
          axis.title.x = element_text(size = axis_title_size,
                                      margin = unit(c(4, 0, 0, 0), 'mm')),
          axis.title.y = element_text(size = axis_title_size,
                                      margin = unit(c(0, 4, 0, 0), 'mm')),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20),
          plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm")) +
    guides(fill = guide_legend(override.aes = list(shape = 21)),
           color = guide_legend(order = 1))

  if(!is.null(plot_title)) p <- p + ggtitle(plot_title)

  if(plot_save == TRUE) {
    figfile <- paste0(my_outdir_figs, '/', pca_id, '.png')
    cat('## Saving plot to file:', figfile, '\n')
    ggsave(figfile, width = 8, height = 6)
    system(paste('xdg-open', figfile))
  }

  return(p)
}

