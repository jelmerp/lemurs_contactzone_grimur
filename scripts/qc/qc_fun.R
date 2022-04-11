#### PROCESS SAMTOOLS FLAGSTAT OUTPUT ------------------------------------------
## Process a single file:
read_flagstats <- function(ID, suffix, bamdir) {

  infile_flagstats <- paste0(bamdir, "/", ID, suffix)

  flagstats <- readLines(infile_flagstats)

  n_tot <- as.integer(unlist(strsplit(flagstats[grep("in total", flagstats)], split = " "))[1])
  n_map <- as.integer(unlist(strsplit(flagstats[grep("mapped.*%", flagstats)], split = " "))[1])
  n_pr <- as.integer(unlist(strsplit(flagstats[grep("paired in sequencing", flagstats)], split = " "))[1])
  n_pp <- as.integer(unlist(strsplit(flagstats[grep("properly paired", flagstats)], split = " "))[1])

  return(c(ID, n_tot, n_map, n_pr, n_pp))
}

## Wrapper to process multiple files:
read_flagstats_all <- function(IDs,
                               bamdir = "qc/bam/",
                               suffix = ".samtools-flagstat.txt") {
  flagstats <- data.frame(
    do.call(rbind,
            lapply(IDs, read_flagstats, bamdir = bamdir, suffix = suffix)),
    stringsAsFactors = FALSE
  )
  colnames(flagstats) <- c("ID", "bam_tot", "bam_map", "bam_pr", "bam_pp")
  integer_cols <- c("bam_tot", "bam_map", "bam_pr", "bam_pp")
  flagstats[, integer_cols] <- lapply(flagstats[, integer_cols], as.integer)
  return(flagstats)
}

#### PROCESS FASTQ STATS -------------------------------------------------------
## Process a single file:
read_fqstats <- function(ID, read = "R1", dir_fq) {
  #ID=mgri068_r03_p5a03; read = R1

  infile_fq <- paste0(dir_fq, "/", ID, ".", read, ".fastqstats.txt")
  fqstats <- read.delim(infile_fq, header = FALSE, as.is = TRUE) %>%
    filter(V1 %in% c("reads", "len mean", "%dup", "qual mean")) %>%
    pull(V2)
  return(c(ID, read, fqstats))
}

## Wrapper to process multiple files:
read_fqstats_all <- function(IDs, read = "R1", dir_fq) {
  fqstats <- data.frame(
    do.call(rbind, lapply(IDs, read_fqstats, read = read, dir_fq)),
    stringsAsFactors = FALSE
    )

  colnames(fqstats) <- c("ID", "read", "fq_reads", "fq_meanlen",
                            "fq_pct_dup", "fq_qual")

  fqstats$fq_reads <- as.integer(fqstats$fq_reads)
  fqstats$fq_meanlen <- round(as.numeric(fqstats$fq_meanlen), 2)
  fqstats$fq_pct_dup <- round(as.numeric(fqstats$fq_pct_dup), 2)
  fqstats$fq_qual <- round(as.numeric(fqstats$fq_qual), 4)

  fqstats <- fqstats %>%
    select(ID, read, fq_reads, fq_qual, fq_pct_dup, fq_meanlen)

  return(fqstats)
}

#### PROCESS EA-UTILS BAM STATS ------------------------------------------------
## Process a single file:
read_bamstats <- function(ID, bamdir, mycolnames, suffix) {

  infile_bamstats <- paste0(bamdir, "/", ID, suffix)

  bamstats <- read.delim(infile_bamstats, header = FALSE, as.is = TRUE, nrows = 42)
  bamstats <- bamstats[1:which(bamstats$V1 == "%N"), ]
  bamstats_df <- setNames(data.frame(matrix(ncol = nrow(bamstats), nrow = 1)), bamstats$V1)
  bamstats_df[1, ] <- bamstats$V2
  colnames(bamstats_df) <- colnames(bamstats_df) %>% gsub(" ", "_", .) %>% gsub("%", "pct", .)

  if(ncol(bamstats_df) < 42) {
    complete_df <- setNames(data.frame(matrix(ncol = 42, nrow = 0)), mycolnames)

    bamstats_df <-
      merge(complete_df, bamstats_df,
            by = colnames(bamstats_df), all.y = TRUE) %>%
      select(any_of(mycolnames))
  }

  bamstats_df$ID <- ID
  bamstats_df <- bamstats_df %>% select(ID, mycolnames)

  return(bamstats_df)
}

## Wrapper to process multiple files:
read_bamstats_all <- function(IDs, bamdir, mycolnames,
                              suffix = ".ea-utils-samstats.txt") {
  bamstats <- as.data.frame(do.call(
    rbind,
    lapply(IDs, read_bamstats, bamdir = bamdir, mycolnames = mycolnames,
           suffix = suffix))
    )
  bamstats$reads <- as.integer(bamstats$reads)
  bamstats$pct_ambiguous <- as.numeric(bamstats$pct_ambiguous)
  bamstats$mapq_mean <- as.numeric(bamstats$mapq_mean)

  return(bamstats)
}


#### PLOTS ------------------------------------------------
## Create a boxplot comparing stats across species:
ggbox <- function(
  my_df, my_x, my_y,
  colby = NULL, colby_cols = NULL, colby_labs = NULL, colby_name = NULL,
  shapeby = NULL, shapeby_name = NULL,
  xlabs = NULL, ymax = NULL,
  ytitle = "", xtitle = NULL, ptitle = NULL,
  saveplot = TRUE, openfig = FALSE,
  figdir = "results/qc/figs/", figname = NULL
  ) {

  p <- ggplot(data = my_df) +
    geom_boxplot(aes_string(x = my_x, y = my_y),
                 outlier.color = NA) +
    geom_jitter(aes_string(x = my_x, y = my_y, color = colby, shape = shapeby),
                width = 0.2, size = 2.5) +
    scale_color_manual(values = colby_cols,
                       labels = colby_labs,
                       name = colby_name) +
    scale_shape_discrete(name = shapeby_name) +
    scale_y_continuous(expand = c(0.015, 0.015)) +
    labs(y = ytitle) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size = 12, margin = margin(0.15, 0, 0.3, 0, "cm")),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14, margin = margin(0, 0.3, 0, 0, "cm")),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      plot.title = element_text(hjust = 0.5, size = 18),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text.align = 0
    )

  if(!is.null(ymax)) p <- p + scale_y_continuous(limits = c(NA, ymax),
                                                 expand = c(0.01, 0.01))
  if(!is.null(xlabs)) p <- p + scale_x_discrete(labels = xlabs)
  if(!is.null(xtitle)) p <- p + labs(y = ytitle, x = xtitle)
  if(!is.null(ptitle)) p <- p + ggtitle(ptitle)

  if(saveplot == TRUE) {
    if(is.null(figname)) figfile <- paste0(figdir, "/", my_y, ".png")
    if(!is.null(figname)) figfile <- paste0(figdir, "/", figfile, ".png")
    ggsave(figfile, p, width = 6, height = 6)
    if(openfig == TRUE) system(paste("xdg-open", figfile))
  }
  return(p)
}
