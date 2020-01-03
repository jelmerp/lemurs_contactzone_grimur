################################################################################
##### PREP DATAFRAME FOR DEMOGRAPHY PLOT - SNAPP12 RUNS #####
################################################################################
ttPrepMur3 <- function(Log = NULL, tt = NULL, summary.provided = FALSE,
                       pops.fixed = NULL, pops.fixed.size = 1, popnames = NULL,
                       x.even = FALSE, pop.spacing = 5, x.start = 5,
                       extra.time.root = 100) {

  if(summary.provided == FALSE) {
    tt <- Log %>%
      subset(var %in% c('theta', 'tau')) %>%
      group_by(pop, var) %>%
      dplyr::summarise(cval.mean = mean(cval / 1000)) %>%
      dcast(pop ~ var)
  }

  tt$NeToScale <- 1
  if(!is.null(pops.fixed)) tt$theta[match(pops.fixed, tt$pop)] <- pops.fixed.size

  ## Functions:
  getpos <- function(pop) match(pop, tt$pop)
  get.xmin <- function(pop) tt$x.min[match(pop, tt$pop)]
  get.th <- function(pop) tt$theta[match(pop, tt$pop)]
  get.ta <- function(pop) tt$tau[match(pop, tt$pop)]
  get.xmax <- function(pop) round(get.xmin(pop) + get.th(pop), 2)

  ## x start positions:
  tt$x.min <- NA
  tt$x.min[getpos('murW')] <- x.start
  tt$x.min[getpos('a.mur')] <- get.xmax('murW')
  tt$x.min[getpos('a.murSE')] <- get.xmax('a.mur')
  tt$x.min[getpos('murGan')] <- tt$x.min[getpos('a.murSE')] - tt$theta[getpos('murGan')]
  tt$x.min[getpos('murHZ')] <- get.xmax('a.murSE')
  tt$x.min[getpos('griHZ')] <- get.xmax('murHZ') + pop.spacing
  tt$x.min[getpos('a.gri')] <- get.xmax('griHZ')
  tt$x.min[getpos('griSW')] <- get.xmax('a.gri')
  tt$x.min[getpos('a.root')] <- get.xmax('a.mur')
  Diff2 <- ((get.xmin('a.gri') - get.xmax('a.mur')) / 2) - (get.th('a.root') / 2)
  tt$x.min[getpos('a.root')] <- get.xmax('a.mur') + Diff2

  ## x end positions:
  tt$x.max <- round(tt$x.min + tt$theta, 2)

  ## y positions:
  tt$y.min <- round(ifelse(tt$pop %in% currentpops, 0, tt$tau), 2)
  tt$y.max <- round(ifelse(tt$pop == 'a.root', tt$y.min + extra.time.root,
                           get.ta(getparent(tt$pop))), 2)

  ## Popcols and popnames:
  tt$popcol <- popcols.df$popcol[match(tt$pop, popcols.df$popname.short)]
  if(all(is.na(tt$popcol))) tt$popcol <- 'grey30'

  if(!is.null(popnames)) tt$pop <- popnames

  return(tt)
}


################################################################################
##### DEMOGRAPHY PLOT WRAPPER FOR EASTWEST RUNS #####
################################################################################
dplotwrap.mur3 <- function(runID.focal) {
  #runID.focal <- 'g2m2anc'

  ## Prepare df underlying plot:
  ttp <- ttPrepMur3(subset(Log, runID == runID.focal),
                    x.even = FALSE, pop.spacing = 20, extra.time.root = 25)

  ## Factor ordering for correct legend:
  ttp$pop <- factor(ttp$pop, levels = levels(Log$pop))
  ttp <- arrange(ttp, pop)
  #ttp$popcol <- factor(ttp$popcol, levels = ttp$popcol)

  ## Main plot:
  p <- dplot(tt = ttp, ann.pops = FALSE, x.min = 0, yticks.by = 100,
             popnames.adj.horz = rep(0, nrow(ttp)), popnames.adj.vert = 15,
             popnames.col = popnames.col, popnames.size = 5, x.extra = 5,
             saveplot = FALSE, plot.title = '') +
    geom_segment(aes(y = ttp$y.max[ttp$pop == 'a.mur'],
                     yend = ttp$y.max[ttp$pop == 'a.mur'],
                     x = ttp$x.max[ttp$pop == 'a.mur'],
                     xend = ttp$x.min[ttp$pop == 'a.gri']), colour = 'grey50')

  ## Save plot:
  plotfile.png <- paste0(plotdir, '/demo/', setID, '.', runID.focal, '.demoplot.png')
  ggsave(filename = plotfile.png, plot = p, width = 8, height = 7)
  system(paste("xdg-open", plotfile.png))

  plotfile.pdf <- paste0(plotdir, '/demo/', setID, '.', runID.focal, '.demoplot.pdf')
  ggsave(filename = plotfile.pdf, plot = p, width = 8, height = 7)

  print(p)
  return(p)
}