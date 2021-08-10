#### PREP DATAFRAME FOR DEMOGRAPHY PLOT ----------------------------------------
ttprep <- function(logs = NULL, popcols, poporder,
                   tt = NULL, summary.provided = FALSE, popnames = NULL,
                   pops.fixed = NULL, pops.fixed.size = 1,
                   x.even = FALSE, pop.spacing = 50, x.start = 5,
                   extra.time.root = 50) {

  if(summary.provided == FALSE) {
    tt <- logs %>%
      subset(var %in% c('theta', 'tau')) %>%
      group_by(pop, var) %>%
      summarise(cval.mean = mean(cval / 1000)) %>%
      pivot_wider(names_from = var, values_from = cval.mean,
                  values_fill = list(cval.mean = 0)) %>%
      ungroup() %>%
      mutate(pop = as.character(pop))
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
  tt$x.min[getpos('A')] <- get.xmax('murW')
  tt$x.min[getpos('B')] <- get.xmax('A')
  tt$x.min[getpos('murE')] <- tt$x.min[getpos('B')] - tt$theta[getpos('murE')]
  tt$x.min[getpos('murC')] <- get.xmax('B')
  tt$x.min[getpos('griC')] <- get.xmax('murC') + pop.spacing
  tt$x.min[getpos('C')] <- get.xmax('griC')
  tt$x.min[getpos('griW')] <- get.xmax('C')
  tt$x.min[getpos('root')] <- get.xmax('B')
  Diff2 <- ((get.xmin('C') - get.xmax('A')) / 2) - (get.th('root') / 2)
  tt$x.min[getpos('root')] <- get.xmax('A') + Diff2

  ## x end positions:
  tt$x.max <- round(tt$x.min + tt$theta, 2)

  ## y positions:
  tt$y.min <- round(ifelse(tt$pop %in% currentpops, 0, tt$tau), 2)
  tt$y.max <- round(ifelse(tt$pop == 'root', tt$y.min + extra.time.root,
                           get.ta(getparent(tt$pop))), 2)

  ## Popcols and popnames:
  tt$popcol <- popcols$col[match(tt$pop, popcols$pop)]
  if(all(is.na(tt$popcol))) tt$popcol <- 'grey30'

  tt$popalpha <- popcols$alpha[match(tt$pop, popcols$pop)]
  if(all(is.na(tt$popalpha))) tt$popcol <- 1

  if(!is.null(popnames)) tt$pop <- popnames

  ## Finalize:
  tt <- tt %>%
    mutate(pop = factor(pop, levels = poporder)) %>%
    arrange(pop)

  return(tt)
}


#### DEMOGRAPHY PLOT WRAPPER FOR EASTWEST RUNS ---------------------------------

dplotwrap <- function(
  logs, runID_f, poplist, popcols, poporder,
  extra.time.root = 50, pop.spacing = 50,
  pops_to_conn = 'root',
  ann.pops.anc = TRUE, popnames.adj.vert = 40,
  x.min = 0, yticks.by = 200,
  ...)
  {

  ## Prepare df underlying plot:
  ttp <- ttprep(filter(logs, runID == runID_f),
                popcols = popcols, poporder = poporder,
                extra.time.root = extra.time.root, pop.spacing = pop.spacing)

  ## Main plot:
  p <- dplot(tt = ttp,
             pops_to_conn = pops_to_conn,
             ann.pops.anc = ann.pops.anc, popnames.adj.vert = popnames.adj.vert,
             x.min = x.min, yticks.by = yticks.by, ...)

  return(p)
}

getparent <- function(child) as.character(poplist[child])