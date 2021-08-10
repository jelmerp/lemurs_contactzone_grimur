library(tidyverse)

#### FUNCTIONS -----------------------------------------------------------------
preplogs <- function(logdir, bpp_prefix,
                      setID, runID, rep,
                      lookup, t_scale,
                      burnin, subsample_by,
                      gentime, gentime_sd,
                      mutrate_gen, mutrate_var) {

  logfiles <- list.files(path = logdir, pattern = bpp_prefix, full.names = TRUE)
  Log <- NULL

  for(logfile in logfiles) {
    rep <- gsub('//', '/', logfile) %>%
      gsub(logdir, '', .) %>%
      gsub(bpp_prefix, '', .) %>%
      gsub('.txt', '', .) %>%
      gsub('r', '', .) %>%
      as.integer()
    cat("Rep:", rep)

    onelog <-
      read.table(logfile, header = TRUE, as.is = TRUE) %>%
      preplog(., setID = setID, runID = runID, rep = rep,
               lookup = lookup, t_scale = t_scale,
               burnin = burnin, subsample_by = subsample_by,
               gentime = gentime, gentime_sd = gentime_sd,
               mutrate_gen = mutrate_gen, mutrate_var = mutrate_var)

    if(is.null(Log)) Log <- onelog else Log <- rbind(Log, onelog)

    cat(" nrow-single:", nrow(onelog), " nrow-total:", nrow(Log), '\n')
  }

  return(Log)
}

## prep Log:
preplog <- function(Log,
                    runID, setID, rep,
                    lookup, t_scale,
                    burnin = 100000, subsample_by = 50,
                    gentime, gentime_sd,
                    mutrate_gen, mutrate_var) {

  Log <- Log %>%
    rm_burnin(., burnin) %>%
    subsample_log(., subsample_by) %>%
    select(-lnL) %>%
    pivot_longer(-Gen, names_to = 'var_pop') %>%
    separate(var_pop, sep = '_', into = c('var', 'pop'), extra = 'merge') %>%
    rename(val = value, Sample = Gen) %>%
    mutate(setID = setID,
           runID = runID,
           rep = rep,
           pop = as.character(poprename(pop, lookup)),
           migpattern = ifelse(var == 'phi', pop, NA),
           pop = ifelse(var == 'phi', NA, pop),
           migfrom = gsub('(.*)_2_.*', '\\1', migpattern),
           migto = gsub('.*_2_(.*)', '\\1', migpattern),
           cval = NA,
           cn = 'aa') %>%
    select(setID, runID, rep, Sample,
           pop, migpattern, migfrom, migto,
           cn, var, val, cval)

  ## Add converted values:
  gentime_dist <- get_gentime_dist(nrow(Log), gentime, gentime_sd)
  mutrate_dist <- get_mutrate_dist(nrow(Log), mutrate_gen, mutrate_var)
  Log <- Log %>%
    convert_theta(., t_scale, mutrate_dist) %>%
    convert_tau(., t_scale, mutrate_dist, gentime_dist)

  return(Log)
}

## subsample:
subsample_log <- function(Log, subsample_by) {
  Log <- Log[seq(from = 1, to = nrow(Log), by = subsample_by), ]
}


## remove burnin:
rm_burnin <- function(Log, burnin) {
  if(burnin > 0) {
    Log <- Log %>% filter(Gen > burnin)
    return(Log)
  } else {
    cat("## rm_burnin(): Not removing burnin...")
  }
  if(nrow(Log) < 1000) {
    warning('\n## rm_burnin(): SKIPPING FILE: <1K ROWS LEFT AFTER REMOVING BURN-IN\n')
    return(NULL)
  }
}

## gentime_dist:
get_gentime_dist <- function(n, gentime, gentime_sd) {
  gentime_dist <- rlnorm(n, meanlog = log(gentime), sdlog = log(gentime_sd))
  return(gentime_dist)
}

## mutrate_dist:
get_mutrate_dist <- function(n, mutrate_gen, mutrate_var,
                             mutrate_scale = 1e-8) {
  my_shape <- mutrate_gen^2 / mutrate_var
  my_rate <- mutrate_gen / mutrate_var
  mutrate_dist <- rgamma(n, shape = my_shape, rate = my_rate) * mutrate_scale
  return(mutrate_dist)
}

## Add converted demographic values:
convert_theta <- function(Log, t_scale, mutrate_gen_dist) {
  my_theta  <- Log$val[Log$var == 'theta'] * t_scale
  my_mutrate_gen <- mutrate_gen_dist[Log$var == 'theta']
  Log$cval[Log$var == 'theta'] <- my_theta / (4 * my_mutrate_gen)
  return(Log)
}

convert_tau <- function(Log, t_scale, mutrate_gen_dist, gentime_dist) {
  my_tau <- Log$val[Log$var == 'tau'] * t_scale
  my_mutrate_yr <- mutrate_gen_dist[Log$var == 'tau'] / gentime_dist[Log$var == 'tau']
  Log$cval[Log$var == 'tau'] <- my_tau / my_mutrate_yr
  return(Log)
}

## Rename pops:
poprename <- function(pop, lookup) {
  pop <- gsub('[0-9]+', '', pop)
  lookup$val[match(pop, lookup$key)]
}
