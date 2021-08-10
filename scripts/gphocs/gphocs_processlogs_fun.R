## Packages:
if(!'pacman' %in% rownames(installed.packages())) install.packages('pacman')
packages <- c('data.table', 'TeachingDemos', 'tidyverse')
pacman::p_load(char = packages, install = TRUE)


#### Step 1: cutlog(): Open one logfile and cut off burn-in and/or last samples ---------
cut_log <- function(
  logfile, logdir,
  burnin, last_sample, subsample,
  return_log = TRUE, write_log = TRUE
  ) {

  Log <- read.table(paste0(logdir, '/', logfile), header = TRUE)
  cat('\n\n## cut_log():', logfile, "\tnrows:", nrow(Log), '\n')

  ## Remove burn-in:
  if(burnin > 0) Log <- Log[-which(Log$Sample < burnin), ]
  if(nrow(Log) < 1000) {
    warning('\n## cut_log(): SKIPPING FILE: <1K ROWS LEFT AFTER REMOVING BURN-IN\n')
    return(NULL)
  }

  ## Remove final samples if "last_sample" arg is specified:
  if(!is.null(last_sample) && any(Log$Sample > last_sample)) {
    Log <- Log[-which(Log$Sample > last_sample), ]
  }

  ## Subsample from the Log::
  Log <- Log[seq(from = 1, to = nrow(Log), by = subsample), ]
  cat("## cut_log(): Last sample:", max(Log$Sample), '\n')

  if(write_log == TRUE) {
    if(!dir.exists(paste0(logdir, '/cut/'))) dir.create(paste0(logdir, '/cut/'))
    write.table(Log, paste0(logdir, '/cut/', logfile),
                sep = '\t', quote = FALSE, row.names = FALSE)
  }

  if(return_log == TRUE) return(Log)
}


#### Step 2: preplog(): Melt dataframe and prep vars --------------------------------------
preplog <- function(
  Log, logfile, setID, poplist, lookup,
  mutrate_var, gentime_sd,
  gdi = TRUE,
  rename_pops_before = FALSE, # Rename pops (using "lookup" df) before converting to long format, i.e. in column names
  rename_pops_after = FALSE # Rename pops after converting to long format, i.e. in df itself
  ) {

  #logfile <- 'hz.mur3gri2c_g2m2_190301_rep1.20190301-1204.log'
  #Log <- read.table(paste0(logdir, '/', logfile), header = TRUE)

  ## Rename columns:
  colnames(Log) <- gsub('\\.\\.', 2, colnames(Log))
  colnames(Log)[(ncol(Log)-2):ncol(Log)] <- c('mut_NA', 'dataLd_NA', 'FullLd_NA')

  ## Rename pops:
  if(rename_pops_before == TRUE) {
    varcols_tt <- grep('theta_|tau_', colnames(Log))
    for(col_idx in varcols_tt) for(pop_idx in 1:nrow(lookup)) {
      pop <- gsub('theta_|tau_|m_(.*)', '\\1', colnames(Log)[col_idx])
      long <- lookup$popname.long[pop_idx]
      short <- lookup$popname.short[pop_idx]
      if(pop == long) colnames(Log)[col_idx] <- sub(long, short, colnames(Log)[col_idx])
    }
    varcols_m <- grep('^m_', colnames(Log))
    for(col_idx in varcols_m) for(pop_idx in 1:nrow(lookup)) {
      pop1 <- gsub('m_(.*)2.*', '\\1', colnames(Log)[col_idx])
      pop2 <- gsub('m_.*2(.*)', '\\1', colnames(Log)[col_idx])
      long <- lookup$popname.long[pop_idx]
      short <- lookup$popname.short[pop_idx]
      if(pop1 == long) colnames(Log)[col_idx] <- sub(long, short, colnames(Log)[col_idx])
      if(pop2 == long) colnames(Log)[col_idx] <- sub(long, short, colnames(Log)[col_idx])
    }

    cat('## cutpreplog(): Renamed colnames:\n', colnames(Log), '\n')
  }

  ## Get migration:
  migcols <- grep('m_', colnames(Log))
  for(migcol in migcols) {
    m <- Log[, migcol]

    migpattern <- unlist(strsplit(colnames(Log)[migcol], split = '_'))[2]
    migto <- unlist(strsplit(migpattern, split = '2'))[2]
    migto.column <- grep(paste0('theta_', migto, '$'), colnames(Log))
    th <- Log[, migto.column]

    ## Population migration rate:
    colname.mig2 <- paste0('2Nm_', migpattern)
    Log$newcolumn1 <- (m * m_scale) * (th * t_scale) / 4
    colnames(Log)[grep('newcolumn1', colnames(Log))] <- colname.mig2

    ## Proportion of migrants:
    my_shape <- mutrate_gen^2 / mutrate_var
    my_rate <- mutrate_gen / mutrate_var
    mutrate_gen_dist <- rgamma(nrow(Log), shape = my_shape, rate = my_rate) * 1e-8

    colname.mig3 <- paste0('mprop_', migpattern)
    Log$newcolumn2 <- (m * m_scale) * mutrate_gen_dist * 100
    colnames(Log)[grep('newcolumn2', colnames(Log))] <- colname.mig3
  }

  ## Get gdi:
  if(gdi == TRUE) {
    childpops <- names(poplist)
    parentpops <- as.character(poplist)

    get_parent <- function(pop, level = 1) {
      if(level == 1) parentpop <- poplist[[pop]]
      if(level == 2) parentpop <- poplist[[poplist[[pop]]]]
      return(parentpop)
    }

    get_gdi <- function(pop, level) {
      parentpop <- get_parent(pop, level = level)
      theta <- Log[, paste0('theta_', pop)]
      tau <- Log[, paste0('tau_', parentpop)]
      gdi <- 1 - exp((-2 * tau) / theta)
    }

    for(pop in childpops) {
      Log$gdi_new <- get_gdi(pop, level = 1)
      colnames(Log)[grep('gdi_new', colnames(Log))] <- paste0('gdi_', pop)
    }

    for(pop in childpops) {
      if(!is.null(get_parent(pop, level = 2))) {
        Log$gdi_new <- get_gdi(pop, level = 2)
        colnames(Log)[grep('gdi_new', colnames(Log))] <- paste0('gdi2_', pop)
      }
    }
  }

  ## Get runID:
  runID <- logfile %>%
    gsub(paste0(setID, '_'), '', .) %>%
    strsplit(., split = '_') %>%
    unlist(.) %>%
    .[1]
  cat('## Run ID:', runID, '\n')

  ## Convert to long format:
  mlog <- Log %>%
    pivot_longer(-Sample, names_to = 'var_pop') %>%
    separate(var_pop, sep = '_', into = c('var', 'pop'), extra = 'merge') %>%
    separate(pop, sep = '2', into = c('migfrom', 'migto'),
             remove = FALSE, fill = 'right') %>%
    rename(val = value) %>%
    mutate(migfrom = ifelse(is.na(migto), NA, migfrom),
           migfrom = factor(migfrom),
           migto = factor(migto),
           setID = setID,
           runID = runID) %>%
    select(setID, runID, Sample, var, val, pop, migfrom, migto)

  mlog$migfrom[which(is.na(mlog$migto))] <- NA ## set migfrom to NA for non-migration vars
  mlog$pop[which(!is.na(mlog$migto))] <- NA ## set pop to NA for migration var (i.e. when migto is not NA)

  ## Add run rep(licate) ID:
  mlog$rep <- factor(as.integer(gsub('.*rep([0-9]).*log', '\\1', logfile)))

  ## Migration type in run (single migration band vs multiple vs none):
  mlog$migtype.run <- factor(
    ifelse(grepl('multmig', runID), 'mult',
           ifelse(grepl('noMig', runID, ignore.case = TRUE), 'none', 'single'))
  )

  ## Rename pops:
  if(rename_pops_after == TRUE) {
    if(rename_pops_before == FALSE) {
      warning('## preplog(): These pops were not found in the lookup:',
              unique(mlog$pop)[!unique(mlog$pop) %in% lookup$popname.long], '\n')
      mlog$pop <- poprename(mlog$pop, lookup)
    }
    mlog$migfrom <- poprename(mlog$migfrom, lookup)
    mlog$migto <- poprename(mlog$migto, lookup)
  }

  ## Add converted values:
  mlog <- add_cvalue(mlog, mutrate_var = mutrate_var,
                     gentime_sd = gentime_sd, poplist = poplist)

  return(mlog)
}


#### cutpreplog(): Single-log wrapper function: cut_log() then preplog() --------------------
cutpreplog <- function(
  logfile, logdir,
  setID, runID, poplist, lookup,
  burnin, last_sample, subsample, cut = TRUE,
  mutrate_var, gentime_sd,
  gdi = TRUE,
  rename_pops_before, rename_pops_after,
  return_log = TRUE, write_log = TRUE
  ) {

  if(cut == TRUE) {
    Log <- cut_log(
      logfile = logfile, logdir = logdir, burnin = burnin,
      last_sample = last_sample, subsample = subsample,
      return_log = return_log, write_log = write_log
      )
  }

  if(cut == FALSE) {
    cat('\n\n## cutpreplog(): Not cutting log, just reading it in.')
    cat('## cutpreplog(): ', logfile, "\tnrows:", nrow(Log))
    Log <- read.table(paste0(logdir, '/cut/', logfile), header = TRUE)
  }

  if(!is.null(Log)) {
    cat('## cutpreplog(): Column names of Log:\n', colnames(Log), '\n')

    Log <- preplog(
      Log = Log, logfile = logfile,
      setID = setID, poplist = poplist, lookup = lookup,
      rename_pops_before = rename_pops_before, rename_pops_after = rename_pops_after,
      mutrate_var = mutrate_var, gentime_sd = gentime_sd,
      gdi = gdi
      )
    return(Log)
  }
}


#### getlogs(): wrapper function to merge all log files into a single df -------
getlogs <- function(
  setID, logdir, poplist, lookup = NULL,
  cut = TRUE, burnin = 100000, last_sample = NULL, subsample = 50,
  rename_pops_before = FALSE, rename_pops_after = FALSE,
  mutrate_var = 0, gentime_sd = 0, gdi = FALSE
  ) {

  ## Collect all logs:
  logfiles <- list.files(logdir, pattern = paste0('.*', setID, '.*\\.log'))
  cat("## getlogs(): Found the following logfiles for", setID, ":\n")
  print(logfiles)

  Log <- lapply(logfiles, cutpreplog,
                rename_pops_before = rename_pops_before, rename_pops_after = rename_pops_after,
                mutrate_var = mutrate_var, gentime_sd = gentime_sd,
                logdir = logdir, setID = setID, burnin = burnin,
                last_sample = last_sample, subsample = subsample, cut = cut,
                gdi = gdi,
                poplist = poplist, lookup = lookup)
  Log <- compact(Log) # Remove NULL entries
  if(length(Log) > 1) Log <- do.call(rbind, Log) else Log <- Log[[1]]
  cat('### getlogs(): Done processing individuals logs.\n')

  ## Edit final merged logfiles:
  Log$setID <- factor(setID)
  Log$cn <- factor('aa')
  Log$migpattern <- paste0(Log$migfrom, '_2_', Log$migto)
  Log$migpattern <- gsub('NA_2_NA', NA, Log$migpattern)

  Log <- Log %>%
    select(setID, runID, rep, Sample, var, val, cval, pop,
           migfrom, migto, migpattern, migtype.run, cn)

  return(Log)
}

#### Helper functions ----------------------------------------------------------
## add_cvalue(): add converted demographic values
add_cvalue <- function(Log, gentime_sd, mutrate_var, poplist) {

  cat('## add_cvalue(): Adding converted values...\n')
  Log$cval <- NA

  ## Total migration rate:
  smr <<- Log %>%
    filter(var == 'tau') %>%
    group_by(pop) %>%
    summarise(tau = mean(val))

  focalpops <- as.character(unique(Log$migto))
  focalpops <- focalpops[!is.na(focalpops)]

  if(length(focalpops) > 1)
    newmig <- do.call(rbind, lapply(focalpops, getlifespan,
                                    Log = Log, smr = smr, poplist = poplist))
  if(length(focalpops) == 1)
    newmig <- getlifespan(focalpops, Log = Log, smr = smr, poplist = poplist)

  if(length(focalpops) >= 1)
    Log$cval[newmig$frows] <- (Log$val[newmig$frows] * m_scale) * (newmig$lifespan * t_scale)

  ## generation time and mut rate dist:
  gentime_dist <- rlnorm(nrow(Log), meanlog = log(gentime), sdlog = log(gentime_sd))

  my_shape <- mutrate_gen^2 / mutrate_var
  my_rate <- mutrate_gen / mutrate_var
  mutrate_gen_dist <- rgamma(nrow(Log), shape = my_shape, rate = my_rate) * 1e-8
  #mutrate_gen_dist <- rnorm(nrow(Log), mean = mutrate_gen, sd = mutrate_var) # old way using normal distribution

  ## theta:
  my_theta  <- Log$val[Log$var == 'theta'] * t_scale
  my_mutrate_gen <- mutrate_gen_dist[Log$var == 'theta']
  Log$cval[Log$var == 'theta'] <- my_theta / (4 * my_mutrate_gen)

  ## tau:
  my_tau <- Log$val[Log$var == 'tau'] * t_scale
  my_mutrate_yr <- mutrate_gen_dist[Log$var == 'tau'] / gentime_dist[Log$var == 'tau']
  Log$cval[Log$var == 'tau'] <- my_tau / my_mutrate_yr

  return(Log)
}

## getlifespan: Get lifespan for a pop in a specific run
getlifespan <- function(focalpop, Log, smr, poplist) {
  childpops <- names(poplist)
  parentpops <- as.character(poplist)
  currentpops <- setdiff(childpops, parentpops)

  if(is.null(focalpop)) focalpop <- unlist(strsplit(migRun, split = '2'))[2]
  tau.parent <- smr$tau[smr$pop == poplist[[focalpop]]]

  if(focalpop %in% currentpops) {
    lifespan <- tau.parent
  } else {
    lifespan <- tau.parent - smr$tau[smr$pop == focalpop]
  }

  frows <- which(Log$var == 'm' & Log$migto == focalpop)

  lifespan.df <- data.frame(frows, lifespan)
  return(lifespan.df)
}

## poprename(): Convert population name
poprename <- function(pop, lookup) {
  lookup$popname.short[match(pop, lookup$popname.long)]
}
