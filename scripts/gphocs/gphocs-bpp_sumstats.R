# SET-UP -----------------------------------------------------------------------

## Packages
library(here)

## Input files
infile_logs <- here('results/gphocs/output/gphocs-bpp_mergedlogs.txt')

## Load GPHOCS and BPP logs
logs <- as.data.frame(fread(infile_logs, stringsAsFactors = FALSE))


# SUMMARIZE GPHOCS AND BPP RESULTS ---------------------------------------------

## Divergence times
(tau.sum <- filter(logs, var == 'tau') %>%
    group_by(runID, pop) %>%
    summarise(tau = round(mean(cval) / 1000, 1),
              min = round(hpd.min(cval) / 1000, 1),
              max = round(hpd.max(cval) / 1000, 1)) %>%
    arrange(pop))

## Population sizes
(th.sum <- filter(logs, var == 'theta') %>%
    group_by(runID, pop) %>%
    summarise(Ne = round(mean(cval) / 1000),
              min = round(hpd.min(cval) / 1000),
              max = round(hpd.max(cval) / 1000)) %>%
    arrange(pop))

th.sum %>% filter(runID == "g2m2anc") %>% print(n=100)

## Migration rates
(m.sum <- filter(logs, var == 'mprop', runID == 'g2m2anc') %>%
        group_by(migfrom, migto, var) %>%
        summarise(m.prop = round(mean(val), 4)))