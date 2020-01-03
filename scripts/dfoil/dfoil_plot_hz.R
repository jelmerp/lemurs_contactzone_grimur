################################################################################
##### SET-UP #####
################################################################################
setwd('/home/jelmer/Dropbox/sc_lemurs/hybridzone/')
library(plyr); library(reshape2); library(tidyverse)
source('/home/jelmer/Dropbox/sc_lemurs/scripts/dfoil/dfoil_plot_fun.R')

###########################################################################
##### PROCESS RESULTS #####
###########################################################################
fileID.base <- 'r03.wOutgroups.indsel.mac3.FS6.'

## HZ: mgri_sw -- mgri_hz -- mmur_hz -- mmur_w
(dfoil <- read.dfoil.out(fileID = paste0(fileID.base, 'HZ'), id.short = 'HZ', alt = FALSE))
plot.dfoil(dfoil, fileID.suffix = 'HZ') # --++
(dfoil <- read.dfoil.out(fileID = paste0(fileID.base, 'HZ'), id.short = 'HZ', alt = TRUE))
plot.dfoil(dfoil, fileID.suffix = 'HZ.alt') # 0(+)(+)+

## SE: mgri_sw -- mgri_se --mmur_se -- mmur_w
(dfoil <- read.dfoil.out(fileID = paste0(fileID.base, 'SE'), id.short = 'SE', alt = FALSE))
plot.dfoil(dfoil, fileID.suffix = 'SE') # --++
(dfoil <- read.dfoil.out(fileID = paste0(fileID.base, 'SE'), id.short = 'SE', alt = TRUE))
plot.dfoil(dfoil, fileID.suffix = 'SE.alt') # 0+-+ # P3->P2 except for last +

## HZgan: mgri_sw -- mgri_hz -- mmur_hz -- mmur_gan
(dfoil <- read.dfoil.out(fileID = paste0(fileID.base, 'HZgan'), id.short = 'SE', alt = FALSE))
plot.dfoil(dfoil, fileID.suffix = 'HZgan') # 00++ / --++
(dfoil <- read.dfoil.out(fileID = paste0(fileID.base, 'HZgan'), id.short = 'SE', alt = TRUE))
plot.dfoil(dfoil, fileID.suffix = 'HZgan.alt')

## SEgan: mgri_sw -- mgri_se -- mmur_se -- mmur_gan
(dfoil <- read.dfoil.out(fileID = paste0(fileID.base, 'SEgan'), id.short = 'SE', alt = FALSE))
plot.dfoil(dfoil, fileID.suffix = 'SEgan') # 00++
(dfoil <- read.dfoil.out(fileID = paste0(fileID.base, 'SEgan'), id.short = 'SE', alt = TRUE))
plot.dfoil(dfoil, fileID.suffix = 'SEgan.alt') # +00-

## gan: mgri_sw -- mgri_hz -- mmur_gan -- mmur_w
(dfoil <- read.dfoil.out(fileID = paste0(fileID.base, 'gan'), id.short = 'SE', alt = FALSE))
plot.dfoil(dfoil, fileID.suffix = 'gan') # --++ (closest to P1->P4)
(dfoil <- read.dfoil.out(fileID = paste0(fileID.base, 'gan'), id.short = 'SE', alt = TRUE))
plot.dfoil(dfoil, fileID.suffix = 'gan.alt') # close to -0++ = P4->P1


## Allele pattern counts:
#dfoil.in <- read.dfoil.in(fileID = fileID, id.short = 't10.maf0.15')
#dfoil.in %>% select(BBABA, BBBAA)
#dfoil.in %>% select(ABABA, ABBAA)
#dfoil.in %>% select(BAABA, BABAA)
