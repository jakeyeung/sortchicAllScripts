# Jake Yeung
# Date of Creation: 2022-03-10
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/18-update_public_data_add_bone_marrow.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# Add bone marrow to public data  -----------------------------------------


# show public data again 


# Load inputs -------------------------------------------------------------

indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/public_data"
inf.all <- file.path(indir, "Bartosovic_reanalysis_many_datasets_frip_with_Zeller_more_scCT.2021-06-25.rds")
dat.all <- readRDS(inf.all)


# add Ku  -----------------------------------------------------------------

inf.ku <- file.path(indir, "total_cuts_peak_cuts_Ku_2021.2021-09-05.txt")
dat.ku <- fread(inf.ku)


ggplot(dat.ku, aes(x = interaction(mark, dataset), y = cuts_total)) +
  geom_boxplot() + 
  scale_y_log10() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.ku.renamed <- dat.ku %>%
  dplyr::select(cell, cuts_total, cuts_in_peaks, mark, dataset) %>%
  dplyr::rename(allcounts = cuts_total, 
                jset = mark, 
                paper = dataset, 
                peakcounts = cuts_in_peaks) %>%
  rowwise() %>%
  mutate(jset = paste0("WBC_", jset))

dat.merged <- rbind(dat.all, dat.ku.renamed) 

dat.merged.sum <- dat.merged %>%
  group_by(paper, jset) %>%
  summarise(ncells = length(cell))


cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")


ggplot(subset(dat.merged, jset != "WBC_H3K4me3"), aes(x = jset, y = allcounts, fill = paper)) +
  geom_boxplot() + 
  scale_y_log10() + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  ylab("Unique counts") + 
  theme(aspect.ratio=1.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), legend.position = "bottom")

ggplot(subset(dat.merged, jset != "WBC_H3K4me3"), aes(x = jset, y = peakcounts / allcounts, fill = paper)) +
  geom_boxplot() + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  xlab("") + 
  ylab("FRiP") + 
  theme(aspect.ratio=1.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), legend.position = "bottom")


# Add zeller --------------------------------------------------------------


