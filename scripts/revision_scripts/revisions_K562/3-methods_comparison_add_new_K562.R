# Jake Yeung
# Date of Creation: 2022-03-23
# File: ~/projects/scchic/scripts/revision_scripts/revisions_K562/3-methods_comparison_add_new_K562.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

# 
# outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_cluster/K562_new/plots"
# outpdf <- file.path(outdir, paste0("total_cuts_frip_with_Ku_2021_new_K562.", Sys.Date(), ".pdf"))
# pdf(outpdf, useDingbats = FALSE)

# Load inputs -------------------------------------------------------------

indir.public <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/public_data"
inf.all <- file.path(indir.public, "Bartosovic_reanalysis_many_datasets_frip_with_Zeller_more_scCT.2021-06-25.rds")
dat.all <- readRDS(inf.all)



# add Ku  -----------------------------------------------------------------

inf.ku <- file.path(indir.public, "total_cuts_peak_cuts_Ku_2021.2021-09-05.txt")
dat.ku <- fread(inf.ku)

pdf(paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_cluster/K562_new/plots/compare_with_other_methods_with_new_K562.", Sys.Date(), ".pdf"), useDingbats = FALSE)

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

ggplot(dat.merged, aes(x = jset, y = allcounts, fill = paper)) +
  geom_boxplot() + 
  scale_y_log10() + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  ylab("Unique counts") + 
  theme(aspect.ratio=1.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), legend.position = "bottom")


ggplot(dat.merged, aes(x = jset, y = peakcounts / allcounts, fill = paper)) +
  geom_boxplot() + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  xlab("") + 
  ylab("FRiP") + 
  theme(aspect.ratio=1.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), legend.position = "bottom")

ggplot(dat.merged.sum, aes(x = jset, y = ncells, color = paper)) +
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  xlab("") + 
  ylab("Number of Cells") + 
  ggtitle(paste(dat.merged.sum$jset, dat.merged.sum$ncells, collapse = "\n")) + 
  theme(aspect.ratio=1.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), legend.position = "bottom")

# remove H3K4me3 Ku so everything is K27me3

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

ggplot(dat.merged.sum, aes(x = jset, y = ncells, color = paper)) +
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  xlab("") + 
  ylab("Number of Cells") + 
  ggtitle(paste(dat.merged.sum$jset, dat.merged.sum$ncells, collapse = "\n")) + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), legend.position = "bottom")


# Add new K562  ----------------------------------------------------------

inrds <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_cluster/K562_new/plots/qc_plots_K562_with_frip.2022-03-23.rds"
assertthat::assert_that(file.exists(inrds))

dat.meta.new.lst <- readRDS(inrds) 

dat.meta.new.mark <- dat.meta.new.lst[["k27me3"]] %>%
  filter(is.good) %>%
  dplyr::select(samp, total_cuts, total_cuts_in_peaks) %>%
  dplyr::rename(cell = samp, 
                allcounts = total_cuts, 
                peakcounts = total_cuts_in_peaks)  %>%
  dplyr::mutate(jset = "ZellerNew", 
                paper = "Zeller")

  

dat.merged.k27me3.with_new <- subset(dat.merged, jset != "WBC_H3K4me3") %>%
  bind_rows(., dat.meta.new.mark)

dat.merged.k27me3.with_new.sum <- dat.merged.k27me3.with_new %>%
  group_by(paper, jset) %>%
  summarise(ncells = length(cell))

ggplot(dat.merged.k27me3.with_new, aes(x = jset, y = peakcounts / allcounts, fill = paper)) +
  geom_boxplot() + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  xlab("") + 
  ylab("FRiP") + 
  theme(aspect.ratio=1.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), legend.position = "bottom")

ggplot(dat.merged.k27me3.with_new, aes(x = jset, y = log10(allcounts), fill = paper)) +
  geom_boxplot() + 
  scale_y_log10() + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  ylab("Unique counts") + 
  scale_y_continuous(breaks = seq(6), labels = seq(6)) +
  theme(aspect.ratio=1.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), legend.position = "bottom")

# without Zeller

ggplot(subset(dat.merged.k27me3.with_new, jset != "Zeller"), aes(x = jset, y = peakcounts / allcounts, fill = paper)) +
  geom_boxplot() + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  xlab("") + 
  ylab("FRiP") + 
  theme(aspect.ratio=1.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), legend.position = "bottom")

ggplot(subset(dat.merged.k27me3.with_new, jset != "Zeller"), aes(x = jset, y = log10(allcounts), fill = paper)) +
  geom_boxplot() + 
  scale_y_log10() + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  ylab("Unique counts") + 
  theme(aspect.ratio=1.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), legend.position = "bottom")

ggplot(subset(dat.merged.k27me3.with_new, jset != "Zeller"), aes(x = jset, y = log10(allcounts), fill = paper)) +
  geom_boxplot() + 
  scale_y_log10() + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  ylab("Unique counts") + 
  scale_y_continuous(breaks = seq(6), labels = seq(6)) +
  theme(aspect.ratio=1.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), legend.position = "bottom")

ggplot(dat.merged.k27me3.with_new.sum, aes(x = jset, y = ncells, color = paper)) +
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  xlab("") + 
  ylab("Number of Cells") + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), legend.position = "bottom")


ggplot(dat.merged.k27me3.with_new.sum, aes(x = jset, y = ncells, color = paper)) +
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  xlab("") + 
  ylab("Number of Cells") + 
  ggtitle(paste(dat.merged.k27me3.with_new.sum$jset, dat.merged.k27me3.with_new.sum$ncells, collapse = "\n")) + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), legend.position = "bottom")


dev.off()




