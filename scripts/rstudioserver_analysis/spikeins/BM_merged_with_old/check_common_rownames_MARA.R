# Jake Yeung
# Date of Creation: 2020-11-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/check_common_rownames_MARA.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# H3K4me3 has lots of common names ----------------------------------------

jmark <- "H3K4me3"
# inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged3_Peaks/H3K4me3/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.H3K4me3.txt"
inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged3_Peaks/", jmark, "/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, ".txt")
dat <- fread(inf)

dat$chromo <- sapply(dat$Gene.ID, GetChromo)

dat.sum <- dat %>%
  group_by(chromo) %>%
  summarise(npeaks = length(Gene.ID))

ggplot(dat.sum, aes(x = chromo, y = npeaks)) + 
  geom_col() + 
  theme_bw() + 
  ggtitle(jmark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


inf.E <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged3_Peaks/count_mats_peaks_norm/ldaOut.count_mat_from_hiddendomains.", jmark, ".filtNAcells_topbins.K-30.keepNbins_0.txt")
E <- fread(inf.E)

common.rows <- intersect(dat$Gene.ID, E$Gene.ID)
notcommon.rows <- setdiff(E$Gene.ID, dat$Gene.ID)

print(length(common.rows))

hubprefix <- "/home/jyeung/hub_oudenaarden"

inf.nochromo <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds/hd_merged.H3K4me1.FromR.maxcount_40_60_R/merged.H3K4me1.cutoff_analysis.merged.nochr.bed")
inf.withchromo <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds/hd_merged.H3K4me1.FromR.maxcount_40_60_R/merged.H3K4me1.cutoff_analysis.merged.withchr2.bed")

bed.nochromo <- fread(inf.nochromo)
bed.withchromo <- fread(inf.withchromo)

# H3K4me1 not much  -------------------------------------------------------










