# Jake Yeung
# Date of Creation: 2021-04-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/dendrograms/10-annotate_GC_dist2TSS_DE_bins.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"


# Load bins ---------------------------------------------------------------

infs.highbins <- lapply(jmarks, function(jmark){
  file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/High_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.", jmark, ".2021-01-30.txt"))
})

infs.debins <- lapply(jmarks, function(jmark){
  file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/tables_top_6085_four_marks_dynamic_bins/top_6085_bins_nearest_gene_gc.", jmark, ".2021-02-17.txt"))
})

dat.highbins.lst <- lapply(infs.highbins, function(inf){
  fread(inf)
})

dat.debins.lst <- lapply(infs.debins, function(inf){
  fread(inf)
})

dat.debins.long <- lapply(jmarks, function(jmark){
  jdat <- dat.debins.lst[[jmark]]
  jdat$mark <- jmark
  return(jdat)
}) %>%
  bind_rows()

dat.debins.long$mark <- factor(dat.debins.long$mark, levels = jmarks)

outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps/DE_bins_gc_and_distance_to_gene.", Sys.Date(), ".pdf")

pdf(outpdf, useDingbats = FALSE)

ggplot(dat.debins.long, aes(y = gc, x = mark)) + 
  geom_boxplot(alpha = 0.25) + 
  theme_bw() + 
  xlab("") + ylab("gc") + 
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.debins.long, aes(y = log10(abs(distanceToTSS)), x = mark)) + 
  geom_boxplot(alpha = 0.25) + 
  theme_bw() + 
  xlab("") + ylab("dist to TSS") + 
  geom_hline(yintercept = log10(25000), linetype = "dotted") + 
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()




