# Jake Yeung
# Date of Creation: 2021-03-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/dendrograms/9-compare_HSPCs_nonHSPCs_log2FC_fits.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

jout.lst <- lapply(jmarks, function(jmark.tmp){
  print(jmark.tmp)
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.HSPCs_vs_nonHSPCs/poisson_fit_bins.", jmark.tmp, ".2021-02-04.newannot2.witherrors.MoreBins.HSPCs_vs_nonHSPCs.total.RData")
  load(inf, v=T)
  jout <- SummarizeParamsPvalues(jfits.lst, jmark = jmark.tmp)
  return(jout)
})



# Load DE bins ------------------------------------------------------------


infs.de.lst <- lapply(jmarks, function(jmark){
  inf.tmp <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/DE_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.jmidbug_fixed.", jmark, ".2021-02-19.txt")
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

dat.de.bins.lst <- lapply(infs.de.lst, function(jinf){
  fread(jinf)
})

dat.params.filt.long <- lapply(jmarks, function(jmark){
  subset(jout.lst[[jmark]], bin %in% dat.de.bins.lst[[jmark]]$CoordOriginal & abs(estimate) < 5)
})  %>%
  bind_rows()

dat.params.filt.long$mark <- factor(dat.params.filt.long$mark, levels = jmarks)

# Show log2FCs ------------------------------------------------------------

outpdf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps/HSPC_vs_nonHSPCs_log2FC_dynamic_bins.pdf"
pdf(outpdf, useDingbats = FALSE)
ggplot(dat.params.filt.long %>% bind_rows(), aes(x = estimate / log(2), fill = mark)) + 
  geom_density(alpha = 0.25, fill = "grey25") + 
  facet_wrap(~mark, nrow = 1) + 
  theme_bw() +  
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  xlab("log2FC non-HSPCs vs HSPCs")

dev.off()


