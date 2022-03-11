# Jake Yeung
# Date of Creation: 2022-03-05
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/11-pseudotime_downstream_k27_k9_effect_sizes.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jctypes <- c("Eryths", "Bcells", "Granulocytes"); names(jctypes) <- jctypes

dat.hits.merge.long <- lapply(jctypes, function(jctype){
  print(jctype)
  
  inf.hits.k9 <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/ptimefits_plots_top_hits/fits_with_annotations.", jctype, ".k9me3.2022-02-06.neg_slope.txt")
  inf.hits.k27 <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/ptimefits_plots_top_hits/fits_with_annotations.", jctype, ".k27me3.2022-02-05.neg_slope.txt")
  
  dat.hits.k9 <- fread(inf.hits.k9) %>%
    mutate(mark = "k9me3")
  dat.hits.k27 <- fread(inf.hits.k27) %>%
    mutate(mark = "k27me3")
  dat.hits.merge <- rbind(dat.hits.k9, dat.hits.k27) %>%
    mutate(ctype = jctype)
  return(dat.hits.merge)
}) %>%
  bind_rows()

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

outpdf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/ptimefits_hits_k27_k9_together/summary_plot_effect_sizes.", Sys.Date(), ".pdf")
pdf(outpdf, useDingbats = FALSE)
m.dens <- ggplot(dat.hits.merge.long, aes(x = estimate, fill = mark)) + 
  geom_density(alpha = 0.5) + 
  ggtitle("Effect size estimates") + 
  scale_fill_manual(values = cbPalette) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  facet_wrap(~ctype) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") 
print(m.dens)


m.boxplot <- ggplot(dat.hits.merge.long, aes(x = ctype, y = estimate, fill = mark)) + 
  geom_violin() + 
  ggtitle("Effect size estimates") + 
  scale_fill_manual(values = cbPalette) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") 
print(m.boxplot)
dev.off()
