# Jake Yeung
# Date of Creation: 2022-04-17
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/21-show_sorted_ctypes_on_UMAP.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)

# jmarksnew <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarksnew) <- jmarksnew
jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks

# jmarknew <- "k4me1"
# jmarksold <- "H3K4me1"; names(jmarksold) <- jmarksold
# jmarkold <- "H3K4me1"

# Load meta ----------------------------------------------------------------

dat.meta.lst <- lapply(jmarks, function(jmark){
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths2/metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.", jmark, ".txt")
  dat.meta <- fread(inf.meta)
})

dat.meta.colors <- subset(dat.meta.lst$k4me1, select = c(ctype.from.LL, colcode))
ctype2col <- hash::hash(dat.meta.colors$ctype.from.LL, dat.meta.colors$colcode)


# Show some ctypes --------------------------------------------------------


# ctypes.show <- c("Granulocytes", "DCs", "Monocytes", "pDCs")

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/ctypes_on_umap"

for (jmark in jmarks){
  ctypes.show <- sort(unique(subset(dat.meta.lst[[jmark]], batch == "New")$ctype.from.LL))
  names(ctypes.show) <- ctypes.show
  
  dat.meta.sub <- dat.meta.lst[[jmark]]
  
  m <- ggplot(dat.meta.sub, aes(x = umap1, y = umap2, color = colcode)) +
    geom_point(alpha = 0.5) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  outpdf <- file.path(outdir, paste0("ctypes_on_umap.", jmark, ".pdf"))
  
  pdf(outpdf, useDingbats = FALSE)
  
  print(m)
  for(jctype in ctypes.show){
    
    jtitle <- paste(jmark, jctype)
    
    dat.meta.sub.label <- dat.meta.sub %>%
      rowwise() %>% 
      mutate(highlight.ctype = ctype == jctype & batch == "New",
             col.highlight = ifelse(highlight.ctype, colcode, "gray75"), 
             col.byctype = ifelse(highlight.ctype, ctype2col[[ctype]], "gray75"),
             jalpha = ifelse(highlight.ctype, 1, 0.5)) %>%
      arrange(highlight.ctype)
    
    
    if (all(!dat.meta.sub.label$highlight.ctype)){
      print(paste("Nothing labeled for", jctype, "skipping"))
      next
    }
    
    m2 <- ggplot(dat.meta.sub.label, aes(x = umap1, y = umap2, color = col.highlight, alpha = jalpha)) +
      geom_point() + 
      theme_bw() + 
      ggtitle(jtitle) + 
      scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                            guide = "legend") + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m2)
    
    m2.batch <- ggplot(dat.meta.sub.label, aes(x = umap1, y = umap2, color = col.highlight, alpha = jalpha)) +
      geom_point() + 
      theme_bw() + 
      facet_wrap(~batch) + 
      ggtitle(jtitle) + 
      scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                            guide = "legend") + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m2.batch)
    
    m3 <- ggplot(dat.meta.sub.label, aes(x = umap1, y = umap2, color = col.byctype, alpha = jalpha)) +
      geom_point() + 
      theme_bw() + 
      ggtitle(jtitle) + 
      scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                            guide = "legend") + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m3)
  }
  dev.off()
}





# Save outputs ------------------------------------------------------------


