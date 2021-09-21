# Jake Yeung
# Date of Creation: 2021-02-22
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/dendrograms/7-show_variance_in_GLMPCA.R
# Plot factors for glmpca for all 4 marks

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(glmpca)
library(scchicFuncs)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"

indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage")

dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  fread(file.path(indir.meta, fname)) %>%
    rowwise()
}) 


# Load glmpcas ------------------------------------------------------------


# jmark <- "H3K9me3"
niter <- "100"
inf.glmpca.lst <- lapply(jmarks, function(jmark){
  if (jmark != "H3K27me3"){
    inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/BM_from_matadj/glmpca.", jmark, ".from_matadj.platename_jrep.szname_none.niter_", niter, ".RData"))
    assertthat::assert_that(file.exists(inf.glmpca))
  } else {
    inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/H3K27me3_rep2rep3reseq.peaks.varfilt/glmpca.", jmark, ".bincutoff_0.binskeep_0.platename_plate.szname_none.niter_500.reorder_rownames.dupfilt.suffix_peaks.RData"))
  }
  return(inf.glmpca)
})

glm.outs <- lapply(inf.glmpca.lst, function(inf){
  load(inf, v=T)
  return(glm.out)
})
for (jmark in jmarks){
  plot(glm.outs[[jmark]]$dev, main = jmark)
}



glmpca.factors <- lapply(inf.glmpca.lst, function(inf){
  load(inf, v=T)
  return(glm.out$factors)
})

dat.glmpca.factors <- lapply(jmarks, function(jmark){
  jtmp <- glmpca.factors[[jmark]]
  data.frame(cell = rownames(jtmp), pc1 = jtmp[, 1], pc2 = jtmp[, 2], pc3 = jtmp[, 3], pc4 = jtmp[, 4], pc30 = jtmp[, 30], stringsAsFactors = FALSE) %>%
    left_join(., dat.metas[[jmark]])
})


norm<-function(v){sqrt(sum(v^2))}

colNorms<-function(x){
  #compute the L2 norms of columns of a matrix
  apply(x,2,norm)
}


l2norms.lst <- lapply(jmarks, function(jmark){
  jfactors <-glmpca.factors[[jmark]]
  l2norms <- colNorms(jfactors)
  l2norms.frac <- l2norms / sum(l2norms)
  return(signif(l2norms.frac, digits = 2) * 100)
})


# fraction % total L2 norm



outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps/glmpca_dim_reductions.with_L2_norm.", Sys.Date(), ".pdf")

pdf(outpdf, useDingbats = FALSE)
for (jmark in jmarks){
  jsub <- dat.glmpca.factors[[jmark]]
  m <- ggplot(jsub, aes(x = pc1, y = pc2, color = clustercol)) + 
    geom_point(alpha = 0.5) + 
    xlab(paste0("PC1 ", l2norms.lst[[jmark]][[1]], " % total L2 norm")) + 
    ylab(paste0("PC2 ", l2norms.lst[[jmark]][[2]], " % total L2 norm")) + 
    ggtitle(jmark) + 
    scale_color_identity(guide = "legend",
                         labels = unique(jsub$cluster),
                         breaks = unique(jsub$clustercol)) +
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  m <- ggplot(jsub, aes(x = pc1, y = pc3, color = clustercol)) + 
    geom_point(alpha = 0.5) + 
    ggtitle(jmark) + 
    scale_color_identity(guide = "legend",
                         labels = unique(jsub$cluster),
                         breaks = unique(jsub$clustercol)) +
    xlab(paste0("PC1 ", l2norms.lst[[jmark]][[1]], " % total L2 norm")) + 
    ylab(paste0("PC3 ", l2norms.lst[[jmark]][[3]], " % total L2 norm")) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  m <- ggplot(jsub, aes(x = pc1, y = pc4, color = clustercol)) + 
    geom_point(alpha = 0.5) + 
    ggtitle(jmark) + 
    xlab(paste0("PC1 ", l2norms.lst[[jmark]][[1]], " % total L2 norm")) + 
    ylab(paste0("PC4 ", l2norms.lst[[jmark]][[4]], " % total L2 norm")) + 
    scale_color_identity(guide = "legend",
                         labels = unique(jsub$cluster),
                         breaks = unique(jsub$clustercol)) +
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  m <- ggplot(jsub, aes(x = pc1, y = pc30, color = clustercol)) + 
    geom_point(alpha = 0.5) + 
    ggtitle(jmark) + 
    scale_color_identity(guide = "legend",
                         labels = unique(jsub$cluster),
                         breaks = unique(jsub$clustercol)) +
    xlab(paste0("PC1 ", l2norms.lst[[jmark]][[1]], " % total L2 norm")) + 
    ylab(paste0("PC30 ", l2norms.lst[[jmark]][[30]], " % total L2 norm")) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
}
dev.off()



# Plot glmpca factors -----------------------------------------------------





