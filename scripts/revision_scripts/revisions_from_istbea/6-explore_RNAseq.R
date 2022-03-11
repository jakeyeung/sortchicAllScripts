# Jake Yeung
# Date of Creation: 2022-02-28
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/6-explore_RNAseq.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


inf <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/analysis_scrnaseq_atacseq/Baccin_scRNAseq_bonemarrow_no_niche.2021-09-05.RData"
load(inf, v=T)

# load metadata
indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/analysis_scrnaseq_atacseq"
inf.meta <- file.path(indir.meta, paste0("metadata_CCA_celltype_colored.2022-02-23.txt"))
dat.meta <- fread(inf.meta)



# Plot some genes ---------------------------------------------------------

jmat <- Blood@assays$RNA@counts
# jmat <- Blood@assays$RNA@scale.data

# jgene <- "Yy"
# grep(jgene, rownames(jmat), value = TRUE)

plot(density(colSums(jmat)))

jmat.norm <- sweep(jmat, MARGIN = 2, STATS = colSums(jmat), FUN = "/") * 1000

# plot gene on UMAP

# jgene <- "Elane"
jgene <- "Zbtb"
jgene <- "Zfp14"


jgene <- "Yy1"
jgene <- "Zbtb16"

jgenes <- c("Yy1", "Zbtb16"); names(jgenes) <- jgenes

outpdf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/scrnaseq_marker_genes/marker_genes.", Sys.Date(), ".pdf")
pdf(outpdf, useDingbats = FALSE)

m1 <- ggplot(dat.meta, aes(x = umap1, y = umap2, color = clustercol)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_identity() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m1)

m1 <- ggplot(dat.meta, aes(x = umap1, y = umap2, color = cluster.rename)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~cluster.rename) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m1)

for (jgene in jgenes){
  (jgene.i <- grep(jgene, rownames(jmat), value = FALSE))
  rownames(jmat)[jgene.i]
  
  plot(density(jmat[jgene.i, ]))
  
  if (length(jgene.i) > 1){
    dat.exprs <- data.frame(cell = colnames(jmat.norm), normcounts = colMeans(jmat.norm[jgene.i, ]), stringsAsFactors = FALSE) %>%
      left_join(., dat.meta)
  } else {
    dat.exprs <- data.frame(cell = colnames(jmat.norm), normcounts = jmat.norm[jgene.i, ], stringsAsFactors = FALSE) %>%
      left_join(., dat.meta)
  }
  
  m1 <- ggplot(dat.exprs, aes(x = umap1, y = umap2, color = clustercol)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_identity() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m1)
  
  m <- ggplot(dat.exprs, aes(x = umap1, y = umap2, color = log(normcounts + 1))) + 
    geom_point() + 
    ggtitle(jgene)  + 
    theme_bw() + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}
dev.off()

