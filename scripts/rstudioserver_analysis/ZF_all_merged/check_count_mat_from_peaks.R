# Jake Yeung
# Date of Creation: 2020-08-24
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/check_count_mat_from_peaks.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(irlba)

library(hash)
library(igraph)
library(umap)

# Load data  --------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
hubprefix <- "/home/jyeung/hub_oudenaarden"

inmain <- file.path(hubprefix, "jyeung/data/zebrafish_scchic/count_tables_all/count_tables.HiddenDomains.imputevarfilt.lessstringent.mapq_40.NewCountFilters")
assertthat::assert_that(dir.exists(inmain))
outmain <- file.path(inmain, "rds_mat_for_LDA")
dir.create(outmain)


for (jmark in jmarks){
  print(jmark)
  
  
  inf <- file.path(inmain, paste0(jmark, ".imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.csv"))
  outf <- file.path(outmain, paste0(jmark, ".imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.rds"))
  outpdf <- file.path(outmain, paste0(jmark, ".imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.pdf"))
  
  if (file.exists(outf)){
    print(paste(outf, "exists, skipping"))
    next
  }
  
  inf.annot <- file.path(hubprefix, paste0("jyeung/data/zebrafish_scchic/celltyping/from_louvain/cell_to_cluster_table.", jmark, ".2020-04-13.txt"))
  assertthat::assert_that(file.exists(inf.annot))
  dat.annot <- fread(inf.annot)
  
  assertthat::assert_that(file.exists(inf))
  assertthat::assert_that(!file.exists(outf))
  
  
  mat <- ReadMatTSSFormat(inf, add.coord = TRUE, sort.rnames = TRUE)
  
  lsi.out <- RunLSI(as.matrix(mat))
  
  jsettings <- umap.defaults
  jsettings$n_neighbors <- 30
  jsettings$min_dist <- 0.1
  jsettings$random_state <- 123
  
  umap.out <- DoUmapAndLouvain(lsi.out$u, jsettings)
  
  umap.out$experi <- sapply(umap.out$cell, function(x) ClipLast(x, jsep = "_"))
  
  m1 <- ggplot(umap.out, aes(x = umap1, y = umap2, color = experi)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette)
  
  
  # Get annots --------------------------------------------------------------
  
  
  umap.out.annot <- left_join(umap.out, subset(dat.annot, select = c(cell, cluster, plate)))
  
  
  m2 <- ggplot(umap.out.annot, aes(x = umap1, y = umap2, color = cluster)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette)
  
  pdf(outpdf, useDingbats = FALSE)
  print(m1)
  print(m2)
  dev.off()
  
  # Write tables for LDA  ---------------------------------------------------
  
  saveRDS(mat, file = outf)
  
}



