# Jake Yeung
# Date of Creation: 2020-08-26
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged_for_MARA/2-check_mouse_BM_glmpca.R
# Check mouse BM glmpca

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(glmpca)

library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Load annots BM ----------------------------------------------------------

# jmark <- "H3K4me3"
# jmark <- "H3K27me3"
# jmark <- "H3K4me1"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

for (jmark in jmarks){
  
  
  # jsuffix <- "RemoveLargePeaks2"
  jsuffix <- "RemoveSmallPeaks"
  pdfout <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/MARA_output_BM/MARA_from_GLMPCA.", jmark, ".",  jsuffix, ".pdf")
  
  inf.glmpca <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.Faster.pois.", jsuffix, "KeepBestPlates2.LDAfromPeaks/PZ_", 
                       jmark, 
                       ".AllMerged.KeepBestPlates2.LDAfromPeaks.GLMPCA_var_correction.mergebinsize_100.binskeep_250.covar_cell.var.within.sum.norm.log2.CenteredAndScaled.penalty_1.5.winsorize_FALSE.2020-08-26.RData")
  assertthat::assert_that(file.exists(inf.glmpca))
  
  inf.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables/BM_AllMerged.", jmark, ".cell_cluster_table.txt")
  mdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged2_Peaks_1000/", jmark, "/mara_output/BM_", jmark, ".BM_AllMerged.VarCorrection.Poisson.GLMPCA.", jsuffix, "-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, "/BM_", jmark, ".BM_AllMerged.VarCorrection.Poisson.GLMPCA.", jsuffix)
  assertthat::assert_that(dir.exists(mdir))
  
  load(inf.glmpca, v=T)
  
  # Load GLMPCA output ------------------------------------------------------
  
  dat.umap <- DoUmapAndLouvain(glm.out$factors, jsettings = jsettings) %>%
    left_join(., dat.merge2, by = "cell")
  
  ggplot(dat.umap, aes(x = umap1.x, y = umap2.x, color = cell.var.within.sum.norm)) + 
    geom_point()  + 
    theme_bw() + 
    scale_color_viridis_c(direction = -1) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
  ggplot(dat.umap, aes(x = umap1.y, y = umap2.y, color = cell.var.within.sum.norm)) + 
    geom_point()  + 
    theme_bw() + 
    scale_color_viridis_c(direction = -1) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
  
  # Add cell annots ---------------------------------------------------------
  
  
  dat.annot <- fread(inf.annot)
  
  dat.umap.annot <- left_join(dat.umap, dat.annot, by = "cell")
  
  cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
  ggplot(dat.umap.annot, aes(x = umap1.x, y = umap2.x, color = cluster)) + 
    geom_point()  + 
    theme_bw() + 
    scale_color_manual(values = cbPalette) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
  
  # Add activity  -----------------------------------------------------------
  
  
  mara.out <- LoadMARA(mdir = mdir, make.cnames = FALSE)
  
  zscores.sub <- subset(mara.out$zscores, zscore > 0.7)
  motifs.keep <- zscores.sub$motif
  
  
  
  # Overlay activities onto UMAP  -------------------------------------------
  
  act.mat.clean <- t(subset(mara.out$act.mat, select = -motif))
  colnames(act.mat.clean) <- mara.out$act.mat$motif
  act.mat.clean.dat <- data.frame(cell = rownames(act.mat.clean), act.mat.clean, stringsAsFactors = FALSE) %>% 
    ungroup() %>%
    mutate(cell = gsub("\\.", "-", cell))
  
  
  jmat <- as.matrix(subset(act.mat.clean.dat, select = -cell))
  rownames(jmat) <- act.mat.clean.dat$cell
  
  jmat.sub <- jmat[, motifs.keep]
  
  # overlay 
  
  print(motifs.keep)
  
  
  # (jmotif <- motifs.keep[[2]])
  # (jmotif <- motifs.keep[[15]])
  
  # jmotif <- "Runx2"
  # jmotif <- "Hoxc6"
  # jmotif <- ""
  
  # (jmotif <- motifs.keep[[1]])
  
  
  cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
  
  
  pdf(file = pdfout, useDingbats = FALSE)
  
  m.var <- ggplot(dat.umap.annot, aes(x = umap1.x, y = umap2.x, color = cell.var.within.sum.norm)) + 
    geom_point()  + 
    theme_bw() + 
    scale_color_viridis_c(direction = -1) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m.var)
  
  m.clstr <- ggplot(dat.umap.annot %>% filter(cluster != ""), aes(x = umap1.x, y = umap2.x, color = cluster)) + 
    geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette)
  print(m.clstr)
  
  for (jmotif in motifs.keep){
    print(jmotif)
    jsub.dat <- data.frame(activity = jmat[, jmotif], cell = rownames(jmat.sub), stringsAsFactors = FALSE) %>%
      left_join(., dat.umap.annot)
    m1 <- ggplot(jsub.dat, aes(x = umap1.x, y = umap2.x, color = activity)) + 
      geom_point() + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle(jmotif) + 
      scale_color_viridis_c()
    print(m1)
  }
  dev.off()
  
  
  
  
}




