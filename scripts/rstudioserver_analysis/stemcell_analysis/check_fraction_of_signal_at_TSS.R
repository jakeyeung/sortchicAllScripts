# Jake Yeung
# Date of Creation: 2020-03-22
# File: ~/projects/scchic/scripts/rstudioserver_analysis/stemcell_analysis/check_fraction_of_signal_at_TSS.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jdist <- "10000"

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")



for (jmark in jmarks){
  
  outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/stemcell_analysis/proxfrac.TSS_", jdist, ".", jmark, ".", Sys.Date(), ".pdf")
  outrdata <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/stemcell_analysis/proxfrac.TSS_", jdist, ".", jmark, ".", Sys.Date(), ".RData")
  
  # Load DE genes -----------------------------------------------------------
  
  inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Giladi_et_al_2018/diff_exprs_Giladi_seurat.celltypes_filt.rds"
  dat.de <- readRDS(inf.de)
  
  
  # Load UMAP annots --------------------------------------------------------
  
  inf.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2/PZ_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_var_correction.mergebinsize_1000.binskeep_1000.covar_ncuts.var.log2.CenteredAndScaled.penalty_1.winsorize_TRUE.2020-02-11.RData")
  load(inf.annot, v=T)
  
  inf.annot2 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
  load(inf.annot2, v=T)
  
  
  
  # Load LDA  ---------------------------------------------------------------
  
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_TSS.first_transcript/", jmark, ".countTableTSS.mapq_40.TSS_", jdist, ".blfiltered.csv")
  assertthat::assert_that(file.exists(inf))
  
  mat <- ReadMatTSSFormat(inf, as.sparse = TRUE)
  
  # Plot fraction of signal at TSS ------------------------------------------
  
  dat.tss <- data.frame(cell = colnames(mat), tss.cuts = colSums(mat), stringsAsFactors = FALSE)
  dat.size <- data.frame(cell = names(glm.inits$size.factor), total.cuts = glm.inits$size.factor, stringsAsFactors = FALSE)
  dat.merge <- left_join(dat.tss, dat.size)
  dat.merge <- left_join(dat.merge, dat.umap.glm.fillNAs)
  
  m <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = tss.cuts / total.cuts)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = 1) + ggtitle(paste(jmark, jdist, "around TSS"))
  m.rev <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = tss.cuts / total.cuts)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) + ggtitle(paste(jmark, jdist, "around TSS"))
  
  m.plates <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = tss.cuts / total.cuts)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) + ggtitle(paste(jmark, jdist, "around TSS")) + facet_wrap(~cond)
  m.plates.rev <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = tss.cuts / total.cuts)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = 1) + ggtitle(paste(jmark, jdist, "around TSS")) + facet_wrap(~cond)
  

  # Density of TSS signal across different conditions ---------------------
  m.dens <- ggplot(dat.merge, aes(x = tss.cuts / total.cuts, group = cond, fill = cond)) + geom_density(alpha = 0.33) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = 1) + ggtitle(paste(jmark, jdist, "around TSS")) + 
    xlab(paste0("Fraction of TSS cuts in genome (", jdist, " bp window)"))  + 
    scale_fill_manual(values = cbPalette)

  
  # plot outputs
  
  pdf(outpdf, useDingbats = FALSE)
    print(m)
    print(m.rev)
    print(m.plates)
    print(m.plates.rev)
    print(m.dens) 
  dev.off()
  # save objs for easier loading llater
  if (!file.exists(outrdata)){
    save(dat.de, mat, dat.merge, file = outrdata)
  } else {
    print(paste("outrdata exists, continuing:", outrdata))
  }
}



