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
library(topicmodels)

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
  jexperi <- "AllMerged"
  inf.lda <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-11.", jexperi, ".K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-11.", jexperi, ".K-30.Robj")
  assertthat::assert_that(file.exists(inf.lda))
  load(inf.lda, v=T)
  tm.result <- posterior(out.lda)
  jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
  dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
  dat.var <- CalculateVarAll(dat.impute.log, jchromos)
  
  
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_TSS.first_transcript/", jmark, ".countTableTSS.mapq_40.TSS_", jdist, ".blfiltered.csv")
  assertthat::assert_that(file.exists(inf))
  
  mat <- ReadMatTSSFormat(inf, as.sparse = TRUE)
  
  # Plot fraction of signal at TSS ------------------------------------------
  
  dat.tss <- data.frame(cell = colnames(mat), tss.cuts = colSums(mat), stringsAsFactors = FALSE)
  
  # check fraction of TSS's with a count
  mat.bin <- BinarizeMatrix(mat)
  
  dat.tss.ngenes <- data.frame(cell = colnames(mat.bin), ngenes = colSums(mat.bin), ngenes.frac = colSums(mat.bin) / nrow(mat.bin), stringsAsFactors = FALSE)
  
  dat.size <- data.frame(cell = names(glm.inits$size.factor), total.cuts = glm.inits$size.factor, stringsAsFactors = FALSE)
  
  dat.merge <- left_join(dat.tss, dat.size)
  dat.merge <- left_join(dat.merge, dat.umap.glm.fillNAs)
  dat.merge <- left_join(dat.merge, dat.tss.ngenes)
  dat.merge <- left_join(dat.merge, dat.var)
  
  m <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = tss.cuts / total.cuts)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = 1) + ggtitle(paste(jmark, jdist, "around TSS"))
  m.rev <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = tss.cuts / total.cuts)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) + ggtitle(paste(jmark, jdist, "around TSS"))
  
  m.ngenes <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = log10(ngenes))) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = 1) + ggtitle(paste(jmark, jdist, "around TSS"))
  m.ngenes.rev <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = log10(ngenes))) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) + ggtitle(paste(jmark, jdist, "around TSS"))
  m.ngenes.frac <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = log10(ngenes.frac))) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = 1) + ggtitle(paste(jmark, jdist, "around TSS"))
  m.ngenes.frac.rev <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = log10(ngenes.frac))) + geom_point() + 
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
  
  dat.var.raw <- data.frame(cell = rownames(glm.inits$X.mat), ncuts.var = glm.inits$X.mat[, 1], stringsAsFactors = FALSE)
  
  dat.merge <- left_join(dat.merge, dat.var.raw)
  
  m.scatter <- ggplot(dat.merge, mapping = aes(x = tss.cuts / total.cuts, y = ncuts.var, color = cond)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m.scatter.impute <- ggplot(dat.merge, mapping = aes(x = tss.cuts / total.cuts, y = cell.var.within.sum.norm, color = cond)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # check var 
  head(dat.var)
  
  dat.var.merge <- left_join(dat.merge, dat.var)
  
  p1 <- ggplot(dat.var.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c()
  
  p2 <- ggplot(dat.var.merge, aes(x = umap1, y = umap2, color = ncuts.var)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c()
  
  p3 <- ggplot(dat.var.merge, aes(x = cell.var.within.sum.norm, y = tss.cuts)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  p4 <- ggplot(dat.var.merge, aes(x = cell.var.within.sum.norm, y = tss.cuts / total.cuts, color = cluster)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  p5 <- ggplot(dat.var.merge, aes(x = cell.var.within.sum.norm, y = tss.cuts / total.cuts, color = cluster)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  p6 <- ggplot(dat.var.merge, aes(x = cluster, y = tss.cuts / total.cuts, fill = cluster)) + geom_boxplot()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  p7 <- ggplot(dat.var.merge, aes(x = cluster, y = cell.var.within.sum.norm, fill = cluster)) + geom_boxplot()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  # plot outputs
  
  pdf(outpdf, useDingbats = FALSE)
    print(m)
    print(m.rev)
    print(m.plates)
    print(m.plates.rev)
    print(m.dens) 
    print(m.ngenes)
    print(m.ngenes.frac)
    print(m.ngenes.rev)
    print(m.ngenes.frac.rev)
    print(m.scatter)
    print(m.scatter.impute)
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    print(p5)
    print(p6)
    print(p7)
  dev.off()
  # save objs for easier loading llater
  if (!file.exists(outrdata)){
    save(dat.de, mat, dat.merge, file = outrdata)
  } else {
    print(paste("outrdata exists, continuing:", outrdata))
  }
}



