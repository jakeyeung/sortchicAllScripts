# Jake Yeung
# Date of Creation: 2020-03-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/stemcell_analysis/1-check_peak_signal.R
# Pseudobulk calling peaks. What is the genome distribution? 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(topicmodels)


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/peaks_analysis"
assertthat::assert_that(dir.exists(outdir))


for (jmark in jmarks){
  print(jmark)
  fname <- paste0("peaks_counts_vs_total.", jmark, ".pdf")
  outf <- file.path(outdir, fname)
  
    
  # Load UMAP annots  -------------------------------------------------------
  
  inf.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
  assertthat::assert_that(file.exists(inf.annot))
  load(inf.annot, v=T)
  
  inf.annot.glms <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2/PZ_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_var_correction.mergebinsize_1000.binskeep_1000.covar_ncuts.var.log2.CenteredAndScaled.penalty_1.winsorize_TRUE.2020-02-11.RData")
  load(inf.annot.glms, v=T)
  
  
  # Load data  --------------------------------------------------------------
  
  inf.bed <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.hiddenDomains_output/hd_merged.", jmark, ".minlength_1000/merged.", jmark, ".minlength_1000.cutoff_analysis.merged.withchr.annotated.bed")
  assertthat::assert_that(file.exists(inf.bed))
  
  dat.annot <- fread(inf.bed)
  colnames(dat.annot) <- c("chromo", "start", "end", "gene", "dist")
  
  inf.mat <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.hiddenDomains_output/count_tables_from_hiddenDomains_marks_merged.rds_format/", jmark, "-BM_AllMerged.merged_by_clusters_no_NAs.rds")
  peak.count.mat <- readRDS(inf.mat)
  
  dat.cellcounts <- data.frame(cell = names(glm.inits$size.factor), ncuts.total = glm.inits$size.factor, stringsAsFactors = FALSE)
  dat.peakcounts <- data.frame(cell = colnames(peak.count.mat), ncuts.peaks = colSums(peak.count.mat), stringsAsFactors = FALSE)
  
  # Plot counts in peaks vs total  ------------------------------------------
  
  dat.mergecounts <- left_join(dat.cellcounts, dat.peakcounts) %>%
    left_join(., dat.umap.glm.fillNAs)
  
  m1 <- ggplot(dat.mergecounts, aes(x = umap1, y = umap2, color = ncuts.peaks / ncuts.total)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = 1) + facet_wrap(~cond) + ggtitle(jmark)
  
  m2 <- ggplot(dat.mergecounts, aes(x = umap1, y = umap2, color = ncuts.peaks / ncuts.total)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = 1) + ggtitle(jmark)
  
  pdf(outf, useDingbats = FALSE)
    print(m1)
    print(m2)
  dev.off()
  
}
