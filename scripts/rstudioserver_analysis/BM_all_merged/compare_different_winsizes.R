# Jake Yeung
# Date of Creation: 2020-01-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/compare_different_winsizes.R
# Compare sliding windows 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)
library(hash)
library(igraph)
library(umap)

library(scchicFuncs)

# Load data  --------------------------------------------------------------

jmark <- "H3K4me3"

jwins <- c("100000_20000", "50000_25000", "20000_10000", "50000_50000", "20000_20000")


# jwin <- "20000_20000"
# jwin <- "50000_50000"

outdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/pdfs_BM_OtherWinSize"

for (jwin in jwins){
  outf <- file.path(outdir, paste0("umaps_", jwin, ".pdf"))
  if (file.exists(outf)){
    next
  }
  if (jwin == "20000_20000"){
    inf <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-12.bsizestepsize_", jwin, "/lda_outputs.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj")
  } else if (jwin == "50000_50000"){
    inf <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-12.bsizestepsize_", jwin, ".NoSliding/lda_outputs.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj")
  } else if (jwin == "50000_25000"){
    inf <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2019-12-30.bsizestepsize_", jwin, "/lda_outputs.count_mat.H3K4me3.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me3.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj")
  } else if (jwin == "20000_10000"){
    inf <- "/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2019-12-30.bsizestepsize_20000_10000/lda_outputs.count_mat.H3K4me3.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me3.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj"
  } else if (jwin == "100000_20000"){
    inf <- "/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2019-12-16/lda_outputs.B6BM_AllMerged_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-16.K-50.binarize.FALSE/ldaOut.B6BM_AllMerged_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-16.K-50.Robj"
  } else {
    stop(paste("Unknown jwin:", jwin))
  }
  
  "/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-12.bsizestepsize_20000_20000/lda_outputs.count_mat.H3K4me3.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me3.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj"
  
  # inf.check <- "/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-12.bsizestepsize_20000_20000/lda_outputs.count_mat.H3K4me3.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me3.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj"
  
  assertthat::assert_that(file.exists(inf))
  
  load(inf, v=T)
  
  topics.mat <- posterior(out.lda)$topics
  terms.mat <- posterior(out.lda)$terms
  
  jsettings <- umap.defaults
  jsettings$n_neighbors <- 30
  jsettings$min_dist <- 0.1
  jsettings$random_state <- 123
  umap.out <- umap(topics.mat, config = jsettings)
  dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
  dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
  
  cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
  
  
  dat.impute.log <- log2(t(topics.mat %*% terms.mat))
  jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
  
  dat.var <- CalculateVarAll(dat.impute.log, jchromos)
  
  dat.merge <- left_join(dat.umap.long, dat.var)
  
  m.louv <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
    ggtitle(jwin)
  m.var <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) +  
    scale_color_viridis_c(direction = -1) + ggtitle(jwin)
  
  pdf(file = outf, useDingbats = FALSE)
    print(m.louv)
    print(m.var)
  dev.off()
}


