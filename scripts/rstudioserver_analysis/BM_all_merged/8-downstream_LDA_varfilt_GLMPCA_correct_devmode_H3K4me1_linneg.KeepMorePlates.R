# Jake Yeung
# Date of Creation: 2020-02-02
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/8-downstream_LDA_varfilt_GLMPCA_correct_devmode_H3K4me1_linneg.R
# 

rm(list=ls())
jstart <- Sys.time()


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)

library(scchicFuncs)
library(hash)
library(igraph)
library(umap)

library(glmpca)

library(devtools)
dev_mode(on = TRUE)
devtools::install_github("willtownes/glmpca")
library(glmpca)

# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmark <- "H3K4me1"
niter <- 500
topn  <- 150
jbins.keep <- 1000

# jdate <- "2020-02-02"
KeepMorePlates <- TRUE


# for (jmark in jmarks){
  
  # Load annotations --------------------------------------------------------
  
  
  # inf.annot <- "/home/jyeung/hpc/scchic/public_data/Giladi_et_al_2018/giladi_pseudobulk_datsum_and_ncells.DESeq_and_QuantNorm.RData"
  inf.annot <- "/home/jyeung/hpc/scChiC/public_data/Giladi_et_al_2018/giladi_pseudobulk_datsum_and_ncells.DESeq_and_QuantNorm.RData"
  load(inf.annot, v=T)
  
  dat.sum.long <- data.table::melt(dat.sum.norm.quantnorm)
  colnames(dat.sum.long) <- c("gene", "celltype", "exprs")
  
  dat.sum.long <- dat.sum.long %>%
    group_by(gene) %>%
    mutate(zscore = scale(exprs, center = TRUE, scale = TRUE)) %>%
    filter(!is.nan(zscore))
  
  
  
  print(jmark)
  print("Current time elapsed:")
  print(Sys.time() - jstart)
  
  # Load data  --------------------------------------------------------------
  
  if (KeepMorePlates){
    inf <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-31.var_filt_keepPlates/lda_outputs.BM_", jmark, ".varcutoff_0.3.keepplates.WithLinneg.K-30.binarize.FALSE/ldaOut.BM_", jmark, ".varcutoff_0.3.keepplates.WithLinneg.K-30.Robj")
  } else {
    inf <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-31.var_filt_keepPlates/lda_outputs.BM_", jmark, ".varcutoff_0.3.platesRemoved.SmoothBinSize_1000.AllMerged.WithLinneg.K-30.binarize.FALSE/ldaOut.BM_", jmark, ".varcutoff_0.3.platesRemoved.SmoothBinSize_1000.AllMerged.WithLinneg.K-30.Robj")
  }
  load(inf, v=T)
  
  bname <- ClipLast(basename(inf), jsep = ".Robj")
  outbase <- paste0(bname, ".", topn, ".", Sys.Date(), ".binskeep_", jbins.keep, ".devmode")
  outname <- paste0(outbase, ".RData")
  outname.pdf <- paste0(outbase, ".pdf")
  # outdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/GLMPCA_outputs"
  outdir <- "/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs"
  outf <- file.path(outdir, outname)
  outf.pdf <- file.path(outdir, outname.pdf)
  pdf(outf.pdf, useDingbats = FALSE)
  
  tm.result <- posterior(out.lda)
  colnames(tm.result$topics) <- paste("topic", colnames(tm.result$topics), sep = "_")
  rownames(tm.result$terms) <- paste("topic", rownames(tm.result$terms), sep = "_")
  
  topics.mat <- tm.result$topics
  
  # Plot it all -------------------------------------------------------------
  
  jsettings <- umap.defaults
  jsettings$n_neighbors <- 30
  jsettings$min_dist <- 0.1
  jsettings$random_state <- 123
  
  
  umap.out <- umap(topics.mat, config = jsettings)
  dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
  dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
  cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
  ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)
  
  
  # Plot variance -----------------------------------------------------------
  
  (jchromos <- paste("chr", c(seq(19)), sep = ""))
  dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
  
  dat.var <- CalculateVarAll(dat.impute.log, jchromos) %>%
    rowwise() %>%
    mutate(plate = ClipLast(as.character(cell), jsep = "_"))
  
  dat.var.merge <- left_join(dat.umap.long, dat.var)
  
  m.plates <- ggplot(dat.var.merge, aes(x = umap1, y = umap2, color = plate)) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_point() + scale_color_manual(values = cbPalette)
  
  m.var <- ggplot(dat.var.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_point() + scale_color_viridis_c(direction = -1)
  
  m.var.plates <- ggplot(dat.var.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_point() + scale_color_viridis_c(direction = -1) + facet_wrap(~plate)
  
  print(m.plates)
  print(m.var)
  print(m.var.plates)
  
  
  
  # Calculate raw varaince and compare with imputed variance  ---------------
  
  binsize <- 50000
  mergesize <- 1000
  bigbinsize <- 50000 * mergesize
  
  
  dat.var.raw <- CalculateVarRaw(count.mat, merge.size = mergesize, chromo.exclude.grep = "^chrX|^chrY", jpseudocount = 1, jscale = 10^6, calculate.ncuts = TRUE)
  
  dat.merge2 <- left_join(dat.var.merge, dat.var.raw)
  
  ggplot(dat.merge2, aes(x = ncuts.var, y = cell.var.within.sum.norm)) + geom_point() + 
    scale_x_log10() + scale_y_log10()
  ggplot(dat.merge2, aes(x = ncuts, y = ncuts.var)) + geom_point() + 
    scale_x_log10() + scale_y_log10()
  
  # set up GLMPCA
  if (!file.exists(outf)){
    print("Correcting variance with GLMPCA")
    glm.inits <- InitGLMPCAfromLDA(as.matrix(count.mat), tm.result, dat.merge2, covar.cname = "ncuts.var", bins.keep = jbins.keep, do.log = FALSE, svd.on.Yinit = TRUE, use.orig.sz = TRUE)
    system.time(
      glm.out <- glmpca(Y = glm.inits$Y.filt, L = glm.inits$ntopics, fam = "mult", 
                        ctl = list(maxIters = niter, eps = 1e-4), 
                        penalty = 1, verbose = TRUE, init = list(factors = glm.inits$U.init, loadings = glm.inits$V.init), 
                        X = glm.inits$X.mat, Z = NULL, sz = glm.inits$size.factor)
    )
    print(traceback())
    save(glm.out, glm.inits, dat.merge2, file = outf)
  } else {
    print("Skipping GLM because outf already detected (useful to just make plots again)")
  }
  
  
  print(Sys.time() - jstart)
  
  dev_mode(on = FALSE)
  
  dev.off()
# }


