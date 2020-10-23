# Jake Yeung
# Date of Creation: 2020-08-26
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged_for_MARA/1-correct_LDA_from_peaks_by_variance.R
# 


rm(list=ls())

jstart <- Sys.time() 

library(parallel)

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
library(DescTools)
library(glmpca)

# Constants ---------------------------------------------------------------

# jcovar.cname <- "ncuts.var.log2.CenteredAndScaled"
jcovar.cname <- "cell.var.within.sum.norm.CenteredAndScaled"
winsorize <- FALSE

jcond <- "AllMerged"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmarks <- c("H3K4me1"); names(jmarks) <- jmarks
jvarcutoffs <- c(0.75, 2, 1, 0.5)
names(jvarcutoffs) <- jmarks
ncores <- length(jmarks)

niter <- 250
jbins.keep <- 250
# calculating var raw
binsize <- 50000
mergesize <- 1000
bigbinsize <- 50000 * mergesize
jpenalty <- 1.5


jdate <- "2020-04-29"
jsuffix <- "ZF_AllMerged.imputevarfilt.lessstringent.LDAfromPeaks"

hubprefix <- "/home/jyeung/hub_oudenaarden"

outdir <- paste0(hubprefix, "/jyeung/data/zebrafish_scchic/GLMPCA_outputs.", jsuffix)
dir.create(outdir)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

inmain.bins <- paste0(hubprefix, "/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZF_AllMerged.winsize_50000.imputevarfilt.lessstringent")
inmain.peaks <- paste0(hubprefix, "/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZFWKM_peaks")
assertthat::assert_that(dir.exists(inmain.bins))
assertthat::assert_that(dir.exists(inmain.peaks))

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
(jchromos <- paste("chr", c(seq(25)), sep = ""))


# 


# jmark <- jmarks[[1]]
jout <- mclapply(jmarks, function(jmark){
  
  print(paste("Running for:", jcond, jmark))
  
  # Setup output paths ------------------------------------------------------
  
  outbase <- paste0("ZF_", jmark, ".", jcond, ".", jsuffix, ".GLMPCA_var_correction.mergebinsize_", mergesize, ".binskeep_", jbins.keep, ".covar_", jcovar.cname, ".penalty_", jpenalty, ".winsorize_", winsorize, ".", jdate)
  outname <- paste0(outbase, ".RData")
  outname.pdf <- paste0(outbase, ".pdf")
  outf <- file.path(outdir, outname)
  outf.pdf <- file.path(outdir, outname.pdf)
  
  if (file.exists(outf)){
    return(NULL)
  }
  
  # Load data  --------------------------------------------------------------
  
  jvarcutoff <- jvarcutoffs[[jmark]]
  infname.bins <- paste0("lda_outputs.counts_table_var_filt.", jmark, ".imputevar_", jvarcutoff, ".K-30.binarize.FALSE/ldaOut.counts_table_var_filt.", jmark, ".imputevar_", jvarcutoff, ".K-30.Robj")
  infname.peaks <- paste0("lda_outputs.", jmark, ".imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-30.binarize.FALSE/ldaOut.", jmark, ".imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-30.Robj")
  inf.bins <- file.path(inmain.bins, infname.bins)
  assertthat::assert_that(file.exists(inf.bins))
  inf.peaks <- file.path(inmain.peaks, infname.peaks)
  assertthat::assert_that(file.exists(inf.peaks))
  
  # load raw counts from peaks 
  load(inf.peaks, v=T)
  count.mat.peaks <- as.matrix(count.mat)
  tm.result.peaks <- AddTopicToTmResult(posterior(out.lda), jsep = "_")
  topics.mat.peaks <- tm.result.peaks$topics
  
  load(inf.bins, v=T)
  count.mat.bins <- as.matrix(count.mat)
  tm.result.bins <- posterior(out.lda)
  tm.result.bins <- AddTopicToTmResult(posterior(out.lda), jsep = "_")
  # colnames(tm.result.bins$topics) <- paste("topic", colnames(tm.result.bins$topics), sep = "_")
  # rownames(tm.result.bins$terms) <- paste("topic", rownames(tm.result.bins$terms), sep = "_")
  # topics.mat.bins <- tm.result.bins$topics
  
  
  
  # start -----------------
  pdf(outf.pdf, useDingbats = FALSE)
  
  print(jmark)
  print("Current time elapsed:")
  print(Sys.time() - jstart)
  
  
  
  # Plot it all -------------------------------------------------------------
  
  
  umap.out <- umap(topics.mat.peaks, config = jsettings)
  dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
  dat.umap.long <- DoLouvain(topics.mat.peaks, jsettings, dat.umap.long)
  ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)
  
  
  # Plot variance -----------------------------------------------------------
  
  
  dat.impute.log <- log2(t(tm.result.bins$topics %*% tm.result.bins$terms))
  
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
  
  # PCA on topics, show variance
  pca.out <- prcomp(topics.mat.peaks, center = TRUE, scale. = TRUE)
  dat.pca <- data.frame(cell = rownames(pca.out$x), PC1 = pca.out$x[, 1], PC2 = pca.out$x[, 2], stringsAsFactors = FALSE) %>%
    left_join(dat.var.merge)
  
  m.pca.var <- ggplot(dat.pca, aes(x = PC1, y = PC2, color = cell.var.within.sum.norm)) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_point() + scale_color_viridis_c(direction = -1)
  
  m.pca.var.plates <- ggplot(dat.pca, aes(x = PC1, y = PC2, color = cell.var.within.sum.norm)) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_point() + scale_color_viridis_c(direction = -1) + facet_wrap(~plate)
  
  print(m.plates)
  print(m.var)
  print(m.var.plates)
  print(m.pca.var)
  print(m.pca.var.plates)
  
  
  
  # Calculate raw varaince and compare with imputed variance  ---------------
  
  dat.var.raw <- CalculateVarRaw(count.mat.bins, merge.size = mergesize, chromo.exclude.grep = "^chrX|^chrY", jpseudocount = 1, jscale = 10^6, calculate.ncuts = TRUE)
  # center and scale ncuts.var
  dat.var.raw$ncuts.var.log2 <- log2(dat.var.raw$ncuts.var)
  dat.var.raw$ncuts.var.CenteredAndScaled <- (dat.var.raw$ncuts.var - mean(dat.var.raw$ncuts.var)) / sd(dat.var.raw$ncuts.var)
  dat.var.raw$ncuts.var.log2.CenteredAndScaled <- (dat.var.raw$ncuts.var.log2 - mean(dat.var.raw$ncuts.var.log2)) / sd(dat.var.raw$ncuts.var.log2)
  
  dat.merge2 <- left_join(dat.var.merge, dat.var.raw)
  dat.merge2$cell.var.within.sum.norm.CenteredAndScaled <- (dat.merge2$cell.var.within.sum.norm - mean(dat.merge2$cell.var.within.sum.norm)) / sd(dat.merge2$cell.var.within.sum.norm)
  
  # winsorize?
  if (winsorize){
    dat.merge2[[jcovar.cname]] <- DescTools::Winsorize(dat.merge2[[jcovar.cname]], probs = c(0.01, 0.99))
  }
  
  m1 <- ggplot(dat.merge2, aes(x = ncuts.var, y = cell.var.within.sum.norm)) + geom_point() + 
    scale_x_log10() + scale_y_log10() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  m2 <- ggplot(dat.merge2, aes(x = ncuts, y = ncuts.var)) + geom_point() + 
    scale_x_log10() + scale_y_log10() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  m3 <- ggplot(dat.merge2, aes_string(x = jcovar.cname, y = "cell.var.within.sum.norm")) + geom_point() + scale_y_log10() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  print(m1)
  print(m2)
  print(m3)
  
  dev.off()
  
  if (file.exists(outf)){
    print(paste("Outf exists, skipping...", outf))
    return(NULL)
  }
  
  # set up GLMPCA
  print("Running GLM inits")
  
  glm.inits <- InitGLMPCAfromLDA(count.mat.peaks, tm.result.peaks, dat.merge2, covar.cname = jcovar.cname, bins.keep = jbins.keep, do.log = FALSE, svd.on.Yinit = TRUE, use.orig.sz = TRUE)
  
  
  system.time(
    glm.out <- glmpca(Y = glm.inits$Y.filt, L = glm.inits$ntopics, fam = "mult", 
                      ctl = list(maxIters = niter, eps = 1e-4), 
                      penalty = jpenalty, verbose = TRUE, init = list(factors = glm.inits$U.init, loadings = glm.inits$V.init), 
                      X = glm.inits$X.mat, Z = NULL, sz = glm.inits$size.factor)
  )
  print(traceback())
  
  save(glm.out, glm.inits, dat.merge2, file = outf)
  
  
  print(Sys.time() - jstart)
  return(glm.out)
  
}, mc.cores = ncores)
# })

