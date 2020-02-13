# Jake Yeung
# Date of Creation: 2020-02-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/8-downstream_LDA_varfilt_GLMPCA_correct_devmode_allmarks.KeepBestPlates2.final.R
# description

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

library(glmpca)

# library(devtools)
# dev_mode(on = TRUE)
# devtools::install_github("willtownes/glmpca")
library(glmpca)

# Constants ---------------------------------------------------------------

# jcovar.cname <- "ncuts.var.CenteredAndScaled"
# jcovar.cname <- "ncuts.var.CenteredAndScaled"
jcovar.cname <- "ncuts.var.log2.CenteredAndScaled"
# jcovar.cname <- "cell.var.within.sum.norm.CenteredAndScaled"

jconds <- c("Unenriched", "AllMerged"); names(jconds) <- jconds
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

jcondsmarks <- as.character(levels(interaction(jconds, jmarks, sep = "_")))
names(jcondsmarks) <- jcondsmarks
ncores <- length(jcondsmarks)


niter <- 250
jbins.keep <- 250
# calculating var raw
binsize <- 50000
mergesize <- 1000
bigbinsize <- 50000 * mergesize
jpenalty <- 2


jdate <- "2020-02-11"
ldadate <- "2020-02-11"
jsuffix <- "KeepBestPlates2"

# outdir <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.", jsuffix)
outdir <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.", jsuffix)
# outdir <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.", jsuffix)
dir.create(outdir)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

inmain <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.", jdate, ".var_filt.UnenrichedAndAllMerged.", jsuffix)
assertthat::assert_that(dir.exists(inmain))


# for (jcond in jconds){
  
mclapply(jcondsmarks, function(jcondmark){
  
  jcond <- strsplit(jcondmark, "_")[[1]][[1]]
  jmark <- strsplit(jcondmark, "_")[[1]][[2]]
  
  print(paste("Running for:", jcond, jmark))
  # lapply(jmarks, function(jmark){
  
  # Setup output paths ------------------------------------------------------

  outbase <- paste0("PZ_", jmark, ".", jcond, ".", jsuffix, ".GLMPCA_var_correction.mergebinsize_", mergesize, ".binskeep_", jbins.keep, ".covar_", jcovar.cname, ".penalty_", jpenalty, ".", jdate)
  outname <- paste0(outbase, ".RData")
  outname.pdf <- paste0(outbase, ".pdf")
  outf <- file.path(outdir, outname)
  outf.pdf <- file.path(outdir, outname.pdf)
  if (file.exists(outf)){
    print(paste("Outf exists, skipping...", outf))
    next
  }
  
  # Load data  --------------------------------------------------------------
  
  infname <- paste0("lda_outputs.BM_", jmark, "_varfilt_countmat.", ldadate, ".", jcond, ".K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.", ldadate, ".", jcond, ".K-30.Robj")
  inf <- file.path(inmain, infname)
  assertthat::assert_that(file.exists(inf))
  
  load(inf, v=T)
  count.mat <- as.matrix(count.mat)
  tm.result <- posterior(out.lda)
  colnames(tm.result$topics) <- paste("topic", colnames(tm.result$topics), sep = "_")
  rownames(tm.result$terms) <- paste("topic", rownames(tm.result$terms), sep = "_")
  topics.mat <- tm.result$topics
  

  # start -----------------
  pdf(outf.pdf, useDingbats = FALSE)
 
  print(jmark)
  print("Current time elapsed:")
  print(Sys.time() - jstart)
  
  
  
  # Plot it all -------------------------------------------------------------
  
  
  
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
  
  # PCA on topics, show variance
  pca.out <- prcomp(topics.mat, center = TRUE, scale. = TRUE)
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
  
  
  
  
  dat.var.raw <- CalculateVarRaw(count.mat, merge.size = mergesize, chromo.exclude.grep = "^chrX|^chrY", jpseudocount = 1, jscale = 10^6, calculate.ncuts = TRUE)
  # center and scale ncuts.var
  dat.var.raw$ncuts.var.log2 <- log2(dat.var.raw$ncuts.var)
  dat.var.raw$ncuts.var.CenteredAndScaled <- (dat.var.raw$ncuts.var - mean(dat.var.raw$ncuts.var)) / sd(dat.var.raw$ncuts.var)
  dat.var.raw$ncuts.var.log2.CenteredAndScaled <- (dat.var.raw$ncuts.var.log2 - mean(dat.var.raw$ncuts.var.log2)) / sd(dat.var.raw$ncuts.var.log2)
  
  dat.merge2 <- left_join(dat.var.merge, dat.var.raw)
  dat.merge2$cell.var.within.sum.norm.CenteredAndScaled <- (dat.merge2$cell.var.within.sum.norm - mean(dat.merge2$cell.var.within.sum.norm)) / sd(dat.merge2$cell.var.within.sum.norm)
  
  
  m1 <- ggplot(dat.merge2, aes(x = ncuts.var, y = cell.var.within.sum.norm)) + geom_point() + 
    scale_x_log10() + scale_y_log10()
  m2 <- ggplot(dat.merge2, aes(x = ncuts, y = ncuts.var)) + geom_point() + 
    scale_x_log10() + scale_y_log10()
  m3 <- ggplot(dat.merge2, aes_string(x = jcovar.cname, y = "cell.var.within.sum.norm")) + geom_point() + 
  scale_y_log10()
  
  print(m1)
  print(m2)
  print(m3)
  
  dev.off()
  
  # set up GLMPCA
  
  glm.inits <- InitGLMPCAfromLDA(count.mat, tm.result, dat.merge2, covar.cname = jcovar.cname, bins.keep = jbins.keep, do.log = FALSE, svd.on.Yinit = TRUE, use.orig.sz = TRUE)
  
  
  system.time(
    glm.out <- glmpca(Y = glm.inits$Y.filt, L = glm.inits$ntopics, fam = "mult", 
                      ctl = list(maxIters = niter, eps = 1e-4), 
                      penalty = jpenalty, verbose = TRUE, init = list(factors = glm.inits$U.init, loadings = glm.inits$V.init), 
                      X = glm.inits$X.mat, Z = NULL, sz = glm.inits$size.factor)
  )
  print(traceback())
  
  save(glm.out, glm.inits, dat.merge2, file = outf)
  
  
  print(Sys.time() - jstart)
  
}, mc.cores = ncores)
# })

# }

# dev_mode(on = FALSE)
