# Jake Yeung
# Date of Creation: 2020-02-07
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/9-downstream_GLMPCA_BM_KeepBestPlates.final.R
# 


jstart <- Sys.time()

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(glmpca)
library(scchicFuncs)
library(topicmodels)

library(hash)
library(igraph)
library(umap)

jconds <- c("Unenriched", "AllMerged"); names(jconds) <- jconds
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

jcondsmarks <- as.character(levels(interaction(jconds, jmarks, sep = "_")))
names(jcondsmarks) <- jcondsmarks

niter <- 1000
topn <- 150
jbins.keep <- 1000
# jbins.keep <- 150
# calculating var raw
binsize <- 50000
mergesize <- 1000
bigbinsize <- 50000 * mergesize

jdate <- "2020-02-06"
ldadate <- "2020-02-06"
jsuffix <- "KeepBestPlates"
  
outdir <- "/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.KeepBestPlates"  # load GLM outputs
pdfdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates.downstream"  # write to pdf dir
dir.create(pdfdir)

inmain <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.", ldadate, ".var_filt.UnenrichedAndAllMerged.", jsuffix)


# ncores <- length(jcondsmarks)
for (jcondmark in jcondsmarks){
# parallel::mclapply(jcondsmarks, function(jcondmark){
  
  jcond <- strsplit(jcondmark, "_")[[1]][[1]]
  jmark <- strsplit(jcondmark, "_")[[1]][[2]]
  print(paste("Running for:", jcond, jmark))
  
  # for loading GLMPCA output
  outbase <- paste0("PZ_", jmark, ".", jcond, ".", jsuffix, ".GLMPCA_var_correction.mergebinsize_", mergesize, ".binskeep_", jbins.keep, ".devmode.", jdate)
  # for writing PDF
  outpdf <- file.path(pdfdir, paste0(outbase, ".pdf"))
  if (file.exists(outpdf)){
    next
  }
  
  outf <- file.path(outdir, paste0(outbase, ".RData"))
  assertthat::assert_that(file.exists(outf))
  
  load(outf, v=T)
  
  infname <- paste0("lda_outputs.BM_", jmark, "_varfilt_countmat.", ldadate, ".", jcond, ".K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.", ldadate, ".", jcond, ".K-30.Robj")
  inf <- file.path(inmain, infname)
  assertthat::assert_that(file.exists(inf))
  
  load(inf, v=T)
  
  tm.result <- posterior(out.lda)
  
  dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
  
  jchromos <- paste("chr", c(seq(19)), sep = "")
  
  dat.var <- CalculateVarAll(dat.impute.log, jchromos)
  
  
  # Plot output -------------------------------------------------------------
  
  topics.mat <- glm.out$factors
  
  jsettings <- umap.defaults
  jsettings$n_neighbors <- 30
  jsettings$min_dist <- 0.1
  jsettings$random_state <- 123
  umap.out <- umap(topics.mat, config = jsettings)
  dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
  dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
  cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
  
  ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)
  
  dat.umap.long <- dat.umap.long %>%
    rowwise() %>%
    mutate(plate = ClipLast(as.character(cell), jsep = "_")) %>%
    left_join(., dat.var)
  
  m.louv <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) 
  
  m.louv.plate <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)  + 
    facet_wrap(~plate)
  
  m.var <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) 
  
  m.var.plate <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1)  + 
    facet_wrap(~plate)
  
  
  pdf(outpdf,useDingbats = FALSE)
  print(m.louv)
  print(m.louv.plate)
  print(m.var)
  print(m.var.plate)
  dev.off()
  
# }, mc.cores = ncores)
# for (jcondmark in jcondsmarks){
}
  
  

print(Sys.time() - jstart)


