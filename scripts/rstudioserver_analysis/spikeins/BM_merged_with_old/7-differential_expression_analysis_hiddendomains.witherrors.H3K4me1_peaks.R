# Jake Yeung
# Date of Creation: 2020-12-03
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/7-differential_expression_analysis_hiddendomains.witherrors.H3K4me1_peaks.R
# 




rm(list=ls())

library(parallel)
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

library(JFuncs)

jstart <- Sys.time()

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


FitGlmRowClustersPlate.withse <- function(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, jbin = NULL, returnobj=FALSE, with.se = FALSE){
  # use Offset by size of library
  # https://stats.stackexchange.com/questions/66791/where-does-the-offset-go-in-poisson-negative-binomial-regression
  # fit GLM for a row of a sparse matrix, should save some space?
  
  # pvalue by deviance goodness of fit: https://thestatsgeek.com/2014/04/26/deviance-goodness-of-fit-test-for-poisson-regression/
  # offset is in log because the model says the log counts is equal to RHS
  
  if (!is.null(nrow(jrow))){
    # probably a matrix of many rows, sum them up
    print(paste("Merging", nrow(jrow), "rows"))
    row <- Matrix::colSums(jrow)
  }
  dat <- data.frame(cell = cnames, ncuts = jrow, stringsAsFactors = FALSE) %>%
    left_join(., dat.annots.filt.mark, by = "cell") %>%
    left_join(., ncuts.cells.mark, by = "cell")
  
  # m1.pois <- glm(ncuts ~ 1 + Cluster + offset(ncuts.total), data = dat, family = "poisson")
  m1.pois <- glm(ncuts ~ 1 + Plate + Cluster + offset(log(ncuts.total)), data = dat, family = "poisson")
  mnull.pois <- glm(ncuts ~ 1 + Plate + offset(log(ncuts.total)), data = dat, family = "poisson")
  
  if (with.se){
    
  }
  
  if (!returnobj){
    jsum <- anova(mnull.pois, m1.pois)
    pval <- pchisq(jsum$Deviance[[2]], df = jsum$Df[[2]], lower.tail = FALSE)
    
    if (!with.se){
      out.dat <- data.frame(pval = pval, 
                            dev.diff = jsum$Deviance[[2]],
                            df.diff = jsum$Df[[2]],
                            t(as.data.frame(coefficients(m1.pois))), 
                            stringsAsFactors = FALSE)
    } else {
      estimates <- summary(m1.pois)$coefficients[, "Estimate"]
      names(estimates) <- make.names(paste(names(estimates), ".Estimate", sep = ""))
      stderrors <- summary(m1.pois)$coefficients[, "Std. Error"]
      names(stderrors) <- make.names(paste(names(stderrors), ".StdError", sep = ""))
      out.dat <- data.frame(pval = pval, 
                            dev.diff = jsum$Deviance[[2]],
                            df.diff = jsum$Df[[2]],
                            t(as.data.frame(c(estimates, stderrors))), 
                            stringsAsFactors = FALSE)
    }
    
    if (!is.null(jbin)){
      out.dat$bin <- jbin
      rownames(out.dat) <- jbin
    }
    return(out.dat)
  } else {
    return(list(fit.full = m1.pois, fit.null = mnull.pois, dat.input = dat))
  }
}


# Load LDA (contains countmat)  ---------------------------------------------------------------

ncores <- 8
hubprefix <- "/home/jyeung/hub_oudenaarden"
jtype <- "hiddendomains.H3K4me1_peaks"
# jdist <- "TES"

# outdir <- "/home/jyeung/data/from_rstudioserver/spikein_fits_BM_poisson"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3"

# jmark <- "H3K4me1"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmarks <- c("H3K9me3"); names(jmarks) <- jmarks

indir.counts <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/count_tables_from_peaks.from_sitecount_mat.H3K4me1_peaks")

for (jmark in jmarks){
  outf <- file.path(outdir, paste0("poisson_fit_", jtype, ".",  jmark, ".", Sys.Date(), ".newannot2.witherrors.RData"))
  if (file.exists(outf)){
    print(paste("File exists, skipping", outf))
    next
  }
  # assertthat::assert_that(!file.exists(outf))
  
  # load count matrix from raw
  inf.counts.lst <- list.files(indir.counts, pattern = paste0("BM_round1_round2_merged_", jmark, "_"), full.names = TRUE)
  mats.lst <- lapply(inf.counts.lst, function(inf.tmp){
    print(basename(inf.tmp))
    ReadMatTSSFormat(inf.tmp, add.coord = FALSE, as.sparse = TRUE, sort.rnames = TRUE)
  })
  all.rnames <- gtools::mixedsort(unlist(unique(lapply(mats.lst, rownames))))
  count.mat <- JFuncs::cbind.fill.lst(mats.lst = mats.lst, all.rnames = all.rnames, fill = 0)
  
  # indir <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.from_", jtype))
  # 
  # fname <- paste0("lda_outputs.count_mat_from_", jtype, ".", jmark, ".K-30.binarize.FALSE/ldaOut.count_mat_from_", jtype, ".", jmark, ".K-30.Robj")
  # 
  # load(file.path(indir, fname), v=T)
  # tm.result <- posterior(out.lda)
  # 
  # dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)
  # 
  # ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
  #   geom_point() + 
  #   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # Load metadata -----------------------------------------------------------
  
  # indir.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2")
  # inf.annot <- file.path(indir.annot, paste0("spikeins_mouse.BMround1and2_umaps_and_ratios.colfix.celltyping.2020-11-01.WithRelLevels.mark_", jmark, ".cell_cluster_tables.txt"))
  indir.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs.merge_with_new_BM/clusterfilt.2020-11-04")
  inf.annot <- file.path(indir.annot, paste0("cell_cluster_table.old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.txt"))
  dat.annot <- fread(inf.annot)
  
  if (jmark == "H3K9me3"){
    dat.annot <- dat.annot %>%
      rowwise() %>%
      mutate(plate = ClipLast(x = cell,jsep = "_"))
  }
  
  cells.keep <- colnames(count.mat)
  dat.annot.filt <- subset(dat.annot, cell %in% cells.keep)
  
  # dat.umap.merge <- left_join(dat.umap, subset(dat.annot, select = c(cell, cluster)))
  # ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = cluster)) + 
  #   geom_point() +  
  #   scale_color_manual(values = cbPalette) + 
  #   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
  # Load spikeins -----------------------------------------------------------
  
  
  # Run fits gene by gene ---------------------------------------------------
  
  cells.keep <- subset(dat.annot.filt, !is.na(cluster))$cell
  
  
  print(jmark)
  jmat.mark <- count.mat[, cells.keep]
  dat.annots.filt.mark <- dat.annot.filt %>%
    filter(cell %in% cells.keep) %>%
    mutate(Cluster = ifelse(cluster == "HSPCs", "aHSPCs", cluster)) %>%
    rowwise() %>%
    mutate(batch = IsRound1(cell, mark = jmark)) %>%
    mutate(Plate = ifelse(batch == "Round2", plate, "Round1"))
  
  print(unique(dat.annots.filt.mark$Plate))
  
  ncuts.for.fit.mark <- data.frame(cell = colnames(jmat.mark), ncuts.total = colSums(jmat.mark), stringsAsFactors = FALSE)
  cnames <- colnames(jmat.mark)
  
  jrow.names <- rownames(jmat.mark)
  names(jrow.names) <- jrow.names
  
  
  jfits.lst <- parallel::mclapply(jrow.names, function(jrow.name){
    jrow <- jmat.mark[jrow.name, ]
    jout <- FitGlmRowClustersPlate(jrow, cnames, dat.annots.filt.mark, ncuts.for.fit.mark, jbin = NULL, returnobj = FALSE, with.se = TRUE)
    return(jout)
  }, mc.cores = ncores)
  
  
  # Ssave outputs -----------------------------------------------------------
  # saveRDS(jfits.lst, outf)
  save(jfits.lst, dat.annots.filt.mark, ncuts.for.fit.mark, jmat.mark, file = outf)
  
  print(Sys.time() - jstart)
}

print("Done")
print(Sys.time() - jstart)


# 
# 
# # get error bars out of object
# jout <- FitGlmRowClustersPlate(jrow, cnames, dat.annots.filt.mark, ncuts.for.fit.mark, jbin = NULL, returnobj = FALSE)
# jout.obj <- FitGlmRowClustersPlate(jrow, cnames, dat.annots.filt.mark, ncuts.for.fit.mark, jbin = NULL, returnobj = TRUE)
# jout.obj.test <- FitGlmRowClustersPlate.withse(jrow, cnames, dat.annots.filt.mark, ncuts.for.fit.mark, jbin = NULL, returnobj = FALSE, with.se = TRUE)
# 
# x1 <- t(as.data.frame(coefficients(jout.obj$fit.full)))
# 
# x2 <- summary(jout.obj$fit.full)$coefficients[, "Estimate"]

