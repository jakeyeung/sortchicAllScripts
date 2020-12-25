# Jake Yeung
# Date of Creation: 2020-12-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/3-run_DE_analysis_bins.more.R
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
jtype <- "bins"
# jdist <- "TES"

# outdir <- "/home/jyeung/data/from_rstudioserver/spikein_fits_BM_poisson"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3"
dir.create(outdir)

# jmark <- "H3K4me1"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmarks <- c("H3K9me3"); names(jmarks) <- jmarks


inf.k4me1 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.H3K4me1.txt")
inf.k9me3 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/BM_celltypes.bincutoff_0.binskeep_1000.byplate.szname_none.niter_500.reorder_rownames.dupfilt.2020-11-23.H3K9me3.txt")
assertthat::assert_that(file.exists(inf.k4me1))
assertthat::assert_that(file.exists(inf.k9me3))



inf.lst <- list(H3K4me1 = inf.k4me1,
                H3K9me3 = inf.k9me3)


for (jmark in jmarks){
  outf <- file.path(outdir, paste0("poisson_fit_", jtype, ".",  jmark, ".", Sys.Date(), ".newannot2.witherrors.MoreBins.RData"))
  if (file.exists(outf)){
    next
  }
  # assertthat::assert_that(!file.exists(outf))
  
  
  
  # indir <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.from_", jtype))
  # indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all.dbl_common_rows.cellfilt_binfilt")
  # # fname <- paste0("lda_outputs.count_mat_from_", jtype, ".", jmark, ".K-30.binarize.FALSE/ldaOut.count_mat_from_", jtype, ".", jmark, ".K-30.Robj")
  # fname <- paste0("lda_outputs.count_mat.", jmark, ".match_dbl.cellfilt.binfilt.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".match_dbl.cellfilt.binfilt.K-30.Robj")
  
  inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj"))
  
  # load(file.path(indir, fname), v=T)
  load(inf.lda, v=T)
  
  tm.result <- posterior(out.lda)
  
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)
  
  ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # Load metadata -----------------------------------------------------------
  
  # indir.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2")
  # inf.annot <- file.path(indir.annot, paste0("spikeins_mouse.BMround1and2_umaps_and_ratios.colfix.celltyping.2020-11-01.WithRelLevels.mark_", jmark, ".cell_cluster_tables.txt"))
  # indir.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs.merge_with_new_BM/clusterfilt.2020-11-04")
  # inf.annot <- file.path(indir.annot, paste0("cell_cluster_table.old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.txt"))
  inf.annot <- inf.lst[[jmark]]
  dat.annot <- fread(inf.annot) %>%
    rowwise() %>%
    dplyr::rename(plate.nbr = plate)
  
  # if (jmark == "H3K9me3"){
    dat.annot <- dat.annot %>%
      rowwise() %>%
      mutate(plate = ClipLast(x = cell,jsep = "_"))
  # }
  
  cells.keep <- colnames(count.mat)
  dat.annot.filt <- subset(dat.annot, cell %in% cells.keep)
  
  dat.umap.merge <- left_join(dat.umap, subset(dat.annot, select = c(cell, cluster)))
  
  cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
  ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() +  
    scale_color_manual(values = cbPalette) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
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


