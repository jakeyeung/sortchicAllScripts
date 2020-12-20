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



# Load LDA (contains countmat)  ---------------------------------------------------------------

ncores <- 8
hubprefix <- "/home/jyeung/hub_oudenaarden"
# jtype <- "bins"
# jdist <- "TSS"
jtype <- "TSS"

# outdir <- "/home/jyeung/data/from_rstudioserver/spikein_fits_BM_poisson"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.again"
dir.create(outdir)

# jmark <- "H3K4me1"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmarks <- c("H3K9me3"); names(jmarks) <- jmarks


inf.k4me1 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.H3K4me1.txt")
inf.k4me3 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.H3K4me3.txt")
inf.k27me3 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.H3K27me3_techrepmerged/BM_rep2_rep3reseq_H3K27me3.2020-12-10.txt")
inf.k9me3 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/BM_celltypes.bincutoff_0.binskeep_1000.byplate.szname_none.niter_500.reorder_rownames.dupfilt.2020-11-23.H3K9me3.txt")
assertthat::assert_that(file.exists(inf.k4me1))
assertthat::assert_that(file.exists(inf.k4me3))
assertthat::assert_that(file.exists(inf.k27me3))
assertthat::assert_that(file.exists(inf.k9me3))

inf.lst <- list(H3K4me1 = inf.k4me1,
                H3K4me3 = inf.k4me3,
                H3K27me3 = inf.k27me3,
                H3K9me3 = inf.k9me3)

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

for (jmark in jmarks){
  outf <- file.path(outdir, paste0("poisson_fit_", jtype, ".",  jmark, ".", Sys.Date(), ".newannot2.witherrors.TSS.RData"))
  if (file.exists(outf)){
    print(paste("File exists, skipping", outf))
    next
  }
  
  if (jmark != "H3K27me3"){
    # inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj"))
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.from_TSS.dist_10000.RemoveDupRows/lda_outputs.count_mat_from_TSS.", jmark, ".dist_10000.K-30.binarize.FALSE/ldaOut.count_mat_from_TSS.", jmark, ".dist_10000.K-30.Robj"))
  } else {
    # inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.Robj"))
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_TSS/lda_outputs.PZ-BM-rep3-H3K27me3-rep2rep3reseq.TSS.varfilt.K-30.binarize.FALSE/ldaOut.PZ-BM-rep3-H3K27me3-rep2rep3reseq.TSS.varfilt.K-30.Robj"))
  }
  load(inf.lda, v=T)
  
  tm.result <- posterior(out.lda)
  
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)
  
  ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # Load metadata -----------------------------------------------------------
  
  inf.annot <- inf.lst[[jmark]]
  dat.annot <- fread(inf.annot) %>%
    rowwise() %>%
    dplyr::rename(plate.nbr = plate)
  
  if (jmark == "H3K27me3"){
    # use cluster.fewer
    dat.annot <- dat.annot %>%
      dplyr::rename(cluster.more = cluster,
                    cluster = cluster.fewer)
  }
  
    dat.annot <- dat.annot %>%
      rowwise() %>%
      mutate(plate = ClipLast(x = cell,jsep = "_"))
  
  cells.keep <- colnames(count.mat)
  dat.annot.filt <- subset(dat.annot, cell %in% cells.keep)
  
  dat.umap.merge <- left_join(dat.umap, subset(dat.annot, select = c(cell, cluster)))
  
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


