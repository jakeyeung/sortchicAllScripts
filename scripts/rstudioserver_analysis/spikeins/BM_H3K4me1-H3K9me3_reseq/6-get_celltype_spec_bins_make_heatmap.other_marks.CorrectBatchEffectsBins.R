# Jake Yeung
# Date of Creation: 2020-12-19
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/6-get_celltype_spec_bins_make_heatmap.other_marks.CorrectBatchEffectsBins.R
# description

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(topicmodels)

library(JFuncs)

library(hash)
library(igraph)
library(umap)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(DescTools)

library(topicmodels)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

ncores <- 4
hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

outrdata2 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_H3K4me1_H3K9me3_analysis/batch_corrected_imputed_values.bins.all_marks.mat.namesfix.", Sys.Date(), ".H3K27me3rep2rep3reseq.RData")
assertthat::assert_that(!file.exists(outrdata2))

# Load LDA outputs --------------------------------------------------------

# for K27me3 
inf.lda.lst <- lapply(jmarks, function(jmark){
  if (jmark != "H3K27me3"){
    inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj"))
  } else {
    inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_", jmark, ".cleaned.varfilt_2.K-30.Robj"))
  }
  return(inf.lda.tmp)
})


out.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.lda <- inf.lda.lst[[jmark]]
  # inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj"))
  load(inf.lda, v=T)
  tm.result <- posterior(out.lda)
  return(list(tm.result = tm.result, count.mat = count.mat))
})


# Load meta data  ---------------------------------------------------------

dat.metas <- lapply(jmarks, function(jmark){
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/umaps_final/umaps_final_get_marks.", jmark, ".metadata.2020-12-17.txt")
  fread(inf)
})

# add jrep2 for batch correction?

dat.metas <- lapply(jmarks, function(jmark){
  dat.metas.tmp <- dat.metas[[jmark]]
  if (jmark != "H3K9me3"){
    dat.metas.tmp$jrep2 <- sapply(dat.metas.tmp$jrep, function(x) ifelse(x == "rep1old", "zold", "anew"))
  } else {
    dat.metas.tmp$jrep2 <- sapply(dat.metas.tmp$jrep, function(x) ifelse(x != "rep1old", "zold", "anew"))
  }
  return(dat.metas.tmp)
})
# dat.metas$H3K4me1$jrep2 <- sapply(dat.metas$H3K4me1$jrep, function(x) ifelse(x == "rep1old", "zold", "anew"))  # batch2 is better than batch1
# dat.metas$H3K9me3$jrep2 <- sapply(dat.metas$H3K9me3$jrep, function(x) ifelse(x != "rep1old", "zold", "anew"))  # batch1 is better than batch2

# Select bins  ------------------------------------------------------------

load("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_H3K4me1_H3K9me3_analysis/de_analysis_H3K4me1_H3K9me3.RData", v=T)

k9.bins <- names(which(pvals.lst2 < 1e-10))
k4.bins <- names(which(pvals.lst1 < 1e-100))


# Select bins and correct -------------------------------------------------

imputed.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  log2(t(out.lst[[jmark]]$tm.result$topics %*% out.lst[[jmark]]$tm.result$terms))
})

imputed.long.lst <- lapply(jmarks, function(jmark){
  rnames.keep <- which(rownames(imputed.lst[[jmark]]) %in% k9.bins)
  jmat.filt <- imputed.lst[[jmark]][rnames.keep, ] %>%
    data.table::melt()
  colnames(jmat.filt) <- c("rname", "cell", "log2exprs")
  jmat.filt <- jmat.filt %>%
    left_join(., dat.metas[[jmark]])
  return(jmat.filt)
})


# Correct batch  ----------------------------------------------------------

print("Correcting batch multicore 4")
system.time(
  # dat.adj.lst <- lapply(imputed.long.lst, function(jdat){
  dat.adj.lst <- mclapply(jmarks, function(jmark){
    jdat <- imputed.long.lst[[jmark]]
    if (jmark != "H3K27me3"){
      dat.adj <- jdat %>%
        group_by(rname) %>%
        do(AdjustBatchEffect(.))
    } else {
      dat.adj <- jdat
      dat.adj$plateadj2 <- 0
      dat.adj$clstradj2 <- 0
      dat.adj$log2exprsadj <- dat.adj$log2exprs
    }
    return(dat.adj)
  }, mc.cores = ncores)
)

dat.adj.lst2 <- lapply(dat.adj.lst, function(jdat){
  subset(jdat, select = c(rname, cell, log2exprs, cluster, batch, jrep, jrep2, plateadj2, clstradj2, log2exprsadj))
})
# 
# outrdata <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_H3K4me1_H3K9me3_analysis/batch_corrected_imputed_values_H3K4me1_H3K9me3.RData"
# save(dat.adj.lst2, file = outrdata)


mat.adj.lst <- lapply(dat.adj.lst2, function(dat.adj){
  mat.adj <- data.table::dcast(dat.adj, formula = rname ~ cell, value.var = "log2exprsadj")
})

save(mat.adj.lst, file = outrdata2)


# 
# jtest2 <- subset(imputed.long.lst$H3K9me3, rname == k9.bins[[1]])
# jtest1 <- subset(imputed.long.lst$H3K4me1, rname == imputed.long.lst$H3K4me1$rname[[1]])
# 
# 
# jtest2.adj <- AdjustBatchEffect(jtest2)
# jtest1.adj <- AdjustBatchEffect(jtest1)
# 
# m1.orig <- ggplot(jtest1.adj, aes(x = cluster, y = log2exprs)) + 
#   geom_point() + 
#   facet_wrap(~jrep, ncol = 1) + 
#   geom_boxplot() 
# 
# m1.after <- ggplot(jtest1.adj, aes(x = cluster, y = log2exprsadj)) + 
#   geom_point() + 
#   facet_wrap(~jrep, ncol = 1) + 
#   geom_boxplot() 
# 
# JFuncs::multiplot(m1.orig, m1.after, cols = 2)
# 
# m2.orig <- ggplot(jtest2.adj, aes(x = cluster, y = log2exprs)) + 
#   geom_point() + 
#   facet_wrap(~jrep, ncol = 1) + 
#   geom_boxplot() 
# 
# m2.after <- ggplot(jtest2.adj, aes(x = cluster, y = log2exprsadj)) + 
#   geom_point() + 
#   facet_wrap(~jrep, ncol = 1) + 
#   geom_boxplot() 
# 
# JFuncs::multiplot(m2.orig, m2.after, cols = 2)
