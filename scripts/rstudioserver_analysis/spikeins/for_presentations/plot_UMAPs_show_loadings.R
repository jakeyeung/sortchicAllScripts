# Jake Yeung
# Date of Creation: 2021-04-30
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/for_presentations/plot_UMAPs_show_loadings.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)


# Load scChIX -------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"


indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage")

dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  fread(file.path(indir.meta, fname)) %>%
    rowwise()
}) 

jmark <- "H3K4me1"
ggplot(dat.metas[[jmark]], aes(x = umap1, y = umap2)) + 
  geom_point(color = "grey85") + 
  theme_bw() + 
  theme_minimal() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# niter <- "100"
# inf.glmpca.lst <- lapply(jmarks, function(jmark){
#   if (jmark != "H3K27me3"){
#     inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/BM_from_matadj/glmpca.", jmark, ".from_matadj.platename_jrep.szname_none.niter_", niter, ".RData"))
#     assertthat::assert_that(file.exists(inf.glmpca))
#   } else {
#     inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/H3K27me3_rep2rep3reseq.peaks.varfilt/glmpca.", jmark, ".bincutoff_0.binskeep_0.platename_plate.szname_none.niter_500.reorder_rownames.dupfilt.suffix_peaks.RData"))
#   }
#   return(inf.glmpca)
# })
# 
# load(inf.glmpca.lst[[jmark]], v=T)
# 
# glm.outs <- lapply(inf.glmpca.lst, function(inf){
#   load(inf, v=T)
#   return(glm.out)
# })


# Load topic loadings -----------------------------------------------------
# 
# for (i in seq(ncol(glm.inits$U.init))){
#   print(i)
#   plot(density(glm.outs[[jmark]]$factors[, i]), main = i)
# }
# 
# for (i in seq(ncol(glm.inits$U.init))){
#   print(i)
#   plot(density(glm.inits$U.init[, i]), main = i)
# }
# 
# for (i in seq(ncol(glm.inits$U.init))){
#   print(i)
#   plot(density(glm.inits$U.init[, i]), main = i)
# }
# 
# 
# jdat <- data.frame(cell = rownames(glm.outs[[jmark]]$factors), glm.outs[[jmark]]$factors)
# dat.merge <- left_join(dat.metas[[jmark]], jdat)
# 
# i <- 4
# plot(density(glm.outs[[jmark]]$factors[, i]), main = i)
# ggplot(dat.merge, aes(x = umap1, y = umap2, color = dim4)) + 
#   geom_point() + 
#   theme_bw() + 
#   scale_color_viridis_c() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Load topic weights ------------------------------------------------------


# for K27me3 
inf.lda.lst <- lapply(jmarks, function(jmark){
  if (jmark != "H3K9me3"){
    # inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/tm_results_from_projs/tm_result_old_to_new.", jmark, ".2020-12-28.RData"))
    inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins.from_sitecount_mat.from_same_annot_file/lda_outputs.count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.K-30.binarize.FALSE/ldaOut.count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.K-30.Robj"))
  } else {
    # inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/tm_results_from_projs/tm_result_new_to_old.", jmark, ".2020-12-28.RData"))
    inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins/lda_outputs.count_mat_from_hiddendomains.", jmark, ".filtNAcells_allbins.K-30.binarize.FALSE/ldaOut.count_mat_from_hiddendomains.", jmark, ".filtNAcells_allbins.K-30.Robj"))
  }
  assertthat::assert_that(file.exists(inf.lda.tmp))
  
  return(inf.lda.tmp)
})

out.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.lda <- inf.lda.lst[[jmark]]
  load(inf.lda, v=T)  # out.lda, count.mat
  tm.result <- posterior(out.lda)
  return(list(tm.result = tm.result, count.mat = count.mat))
})

# for (i in 1:30){
#   x <- out.lst$H3K4me1$tm.result$topics[, i]
#   plot(density(log(x / (1 - x))), main = i)
# }

Xraw <- out.lst$H3K4me1$tm.result$topic
colnames(Xraw) <- paste("topic", colnames(Xraw), sep = "_")
X <- log2(Xraw / (1 - Xraw))
dat.lda <- data.frame(cell = rownames(X), X, stringsAsFactors = FALSE)
dat.lda.raw <- data.frame(cell = rownames(Xraw), Xraw, stringsAsFactors = FALSE)
dat.merge.lda <- left_join(dat.metas[[jmark]], dat.lda)
dat.merge.lda.raw <- left_join(dat.metas[[jmark]], dat.lda.raw)

outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/for_presentations/H3K4me1_UMAP_fitting_Gauss_mix.", Sys.Date(), ".pdf")

pdf(outpdf, useDingbats = FALSE)
jtop <- "topic_8"


for (jmarktmp in jmarks){
  m <- ggplot(dat.metas[[jmarktmp]], aes(x = umap1, y = umap2)) + 
    geom_point(color = "grey85") + 
    theme_bw() + 
    theme_minimal() + 
    ggtitle(jmarktmp) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}

ggplot(dat.merge.lda.raw, aes_string(x = "umap1", y = "umap2", color = jtop)) + 
  geom_point() + 
  theme_minimal() + 
  ggtitle(jtop) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge.lda.raw, aes_string(x = "umap1", y = "umap2", color = jtop)) + 
  geom_point() + 
  theme_minimal() + 
  ggtitle(jtop) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


ggplot(dat.merge.lda, aes_string(x = jtop)) + 
  geom_density() + 
  theme_bw() + 
  ggtitle(jtop) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# fit 
library(mixtools)

x <- dat.merge.lda[[jtop]]

jthres <- 0.5
mm <- mixtools::normalmixEM(x = x, lambda = c(0.7, 0.3), mu = c(-6, -1), sigma = c(1, 1), k = 2)
xline <- min(mm$x[which(mm$posterior[, 1] < jthres)])

par(mfrow=c(1,1), mar=c(12, 4.1, 12, 2.1), mgp=c(3, 1, 0), las=0)
plot(mm, whichplots = 2, xlab2 = "log2(th / (1 - th))", main2 = paste("Distribution of", jtop, "scores across cells"))
abline(v = xline)

cells.keep <- subset(dat.merge.lda, cluster == "pDCs")$cell
ggplot(dat.merge.lda %>% arrange(topic_8), aes(x = umap1, y = umap2, color = cell %in% cells.keep)) + 
  geom_point() + 
  ggtitle("pDCs") + 
  theme_minimal() + 
  scale_color_manual(values = c("#42245D", "#32CD32")) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
cells.keep <- subset(dat.merge.lda, topic_8 >= xline)$cell
ggplot(dat.merge.lda %>% arrange(topic_8), aes(x = umap1, y = umap2, color = cell %in% cells.keep)) + 
  geom_point() + 
  ggtitle("Topic8") + 
  theme_minimal() + 
  scale_color_manual(values = c("#42245D", "#32CD32")) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


dev.off()


# 
# ggplot(dat.merge.lda, aes(x = umap1, y = umap2, color = topic_8)) + 
#   geom_point() + 
#   theme_bw() + 
#   scale_color_viridis_c() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# ggplot(dat.merge.lda, aes(x = umap1, y = umap2, color = X8)) + 
#   geom_point() + 
#   theme_bw() + 
#   scale_color_viridis_c() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
