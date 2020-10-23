# Jake Yeung
# Date of Creation: 2020-10-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged/1-LDA_downstream.filter_var.stainlater.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)
library(scchicFuncs)

library(topicmodels)

library(ggrepel)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)


SplitGetLast <- function(x, jsplit = "-"){
  # get last element after splitting
  xsplit <- strsplit(x, jsplit)[[1]]
  xnew <- xsplit[[length(xsplit)]]
  return(xnew)
}



# Functions ---------------------------------------------------------------




# Constnats ---------------------------------------------------------------

inf.spikeins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all.stainlater/spikein_info_BM_round2_all.txt"
assertthat::assert_that(file.exists(inf.spikeins))

lowvar.filt <- 1.5
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all.stainlater/varfilt"
dir.create(outdir)
outpdf <- file.path(outdir, "plot_varfilt_UMAP_before_after.stainlater.pdf")

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#28f9ff", "#88497e", "#bcf5c3", "#86f115")

jmarks <- c("H3K27me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_mouse.BMround2all"
# dir.create(outdir)

# inmain <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_VAN5046")
inmain <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all.stainlater")
assertthat::assert_that(dir.exists(inmain))

# load(inf.lda, v=T)

jsuffix <- ".filt_0.15_0.95_counts_and_l2r.K-30"


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

out.lst <- lapply(jmarks, function(jmark){
  # dname <- paste0("lda_outputs.count_mat_", jmark, jsuffix, ".binarize.FALSE")
  dname <- paste0("lda_outputs.count_mat.", jmark, jsuffix, ".binarize.FALSE")
  indir <- file.path(inmain, dname)
  # fname <- paste0("ldaOut.count_mat_", jmark, "_counts_filt.2020-09-12.K-30.Robj")
  fname <- paste0("ldaOut.count_mat.", jmark, jsuffix, ".Robj")
  inf.lda <- file.path(indir, fname)
  assertthat::assert_that(file.exists(inf.lda))
  print(inf.lda)
  load(inf.lda, v=T)
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result)
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)
  dat.umap$mark <- jmark
  return(list(dat.umap = dat.umap, tm.result = tm.result, count.mat = count.mat))
})

dat.umaps.lst <- lapply(out.lst, function(lst){
  return(lst$dat.umap)
})

tm.result.lst <- lapply(out.lst, function(lst){
  return(lst$tm.result)
})

count.mat.lst <- lapply(out.lst, function(lst){
  return(lst$count.mat)
})

# 
# tm.result.lst <- lapply(jmarks, function(jmark){
#   dname <- paste0("lda_outputs.count_mat.", jmark, jsuffix, ".binarize.FALSE")
#   indir <- file.path(inmain, dname)
#   fname <- paste0("ldaOut.count_mat.", jmark, jsuffix, ".Robj")
#   inf.lda <- file.path(indir, fname)
#   assertthat::assert_that(file.exists(inf.lda))
#   print(inf.lda)
#   load(inf.lda, v=T)
#   tm.result <- posterior(out.lda)
#   tm.result <- AddTopicToTmResult(tm.result)
#   return(tm.result)
# })



for (jmark in jmarks){
  m <- ggplot(dat.umaps.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette) + 
    ggtitle(jmark)
  print(m)
}


dat.umaps.long <- bind_rows(dat.umaps.lst) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "_"), 
         plate = as.numeric(SplitGetLast(experi, jsplit = "-")),
         rowcoord = AddPlateCoordinates(cell)$rowcoord,
         colcoord = AddPlateCoordinates(cell)$colcoord,
         jrep = GetRepBM(experiname = experi), 
         stype = AnnotateSortFromLayoutBMall(plate = plate, rowcoord = rowcoord, colcoord = colcoord, jrep = jrep, jmark = mark))

dat.meta <- subset(dat.umaps.long, select = c(cell, experi, plate, rowcoord, colcoord, stype))

# for (jmark in jmarks){
#   m <- ggplot(dat.umaps.long %>% filter(mark == jmark), aes(x = umap1, y = umap2, color = stype)) + 
#     geom_point() + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     scale_color_manual(values = cbPalette) + 
#     ggtitle(jmark)
#   print(m)
# }




# Variance stuff ----------------------------------------------------------

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var.lst <- lapply(jmarks, function(jmark){
  tm.result <- tm.result.lst[[jmark]]
  dat.impute <- t(log2(tm.result$topics %*% tm.result$terms))
  dat.var <- CalculateVarAll(dat.impute, jchromos)
}) 

dat.var <- dat.var.lst %>%
  bind_rows()

dat.merge <- left_join(dat.umaps.long, dat.var)

# for (jmark in jmarks){
#   m <- ggplot(dat.merge %>% filter(mark == jmark), aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
#     geom_point() + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     scale_color_viridis_c(direction = -1) + 
#     ggtitle(jmark)
#   print(m)
# }


# Add log2ratio -----------------------------------------------------------



dat.spikeins.mat <- fread(inf.spikeins) %>%
  rowwise() %>%
  mutate(plate = as.numeric(SplitGetLast(experi, jsplit = "-")),
         jrep = GetRepBM(experiname = experi), 
         stype = AnnotateSortFromLayoutBMall(plate = plate, rowcoord = rowcoord, colcoord = colcoord, jrep = jrep, jmark = mark))


jsub <- subset(dat.spikeins.mat, select = c(samp, chromocounts, spikeincounts)) %>%
  rowwise() %>%
  mutate(l2r = log2(chromocounts / spikeincounts))

dat.merge <- left_join(dat.merge, jsub, by = c("cell" = "samp"))



# # Annotate from Giladi ----------------------------------------------------
# 
# 
# jinf.tss <- file.path(hubprefix, "jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed")
# assertthat::assert_that(file.exists(jinf.tss))
# 
# inf.annot <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Giladi_et_al_2018/giladi_pseudobulk_datsum_and_ncells.DESeq_and_QuantNorm.RData")
# load(inf.annot, v=T)
# 
# dat.sum.long <- data.table::melt(dat.sum.norm.quantnorm)
# colnames(dat.sum.long) <- c("gene", "celltype", "exprs")
# 
# dat.sum.long <- dat.sum.long %>%
#   group_by(gene) %>%
#   mutate(zscore = scale(exprs, center = TRUE, scale = TRUE)) %>%
#   filter(!is.nan(zscore))
# 
# topn <- 150
# 

# Filter low var ----------------------------------------------------------

dat.merge <- dat.merge %>%
  rowwise() %>%
  mutate(is.highvar = cell.var.within.sum.norm >= lowvar.filt)

cells.keep <- subset(dat.merge, is.highvar)$cell


# Write output ------------------------------------------------------------

for (jmark in jmarks){
  jmat.tmp <- count.mat.lst[[jmark]]
  cols.keep <- colnames(jmat.tmp) %in% cells.keep
  print(dim(jmat.tmp))
  jmat.tmp.filt <- jmat.tmp[, cols.keep]
  print(dim(jmat.tmp.filt))
  outf <- file.path(outdir, paste0("count_mat.", jmark, ".varfilt_", lowvar.filt, ".rds"))
  saveRDS(jmat.tmp.filt, file = outf)
}


# Check k4me1 -------------------------------------------------------------

pdf(outpdf, useDingbats = FALSE)

# jmark <- "H3K4me1"
for (jmark in jmarks){
  m <- ggplot(dat.merge %>% filter(mark == jmark), 
              aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + 
    scale_color_viridis_c(direction = -1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark, "before filt")
  print(m)
  m <- ggplot(dat.merge %>% filter(mark == jmark & is.highvar), 
              aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + 
    scale_color_viridis_c(direction = -1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark, "after filt")
  print(m)
}

dev.off()

