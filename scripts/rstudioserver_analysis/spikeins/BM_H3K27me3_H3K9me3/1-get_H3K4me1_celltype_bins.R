# Jake Yeung
# Date of Creation: 2021-01-18
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_H3K9me3/1-get_H3K4me1_celltype_bins.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)
library(topicmodels)

make.topics.plots <- TRUE
hubprefix <- "/home/jyeung/hub_oudenaarden"

# jmarks <- c("H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

jmark1 <- "H3K4me1"
jmark2 <- "H3K9me3"

# Load .metasmetas  -------------------------------------------------------------

# indir.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned"
indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage")

dat.metas <- lapply(jmarks, function(jmark){
  print(jmark)
  # fname.tmp <- paste0("cell_cluster_table_with_spikeins.", jmark, ".2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt")
  fname.tmp <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  inf <- file.path(indir.meta, fname.tmp)
  jdat <- fread(inf)
  if (jmark == "H3K9me3"){
    jdat$mark <- jmark
  }
  return(jdat)
})


# Load LDA ----------------------------------------------------------------

# jmark <- "H3K27me3"
# jmark2 <- "H3K9me3"


outs.all.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  if (jmark != "H3K27me3"){
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj"))
  } else {
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.Robj"))
  }
  load(inf.lda, v=T)
  return(list(out.lda = out.lda, count.mat = count.mat))
})

out.lda.lst <- lapply(outs.all.lst, function(jout) jout$out.lda)


count.mat.lst <- lapply(outs.all.lst, function(jout){
  jmat <- jout$count.mat
  jmat <- sweep(jmat, MARGIN = 2, STATS = colSums(jmat), FUN = "/")
  # jmat <- BinarizeMatrix(jmat)
})

tm.result.lst <- lapply(out.lda.lst, function(jout) AddTopicToTmResult(posterior(jout)))

# order topics: K27me3 only 
topics.ordered <- OrderTopicsByEntropy(tm.result.lst$H3K27me3)

# Get topics --------------------------------------------------------------

# write 
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K27me3_H3K9me3_analysis"
pdfname <- paste0(jmark1, "_50kb_topics.pdf")

tm.result.lst[[jmark1]]$terms[1:5, 1:5]

jtopics <- topics.ordered$topic
# jmark <- "H3K27me3"


if (make.topics.plots){
  pdf(file.path(outdir, pdfname), useDingbats = FALSE)
}

for (jtopic in jtopics){
  jvec <- tm.result.lst[[jmark1]]$topics[, jtopic]
  jdat <- data.frame(loading = jvec, cell = names(jvec), stringsAsFactors = FALSE)
  dat.metas.tmp <- left_join(dat.metas[[jmark1]], jdat)
  m <- ggplot(dat.metas.tmp, mapping = aes(x = umap1, y = umap2, color = loading)) + 
    geom_point() + 
    ggtitle(jtopic) + 
    theme_minimal() + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}

if (make.topics.plots){
  dev.off() 
}

keeptopn <- 150

# 
# # Define celltype specific topics -----------------------------------------
# 

jtopics.lst <- list("topic7" = "DCs",
                    "topic10" = "Eryths",
                    "topic12" = "Eryths2",
                    "topic20" = "GranulocytesFirst",
                    "topic1" = "Granulocytes",
                    "topic5" = "Basophils",
                    "topic30" = "BcellsSmaller",
                    "topic18" = "BcellsLarger",
                    "topic21" = "HSPCslater",
                    "topic17" = "HSPCs",
                    "topic3" = "NKs",
                    "topic26" = "pDCs"
                    )

jtopics2.lst <- names(jtopics.lst)
names(jtopics2.lst) <- unlist(jtopics.lst)

# get topics
jbins.lst <- lapply(jtopics2.lst, function(jtopic){
  jvec <- sort(tm.result.lst[[jmark1]]$terms[jtopic, ], decreasing = TRUE)
  names(jvec[1:keeptopn])
})



# Order metas  ------------------------------------------------------------

ctypes <- c("Eryths", "Bcells", "NKs", "pDCs", "Granulocytes", "DCs", "Basophils", "HSPCs")
dat.metas.reordered <- lapply(jmarks, function(jmark){
  if (jmark != "H3K9me3"){
    jdat.reordered <- dat.metas[[jmark]] 
    jdat.reordered$cluster <- factor(jdat.reordered$cluster, levels = ctypes)
    jdat.reordered <- jdat.reordered %>%
      arrange(cluster)
  } else {
    jdat.reordered <- dat.metas[[jmark]]
  }
  return(jdat.reordered)
})

# order bins
ctypes.bins <- c("Eryths", "BcellsLarger", "NKs", "pDCs", "Granulocytes", "DCs", "Basophils", "HSPCs"); names(ctypes.bins) <- ctypes.bins

jbins.lst.filt <- jbins.lst[ctypes.bins]


# Write bins to output ----------------------------------------------------

head(ctypes.bins)

outbins <- file.path(outdir, paste0(jmark1, "_cluster_specific_bins_50kb.txt"))

jbins.filt.dat <- lapply(ctypes.bins, function(jctype){
  jbins.tmp <- jbins.lst.filt[[jctype]]
  jchr <- sapply(jbins.tmp, GetChromo)
  jstart <- sapply(jbins.tmp, GetStart)
  jend <- sapply(jbins.tmp, GetEnd)
  dat <- data.frame(Chr = jchr, Start = jstart, End = jend, bin = jbins.tmp, ctype = jctype, stringsAsFactors = FALSE)
}) %>%
  bind_rows()

if (make.topics.plots){
  fwrite(jbins.filt.dat, file = outbins, sep = "\t")
}


# Show heatmap  -----------------------------------------------------------

dat.imputed.lst <- lapply(jmarks, function(jmark){
  jmat <- t(log2(tm.result.lst[[jmark]]$topic %*% tm.result.lst[[jmark]]$term))
  # reorder
  cells.ordered <- dat.metas.reordered[[jmark]]$cell
  jmat <- jmat[, cells.ordered]
})

cells.ordered.lst <- lapply(jmarks, function(jmark){
  dat.metas.reordered[[jmark]]$cell
})

# heatmap(dat.imputed.lst$H3K27me3[unlist(jbins.lst), ], )
bins.keep <- unlist(jbins.lst.filt)

pdfname.heatmap <- paste0(jmark1, "_50kb_topics.heatmap.pdf")
# pdf(file.path(outdir, pdfname.heatmap), useDingbats = FALSE)


# do pseudogenes
jmark.ref <- "H3K4me1"

colsidecolors <- dat.metas.reordered[[jmark.ref]]$colorcode
heatmap3::heatmap3(dat.imputed.lst[[jmark.ref]][bins.keep, ], labCol = "", Rowv = NA, Colv = NA, ColSideColors = colsidecolors, scale = "row", revC = TRUE, main = paste0(jmark.ref))

jmat.tmp <- count.mat.lst[[jmark.ref]][, cells.ordered.lst[[jmark.ref]]]
# cells.tmp <- which(colnames(jmat.tmp) %in% cells.ordered.lst[[jmark.ref]])
# cells.tmp.names <- colnames(jmat.tmp)[cells.tmp]

mat.pseudogenes <- lapply(ctypes.bins, function(jctype){
  jbins.tmp <- which(rownames(jmat.tmp) %in% jbins.lst.filt[[jctype]])
  colMeans(as.matrix(jmat.tmp[jbins.tmp, ]))
}) %>%
  as.data.frame() %>%
  bind_rows() %>%
  t()
colnames(mat.pseudogenes) <- cells.tmp.names

# heatmap3::heatmap3(as.matrix(count.mat.lst[[jmark]][bins.keep, cells.ordered.lst[[jmark]]]), labCol = "", Rowv = NA, Colv = NA, ColSideColors = colsidecolors, scale = "row", revC = TRUE, main = paste0(jmark))
# heatmap3::heatmap3(log2(mat.pseudogenes + 1), labCol = "", Rowv = NA, Colv = NA, ColSideColors = colsidecolors, scale = "row", revC = TRUE, main = paste0(jmark))
heatmap3::heatmap3(mat.pseudogenes, labCol = "", Rowv = NA, Colv = NA, ColSideColors = colsidecolors, scale = "row", revC = TRUE, main = paste0(jmark.ref))

# dev.off()


# Compare with K9me3 bins  ------------------------------------------------

# jmark.ref <- "H3K9me3"
# jmark.ref <- "H3K4me3"
# jmark.ref <- "H3K4me1"
# jmark.ref <- "H3K4me3"
# jmark.ref <- "H3K27me3"
# jmark.ref <- "H3K9me3"

dat.metas.tmp <- dat.metas.reordered[[jmark.ref]] %>%
  filter(cluster != "HSPCs")

cells.keep <- dat.metas.tmp$cell
  
colsidecolors <- dat.metas.tmp$clustercol
print(head(colsidecolors))
bins.keep2 <- rownames(dat.imputed.lst[[jmark.ref]])[rownames(dat.imputed.lst[[jmark.ref]]) %in% bins.keep]
bins.keep2.reordered <- bins.keep[bins.keep %in% bins.keep2]
print(length(bins.keep2.reordered))

# colsidecolors <- dat.metas.reordered[[jmark.ref]]$clustercol
# print(head(colsidecolors))
# bins.keep2 <- rownames(dat.imputed.lst[[jmark.ref]])[rownames(dat.imputed.lst[[jmark.ref]]) %in% bins.keep]
# bins.keep2.reordered <- bins.keep[bins.keep %in% bins.keep2]
# print(length(bins.keep2.reordered))

heatmap3::heatmap3(dat.imputed.lst[[jmark.ref]][bins.keep2.reordered, cells.keep], labCol = "", Rowv = NA, Colv = NA, ColSideColors = colsidecolors, scale = "row", revC = TRUE, main = paste0(jmark.ref))

jmat.tmp <- count.mat.lst[[jmark.ref]][, cells.keep]
cells.tmp <- which(colnames(jmat.tmp) %in% cells.keep)
cells.tmp.names <- colnames(jmat.tmp)[cells.tmp]

mat.pseudogenes <- lapply(ctypes.bins, function(jctype){
  jbins.tmp <- which(rownames(jmat.tmp) %in% jbins.lst.filt[[jctype]])
  colMeans(as.matrix(jmat.tmp[jbins.tmp, cells.tmp.names]))
}) %>%
  as.data.frame() %>%
  bind_rows() %>%
  t()
colnames(mat.pseudogenes) <- cells.tmp.names

# heatmap3::heatmap3(log2(mat.pseudogenes + 1), labCol = "", Rowv = NA, Colv = NA, ColSideColors = colsidecolors, scale = "row", revC = TRUE, main = paste0(jmark))
heatmap3::heatmap3(mat.pseudogenes, labCol = "", Rowv = NA, Colv = NA, ColSideColors = colsidecolors, scale = "row", revC = TRUE, main = paste0(jmark.ref))

# 
# # Get K9me3 counts at these bins ------------------------------------------------
# 
# 
# 
# # H3K9me3 cluster-specific pseduboulks
# 
# dat.metas.k9me3.round2only <- subset(dat.metas$H3K9me3, jrep == "rep1old") 
# cnames.keep.lst <- split(x = dat.metas.k9me3.round2only$cell, f = dat.metas.k9me3.round2only$cluster)
# pbulks.lst <- SumAcrossClusters(outs.all.lst[[jmark2]]$count.mat, cnames.keep.lst = cnames.keep.lst)
# 
# mat.pbulk <- bind_rows(pbulks.lst) %>%
#   as.data.frame()
# rownames(mat.pbulk) <- rownames(count.mat.lst[[jmark2]])
# mat.pbulk <- sweep(mat.pbulk, MARGIN = 2, STATS = colSums(mat.pbulk), FUN = "/")
# 
# dat.pbulk <- as.matrix(mat.pbulk) %>%
#   melt()
# colnames(dat.pbulk) <- c("bin", "ctype", "cuts")
# 
# ggplot(dat.pbulk, aes(x = log2(cuts), fill = ctype)) + 
#   geom_density(alpha = 0.25) + 
#   theme_bw() + 
#   facet_wrap(~ctype) + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# pbulks.avg.filt.byctype <- lapply(jbins.lst.filt, function(jbins){
#   subset(dat.pbulk, bin %in% jbins) %>%
#     group_by(ctype) %>%
#     summarise(cuts.log2 = mean(log2(cuts)))
# })
# 
# print(pbulks.avg.filt.byctype)
# 
# cnames.keep.lst <- split(x = dat.metas[[jmark1]]$cell, f = dat.metas[[jmark1]]$cluster)
# pbulks.lst <- SumAcrossClusters(outs.all.lst[[jmark1]]$count.mat, cnames.keep.lst = cnames.keep.lst)
# 
# mat.pbulk <- bind_rows(pbulks.lst) %>%
#   as.data.frame()
# rownames(mat.pbulk) <- rownames(count.mat.lst[[jmark1]])
# mat.pbulk <- sweep(mat.pbulk, MARGIN = 2, STATS = colSums(mat.pbulk), FUN = "/")
# 
# dat.pbulk <- as.matrix(mat.pbulk) %>%
#   melt()
# colnames(dat.pbulk) <- c("bin", "ctype", "cuts")
# 
# ggplot(dat.pbulk, aes(x = log2(cuts), fill = ctype)) + 
#   geom_density(alpha = 0.25) + 
#   theme_bw() + 
#   facet_wrap(~ctype) + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# # 
# # jmat.pbulk.long$gene <- names(jmat.pbulk[[1]])
# # jmat.pbulk.long <- jmat.pbulk.long %>%
# #   melt()
# # colnames(jmat.pbulk.long) <- c("gene", "ctype", "cuts")
# # 
# # pbulk.long <- do.call()
# # 
# 
# 
# 
# 
# 
# 
