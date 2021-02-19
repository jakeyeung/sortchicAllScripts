# Jake Yeung
# Date of Creation: 2021-02-19
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/19-make_heatmap_K4me1_and_K9me3_dynamic_bins.from_LDA.subset_regions.R
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
library(DescTools)
library(hash)




GetPbulks <- function(dat.metas, jmarks, jsuffixmain, scale.factor = 1000000, merge.ctypes.by.lineage = FALSE){
  out.lst <- lapply(jmarks, function(jmark){
    # jsuffixmain <- "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table.H3K9me3.2021-02-15.txt"
    
    if (jsuffixmain == "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table"){
      jsuffix <- paste0(jsuffixmain, ".", jmark, ".2021-02-15.txt")
    } else {
      jsuffix <- jsuffixmain
      # jsuffix <- paste0("celltype_specific_genes.TSS_10kb.txt")
    }
    print(jsuffix)
    
    # jsuffix <- paste0(jsuffixmain, ".", jmark, ".txt")
    print(jmark)
    
    if (jsuffixmain == "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table"){
      inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM.varfilt.dynamic_bins_genes.top_6085/lda_outputs.count_tables_merged.", jmark, ".", jsuffix, ".K-30.binarize.FALSE/ldaOut.count_tables_merged.", jmark, ".", jsuffix, ".K-30.Robj"))
    } else {
      inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM.varfilt.dynamic_bins_genes/lda_outputs.count_tables_merged.", jmark, ".", jsuffix, ".K-30.binarize.FALSE/ldaOut.count_tables_merged.", jmark, ".", jsuffix, ".K-30.Robj"))
    }
    assertthat::assert_that(file.exists(inf.lda))
    
    load(inf.lda, v=T)
    
    tm.result <- posterior(out.lda)
    
    dat.impute <- t(tm.result$topics %*% tm.result$terms)
    
    if (merge.ctypes.by.lineage){
      cnames.keep.lst <- split(x = dat.metas[[jmark]]$cell, f = dat.metas[[jmark]]$lineage)
    } else {
      cnames.keep.lst <- split(x = dat.metas[[jmark]]$cell, f = dat.metas[[jmark]]$cluster)
    }
    dat.impute.pbulk <- do.call(cbind, SumAcrossClusters(dat.impute, cnames.keep.lst))
    dat.impute.pbulk <- sweep(dat.impute.pbulk, MARGIN = 2, STATS = colSums(dat.impute.pbulk), FUN = "/", check.margin = TRUE)
    
    dat.raw.pbulk <- SumAcrossClusters(count.mat, cnames.keep.lst)
    dat.raw.pbulk <- do.call(cbind, dat.raw.pbulk)
    # normalize
    dat.raw.pbulk <- sweep(dat.raw.pbulk, MARGIN = 2, STATS = colSums(dat.raw.pbulk), FUN = "/", check.margin = TRUE) * scale.factor + 1
    return(list(dat.impute.pbulk = dat.impute.pbulk, dat.raw.pbulk = dat.raw.pbulk))
  })
  return(out.lst)
}



scale.factor <- 10^6

jmarks <- c("H3K9me3", "H3K4me1"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"



# Get dat meats -----------------------------------------------------------

dat.metas <- lapply(jmarks, function(jmark){
  inf.tmp <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  fread(inf.tmp)
})

cname2color <- hash::hash(dat.metas$H3K4me1$cluster, dat.metas$H3K4me1$clustercol)

# Load batch corrected mat  -----------------------------------------------

# indir.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM_dbl_reseq.varfilt.k4_k9_dynamic_bins"
indir.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins"
outs.lst <- lapply(jmarks, function(jmark){
  if (jmark == "H3K4me1"){
    dname <- "ldaAnalysisBins_mouse_spikein_BM_dbl_reseq.varfilt.k4_genes_only/lda_outputs.count_name.H3K4me1.k4_genes_only.2021-02-18.K-30.binarize.FALSE/ldaOut.count_name.H3K4me1.k4_genes_only.2021-02-18.K-30.Robj"
  } else {
    dname <- "ldaAnalysisBins_mouse_spikein_BM_dbl_reseq.varfilt.k9_genes_only/lda_outputs.count_name.H3K9me3.k9_genes_only.2021-02-18.K-30.binarize.FALSE/ldaOut.count_name.H3K9me3.k9_genes_only.2021-02-18.K-30.Robj"
  }
  inf.tmp <- file.path(indir.lda, dname)
  assertthat::assert_that(file.exists(inf.tmp))
  print(jmark)
  # inf.tmp <- file.path(indir.lda, paste0("lda_outputs.count_name.", jmark, ".k4_k9_dynamic_bins.2021-01-30.K-30.binarize.FALSE/ldaOut.count_name.", jmark, ".k4_k9_dynamic_bins.2021-01-30.K-30.Robj"))
  load(inf.tmp, v=T)
  return(list(out.lda = out.lda, count.mat = count.mat))
})

dat.raw.lst <- lapply(outs.lst, function(x){
  x$count.mat
})

out.lda.lst <- lapply(outs.lst, function(x){
  x$out.lda
})

tm.result.lst <- lapply(out.lda.lst, function(x){
  posterior(x)
})

dat.impute.lst <- lapply(tm.result.lst, function(tm.result){
  dat.impute <- t(tm.result$topics %*% tm.result$terms)
})


# Make pseudobulk  --------------------------------------------------------

dat.impute.pbulk.lst <- lapply(jmarks, function(jmark.tmp){
  print(jmark.tmp)
  cnames.keep.lst <- split(x = dat.metas[[jmark.tmp]]$cell, f = dat.metas[[jmark.tmp]]$cluster)
  sum.lst <- SumAcrossClusters(dat.impute.lst[[jmark.tmp]], cnames.keep.lst)
  pbulk <- do.call(cbind, sum.lst)
  # normalize
  pbulk <- sweep(pbulk, MARGIN = 2, STATS = colSums(pbulk), FUN = "/", check.margin = TRUE)
})

# compare with the other H3K9me3


jsuffixmain <- "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table"
out.lst <- GetPbulks(dat.metas, jmarks, jsuffixmain, merge.ctypes.by.lineage = FALSE, scale.factor = scale.factor)

rnames.common <- rownames(dat.impute.pbulk.lst$H3K9me3, out.lst$H3K9me3$dat.raw.pbulk)

print(dim(dat.impute.pbulk.lst$H3K9me3))
print(dim(out.lst$H3K9me3$dat.impute.pbulk))

print(head(dat.impute.pbulk.lst$H3K9me3[rnames.common, ]))
print(head(out.lst$H3K9me3$dat.impute.pbulk[rnames.common, ]))

print(head(dat.raw.pbulk.lst$H3K9me3[rnames.common, ]))
print(head(out.lst$H3K9me3$dat.raw.pbulk[rnames.common, ]))

plot(log2(dat.impute.pbulk.lst$H3K9me3[rnames.common, 1]), log2(out.lst$H3K9me3$dat.impute.pbulk[rnames.common, 1]))
# plot(log2(dat.raw.pbulk.lst$H3K9me3[rnames.common, 1]), log2(out.lst$H3K9me3$dat.raw.pbulk[, 1]))


dat.raw.pbulk.lst <- lapply(jmarks, function(jmark.tmp){
  print(jmark.tmp)
  cnames.keep.lst <- split(x = dat.metas[[jmark.tmp]]$cell, f = dat.metas[[jmark.tmp]]$cluster)
  sum.lst <- SumAcrossClusters(dat.raw.lst[[jmark.tmp]], cnames.keep.lst)
  pbulk <- do.call(cbind, sum.lst)
  # normalize
  pbulk <- sweep(pbulk, MARGIN = 2, STATS = colSums(pbulk), FUN = "/", check.margin = TRUE)
  pbulk <- pbulk * scale.factor + 1
})

# jmark.tmp <- "H3K9me3"
jmark.tmp <- "H3K4me1"

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_differential_expression_outputs/heatmaps_scchix"
# use.imputed <- FALSE
tf.vec <- c(TRUE, FALSE)

colpalette <- colorRampPalette(c("grey1", "grey35", "grey99"))(1024)


for (use.imputed in tf.vec){
  print(use.imputed)
  for (jmark.tmp in jmarks){
    print(jmark.tmp)
    fname <- paste0("heatmaps_LDA_pseudobulk.mark_specific_regions_only.ward_hclust.", jmark.tmp, ".use_imputed.", use.imputed, ".", Sys.Date(), ".pdf")
    pdf(file.path(outdir, fname), useDingbats = FALSE)
    
    if (use.imputed){
      dat.tmp <- dat.impute.pbulk.lst[[jmark.tmp]]
      jmat2.check <- out.lst[[jmark.tmp]]$dat.impute.pbulk
    } else {
      dat.tmp <- dat.raw.pbulk.lst[[jmark.tmp]]
      jmat2.check <- out.lst[[jmark.tmp]]$dat.impute.pbulk
    }
    # jmat2 <- t(apply(dat.tmp, 1, function(jrow) Winsorize(jrow, probs = c(0.01, 0.99))))
    # jmat2 <- apply(jmat2, 2, function(jcol) Winsorize(jcol, probs = c(0.01, 0.99)))
    jmat2 <- dat.tmp
    
    # rnames <- rownames(jmat2)
    # rnames.label <- sapply(rnames, function(x) ifelse(grepl(";chr", x), "K9", "K4"))
    # rnames.color <- sapply(rnames.label, function(x) ifelse(x == "K4", "red", "blue"))
    
    cnames.color <- sapply(colnames(jmat2), function(x) AssignHash(x, cname2color, null.fill = x))
    
    # heatmap3::heatmap3(jmat.lst[[jmark]], Rowv = NA, Colv = NA, scale = "row", ColSideColors = jmetas.pretty.lst[[jmark]]$clustercol, RowSideColors = sapply(bins.keep.common, AssignHash, jhash = bin2col),  revC = TRUE, main = paste0(jmark, " 50kb bins"), margins = c(5, 8), col = colpalette)
    
    # heatmap3::heatmap3(log2(jmat2), Rowv = TRUE, Colv = TRUE, scale = "column",  revC = TRUE, main = paste0(jmark, " 50kb bins"), margins = c(5, 8))
    # heatmap3::heatmap3(log2(jmat2), Rowv = TRUE, Colv = TRUE, scale = "column",  revC = TRUE, main = paste0(jmark.tmp, " 50kb bins and TES. Complete"), margins = c(5, 8), RowSideColors = rnames.color, ColSideColors = cnames.color, col = colpalette, cexRow = 0.5, method = "complete")
    heatmap3::heatmap3(log2(jmat2), Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, main = paste0(jmark.tmp, " DynamicOnly. Complete"), margins = c(5, 8), ColSideColors = cnames.color, col = colpalette, cexRow = 0.5, method = "complete")
    heatmap3::heatmap3(log2(jmat2.check), Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, main = paste0("Check ", jmark.tmp, "50kb bins and TES. Complete"), margins = c(5, 8), ColSideColors = cnames.color, col = colpalette, cexRow = 0.5, method = "complete")
    
    # heatmap3::heatmap3(log2(jmat2), Rowv = TRUE, Colv = TRUE, scale = "column",  revC = TRUE, main = paste0(jmark.tmp, " 50kb bins and TES. Ward"), margins = c(5, 8), RowSideColors = rnames.color, ColSideColors = cnames.color, col = colpalette, cexRow = 0.5, method = "ward.D2")
    heatmap3::heatmap3(log2(jmat2), Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, main = paste0(jmark.tmp, " DynamicOnly. Ward"), margins = c(5, 8), ColSideColors = cnames.color, col = colpalette, cexRow = 0.5, method = "ward.D2")
    heatmap3::heatmap3(log2(jmat2.check), Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, main = paste0("Check ", jmark.tmp, " 50kb bins and TES. Ward"), margins = c(5, 8), ColSideColors = cnames.color, col = colpalette, cexRow = 0.5, method = "ward.D2")
    dev.off()
  }
}

# 
# 
# # Check most weights are in K4 or K9 bins respectively  -------------------
# 
# # jmark.tmp <- "H3K4me1"
# jmark.tmp <- "H3K9me3"


