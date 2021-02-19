# Jake Yeung
# Date of Creation: 2021-02-06
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/19-make_heatmap_K4me1_and_K9me3_dynamic_bins.from_LDA.R
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

jmarks <- c("H3K9me3", "H3K4me1"); names(jmarks) <- jmarks


# Get dat meats -----------------------------------------------------------

dat.metas <- lapply(jmarks, function(jmark){
  inf.tmp <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  fread(inf.tmp)
})

cname2color <- hash::hash(dat.metas$H3K4me1$cluster, dat.metas$H3K4me1$clustercol)

# Load batch corrected mat  -----------------------------------------------

indir.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM_dbl_reseq.varfilt.k4_k9_dynamic_bins"
outs.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.tmp <- file.path(indir.lda, paste0("lda_outputs.count_name.", jmark, ".k4_k9_dynamic_bins.2021-01-30.K-30.binarize.FALSE/ldaOut.count_name.", jmark, ".k4_k9_dynamic_bins.2021-01-30.K-30.Robj"))
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


dat.raw.pbulk.lst <- lapply(jmarks, function(jmark.tmp){
  print(jmark.tmp)
  cnames.keep.lst <- split(x = dat.metas[[jmark.tmp]]$cell, f = dat.metas[[jmark.tmp]]$cluster)
  sum.lst <- SumAcrossClusters(dat.raw.lst[[jmark.tmp]], cnames.keep.lst)
  pbulk <- do.call(cbind, sum.lst)
  # normalize
  pbulk <- sweep(pbulk, MARGIN = 2, STATS = colSums(pbulk), FUN = "/", check.margin = TRUE)
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
    fname <- paste0("heatmaps_LDA_pseudobulk.ward_hclust.", jmark.tmp, ".use_imputed.", use.imputed, ".", Sys.Date(), ".pdf")
    pdf(file.path(outdir, fname), useDingbats = FALSE)
    
    if (use.imputed){
      dat.tmp <- dat.impute.pbulk.lst[[jmark.tmp]]
    } else {
      dat.tmp <- dat.raw.pbulk.lst[[jmark.tmp]]
    }
    jmat2 <- t(apply(dat.tmp, 1, function(jrow) Winsorize(jrow, probs = c(0.01, 0.99))))
    jmat2 <- apply(jmat2, 2, function(jcol) Winsorize(jcol, probs = c(0.01, 0.99)))
    
    rnames <- rownames(jmat2)
    rnames.label <- sapply(rnames, function(x) ifelse(grepl(";chr", x), "K9", "K4"))
    rnames.color <- sapply(rnames.label, function(x) ifelse(x == "K4", "red", "blue"))
    
    cnames.color <- sapply(colnames(jmat2), function(x) AssignHash(x, cname2color, null.fill = x))
    
    # heatmap3::heatmap3(jmat.lst[[jmark]], Rowv = NA, Colv = NA, scale = "row", ColSideColors = jmetas.pretty.lst[[jmark]]$clustercol, RowSideColors = sapply(bins.keep.common, AssignHash, jhash = bin2col),  revC = TRUE, main = paste0(jmark, " 50kb bins"), margins = c(5, 8), col = colpalette)
    
    # heatmap3::heatmap3(log2(jmat2), Rowv = TRUE, Colv = TRUE, scale = "column",  revC = TRUE, main = paste0(jmark, " 50kb bins"), margins = c(5, 8))
    heatmap3::heatmap3(log2(jmat2), Rowv = TRUE, Colv = TRUE, scale = "column",  revC = TRUE, main = paste0(jmark.tmp, " 50kb bins and TES. Complete"), margins = c(5, 8), RowSideColors = rnames.color, ColSideColors = cnames.color, col = colpalette, cexRow = 0.5, method = "complete")
    heatmap3::heatmap3(log2(jmat2), Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, main = paste0(jmark.tmp, " 50kb bins and TES. Complete"), margins = c(5, 8), RowSideColors = rnames.color, ColSideColors = cnames.color, col = colpalette, cexRow = 0.5, method = "complete")
    
    heatmap3::heatmap3(log2(jmat2), Rowv = TRUE, Colv = TRUE, scale = "column",  revC = TRUE, main = paste0(jmark.tmp, " 50kb bins and TES. Ward"), margins = c(5, 8), RowSideColors = rnames.color, ColSideColors = cnames.color, col = colpalette, cexRow = 0.5, method = "ward.D2")
    heatmap3::heatmap3(log2(jmat2), Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, main = paste0(jmark.tmp, " 50kb bins and TES. Ward"), margins = c(5, 8), RowSideColors = rnames.color, ColSideColors = cnames.color, col = colpalette, cexRow = 0.5, method = "ward.D2")
    
    # filter only K4 or K9 bins 
    
    jmat2.k4only <- jmat2[rnames.label == "K4", ]
    jmat2.k9only <- jmat2[rnames.label == "K9", ]
    heatmap3::heatmap3(log2(jmat2.k4only), Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, main = paste0(jmark.tmp, " K4 dynamic only. Complete"), margins = c(5, 8), ColSideColors = cnames.color, col = colpalette, cexRow = 0.5, method = "complete")
    heatmap3::heatmap3(log2(jmat2.k9only), Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, main = paste0(jmark.tmp, " K9 dynamic only. Complete"), margins = c(5, 8), ColSideColors = cnames.color, col = colpalette, cexRow = 0.5, method = "complete")
    
    heatmap3::heatmap3(log2(jmat2.k4only), Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, main = paste0(jmark.tmp, " K4 dynamic only. Ward"), margins = c(5, 8), ColSideColors = cnames.color, col = colpalette, cexRow = 0.5, method = "ward.D2")
    heatmap3::heatmap3(log2(jmat2.k9only), Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, main = paste0(jmark.tmp, " K9 dynamic only. Ward"), margins = c(5, 8), ColSideColors = cnames.color, col = colpalette, cexRow = 0.5, method = "ward.D2")
    
    
    rnames <- rownames(dat.tmp)
    rnames.label <- sapply(rnames, function(x) ifelse(grepl(";chr", x), "K9", "K4"))
    rnames.color <- sapply(rnames.label, function(x) ifelse(x == "K4", "red", "blue"))
    
    dat.tmp.mean <- data.frame(bin = rownames(dat.tmp), mean.exprs = rowMeans(dat.tmp), stringsAsFactors = FALSE) %>%
      rowwise() %>%
      mutate(label = ifelse(grepl(";chr", bin), "K9", "K4"))
    
    m <- ggplot(dat.tmp.mean, aes(x = log2(mean.exprs * 10^5 + 1), fill = label)) + 
      geom_density(alpha = 0.25) + 
      theme_bw() + 
      ggtitle(jmark.tmp) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
    
    # write to file 
    fname.k4only <- file.path(outdir, paste0("K4_regions_matrix_input.use_imputed.", use.imputed, ".mark.", jmark.tmp, ".log2scale.txt"))
    fname.k9only <- file.path(outdir, paste0("K9_regions_matrix_input.use_imputed.", use.imputed, ".mark.", jmark.tmp, ".log2scale.txt"))
    
    write.table(log2(jmat2.k4only), file = fname.k4only, row.names = TRUE, col.names = NA, quote = FALSE)
    write.table(log2(jmat2.k9only), file = fname.k9only, row.names = TRUE, col.names = NA, quote = FALSE)
    
    # fwrite(jmat2.k4only, file = fname.k4only, row.names = TRUE, col.names = TRUE)
    # fwrite(jmat2.k9only, file = fname.k9only, row.names = TRUE, col.names = TRUE)
    
    dev.off()
  }
}

# 
# 
# # Check most weights are in K4 or K9 bins respectively  -------------------
# 
# # jmark.tmp <- "H3K4me1"
# jmark.tmp <- "H3K9me3"


