# Jake Yeung
# Date of Creation: 2021-04-07
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/dendrograms/5-downstream_make_plots_of_dynamic_bins.from_raw_CPM.check_nozscore_no_log2FC.order_rows.DE_bins_same_as_4B.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)
library(hash)
library(topicmodels)

# Funmctions --------------------------------------------------------------


GetPbulks <- function(dat.metas, jmarks, scale.factor = 1000, merge.ctypes.by.lineage = FALSE){
  out.lst <- lapply(jmarks, function(jmark){
    # jsuffixmain <- "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table.H3K9me3.2021-02-15.txt"
    inf.lda <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM.varfilt.dynamic_bins_genes.corrected_DE_tables/lda_outputs.count_tables_merged.", jmark ,".dynamic_bins.50kb.corrected_DE_tables.", jmark, ".2021-04-07.txt.corrected_DE_tables.K-30.binarize.FALSE/ldaOut.count_tables_merged.", jmark, ".dynamic_bins.50kb.corrected_DE_tables.", jmark, ".2021-04-07.txt.corrected_DE_tables.K-30.Robj")
    load(inf.lda, v=T)
    
    tm.result <- posterior(out.lda)
    
    dat.impute <- t(tm.result$topics %*% tm.result$terms)
    rownames(dat.impute) <- sapply(rownames(dat.impute), function(x) strsplit(x, ";")[[1]][[2]])
    rownames(count.mat) <- sapply(rownames(count.mat), function(x) strsplit(x, ";")[[1]][[2]])
    
    # filter bins 
    inf.tmp <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/DE_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.jmidbug_fixed.", jmark, ".2021-02-19.txt")
    dat.de.bins.tmp <- fread(inf.tmp)
    
    dat.impute.filt <- dat.impute[dat.de.bins.tmp$CoordOriginal, ]
    count.mat.filt <- count.mat[dat.de.bins.tmp$CoordOriginal, ]
    
    
    if (merge.ctypes.by.lineage){
      cnames.keep.lst <- split(x = dat.metas[[jmark]]$cell, f = dat.metas[[jmark]]$lineage)
    } else {
      cnames.keep.lst <- split(x = dat.metas[[jmark]]$cell, f = dat.metas[[jmark]]$cluster)
    }
    dat.impute.pbulk <- do.call(cbind, SumAcrossClusters(dat.impute.filt, cnames.keep.lst))
    dat.impute.pbulk <- sweep(dat.impute.pbulk, MARGIN = 2, STATS = colSums(dat.impute.pbulk), FUN = "/", check.margin = TRUE)
    
    dat.raw.pbulk <- SumAcrossClusters(count.mat.filt, cnames.keep.lst)
    dat.raw.pbulk <- do.call(cbind, dat.raw.pbulk)
    # normalize
    print(scale.factor)
    dat.raw.pbulk <- sweep(dat.raw.pbulk, MARGIN = 2, STATS = colSums(dat.raw.pbulk), FUN = "/", check.margin = TRUE) * scale.factor + 1
    return(list(dat.impute.pbulk = dat.impute.pbulk, dat.raw.pbulk = dat.raw.pbulk))
  })
  return(out.lst)
}

make.plots <- TRUE

# Constants ---------------------------------------------------------------

ctypes <- list("Eryths" = "Erythroid", 
               "Bcells" = "Lymphoid", 
               "NKs" = "Lymphoid", 
               "Granulocytes" = "Myeloid",
               "Basophils" = "Myeloid", 
               "pDCs" = "Lymphoid",
               "DCs" = "Myeloid", 
               "HSPCs" = "HSPCs",
               "Erythroid" = "Erythroid",
               "Lymphoid" = "Lymphoid",
               "Myeloid" = "Myeloid")



binsize <- 50000
jscale <- 1000000L
jymax <- 20000

# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"


# Load pbulks from 50kb ---------------------------------------------------


# Load metadata -----------------------------------------------------------

indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage")

dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  fread(file.path(indir.meta, fname)) %>%
    rowwise() %>%
    mutate(lineage = ctypes[[cluster]])
}) 

# cluster to col
jmarkref <- "H3K4me3"
# cluster2col <- hash(dat.metas$H3K4me1$cluster, dat.metas$H3K4me1$clustercol)
cluster2col <- hash(dat.metas[[jmarkref]]$cluster, dat.metas[[jmarkref]]$clustercol)
cluster2col[["Erythroid"]] <- "#0072B2"
cluster2col[["Lymphoid"]] <- "#56B4E9"
cluster2col[["Myeloid"]] <- "#D55E00"


cname2color <- hash::hash(dat.metas[[jmarkref]]$cluster, dat.metas[[jmarkref]]$clustercol)
colpalette <- colorRampPalette(c("grey1", "grey35", "grey99"))(1024)


# Load pbulks  ------------------------------------------------------------

# jsuffixmain <- "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table"
# jsuffixmain2 <- "celltype_specific_genes.TSS_10kb.txt"
# jsuffixmain3 <- "celltype_specific_genes.TSS_TES.txt"

out.lst <- GetPbulks(dat.metas, jmarks, merge.ctypes.by.lineage = FALSE, scale.factor = jscale)
# out.lst.tss <- GetPbulks(dat.metas, jmarks, jsuffixmain2, merge.ctypes.by.lineage = FALSE)
# out.lst.tes <- GetPbulks(dat.metas, jmarks[1:3], jsuffixmain3, merge.ctypes.by.lineage = FALSE)

# jmark <- "H3K4me1"
# jmat2 <- log2(out.lst.tes[[jmark]]$dat.impute.pbulk)
# heatmap3::heatmap3(jmat2, 
#                    Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, 
#                    main = paste("TSS", jmark, jmethod), margins = c(5, 8), cexRow = 0.5, method = jmethod, ColSideColors = cnames.color, col = colpalette)

# # check K4me1
# jinf.check <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps/table.H3K4me1.dynamic_bins.50kb.ward.D2.1000000.2021-02-19.txt")
# jcheck <- read.table(jinf.check, header = TRUE, row.names = 1)
# jref <- out.lst$H3K4me1$dat.raw.pbulk
# jref2 <- GetPbulks(dat.metas, jmarks[["H3K4me1"]], jsuffixmain, merge.ctypes.by.lineage = FALSE, scale.factor = 1000000)[[1]]$dat.raw.pbulk
# jref3 <- GetPbulks(dat.metas, jmarks[["H3K4me1"]], jsuffixmain, merge.ctypes.by.lineage = FALSE, scale.factor = 1000)[[1]]$dat.raw.pbulk
# 
# dim(jcheck)
# dim(jref)
# 
# rnames.check <- rownames(jcheck)
# rnames.ref <- rownames(jref)
# rnames.common <- intersect(rnames.check, rnames.ref)


# Make heatmaps  ----------------------------------------------------------

#     heatmap3::heatmap3(log2(jmat2), Rowv = TRUE, Colv = TRUE, scale = "column",  revC = TRUE, main = paste0(jmark.tmp, " 50kb bins and TES. Complete"), margins = c(5, 8), RowSideColors = rnames.color, ColSideColors = cnames.color, col = colpalette, cexRow = 0.5, method = "complete")

jmethod <- "ward.D2"
# outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps/heatmaps_final.impute.", jmethod, ".pdf")
# pdf(outpdf, useDingbats = FALSE)
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps"

jmethod <- "ward.D2"
if (make.plots){
  pdf(file.path(outdir, paste0("primetime_dynamic_features_and_heatmaps.from_CPM.", jmethod, ".", Sys.Date(), ".noscaling_nolog2.order_rows.corrected_DE_tables.pdf")), useDingbats = FALSE)
}

# jmethod <- "complete"
# jmethod <- "average"
# jmethod <- "single"
# jmethod <- "median"
# jmethod <- "centroid"
# jmethod <- "mcquitty"
for (jmark in jmarks){
  # jmat2 <- log2(out.lst[[jmark]]$dat.impute.pbulk)
  jmat2 <- log2(out.lst[[jmark]]$dat.raw.pbulk)
  # jmat2 <- sweep(jmat2, MARGIN = 1, STATS = rowMeans(jmat2), FUN = "-")
  cnames.color <- sapply(colnames(jmat2), function(x) AssignHash(x, cname2color, null.fill = x))
  heatmap3::heatmap3(jmat2, 
                     Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, 
                     main = paste("50kb", jmark, jmethod), margins = c(5, 8), cexRow = 0.5, method = jmethod, ColSideColors = cnames.color, col = colpalette)
}


for (jmark in jmarks){
  # jmat2 <- log2(out.lst[[jmark]]$dat.impute.pbulk)
  jmat2 <- log2(out.lst[[jmark]]$dat.raw.pbulk)
  # jmat2 <- sweep(jmat2, MARGIN = 1, STATS = rowMeans(jmat2), FUN = "-")
  cnames.color <- sapply(colnames(jmat2), function(x) AssignHash(x, cname2color, null.fill = x))
  heatmap3::heatmap3(jmat2, 
                     Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, 
                     main = paste("50kb", jmark, jmethod), margins = c(5, 8), cexRow = 0.5, method = jmethod, ColSideColors = cnames.color)
}




# dev.off()

# order rows by decreasing HSPC z-score
for (jmark in jmarks){
  # jmat2 <- log2(out.lst[[jmark]]$dat.impute.pbulk)
  jmat2 <- log2(out.lst[[jmark]]$dat.raw.pbulk)
  jmat2.zscore <- t(scale(t(jmat2), center = TRUE, scale = TRUE))
  # order rows manually 
  rows.order <- order(jmat2.zscore[, "HSPCs"], decreasing = TRUE)
  jmat2 <- jmat2[rows.order, ]
  
  cnames.color <- sapply(colnames(jmat2), function(x) AssignHash(x, cname2color, null.fill = x))
  heatmap3::heatmap3(jmat2, 
                     Rowv = NA, Colv = TRUE, scale = "row",  revC = TRUE, 
                     main = paste("50kb", jmark, jmethod), margins = c(5, 8), cexRow = 0.5, method = jmethod, ColSideColors = cnames.color, col = colpalette)
}

# 
# for (jmark in jmarks){
#   # jmat2 <- log2(out.lst[[jmark]]$dat.impute.pbulk)
#   jmat2 <- log2(out.lst[[jmark]]$dat.raw.pbulk)
#   jmat2 <- sweep(jmat2, MARGIN = 1, STATS = rowMeans(jmat2), FUN = "-")
#   cnames.color <- sapply(colnames(jmat2), function(x) AssignHash(x, cname2color, null.fill = x))
#   heatmap3::heatmap3(jmat2, 
#                      Rowv = TRUE, Colv = TRUE, scale = "none",  revC = TRUE, 
#                      main = paste("50kb", jmark, jmethod, "log2 center no scale"), margins = c(5, 8), cexRow = 0.5, method = jmethod, ColSideColors = cnames.color, col = colpalette)
# }
# 
# for (jmark in jmarks){
#   # jmat2 <- log2(out.lst[[jmark]]$dat.impute.pbulk)
#   jmat2 <- out.lst[[jmark]]$dat.raw.pbulk
#   cnames.color <- sapply(colnames(jmat2), function(x) AssignHash(x, cname2color, null.fill = x))
#   heatmap3::heatmap3(jmat2, 
#                      Rowv = TRUE, Colv = TRUE, scale = "none",  revC = TRUE, 
#                      main = paste("50kb", jmark, jmethod, "nolog2 noscale"), margins = c(5, 8), cexRow = 0.5, method = jmethod, ColSideColors = cnames.color, col = colpalette)
# }



# Plot distances to gene ----------------------------------------------------------

infs.highbins <- lapply(jmarks, function(jmark){
  file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/High_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.", jmark, ".2021-01-30.txt"))
})

infs.debins <- lapply(jmarks, function(jmark){
  file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/tables_top_6085_four_marks_dynamic_bins/top_6085_bins_nearest_gene_gc.", jmark, ".2021-02-17.txt"))
})

dat.highbins.lst <- lapply(infs.highbins, function(inf){
  fread(inf)
})

dat.debins.lst <- lapply(infs.debins, function(inf){
  fread(inf)
})


# 
# 
# dat.bins.merge.long <- lapply(jmarks, function(jmark){
#   rbind(subset(dat.highbins.lst[[jmark]], select = c(distanceToTSS, region_coord, mark)) %>% mutate(type = "High"), 
#         subset(dat.debins.lst[[jmark]], select = c(distanceToTSS, region_coord, mark)) %>% mutate(type = "DE"))
# }) %>%
#   bind_rows()
# 
# ggplot(dat.bins.merge.long, aes(x = log10(distanceToTSS + 1), fill = type)) + 
#   geom_density(alpha = 0.25) + 
#   facet_wrap(~mark, nrow = 1) + 
#   geom_vline(xintercept = log10(25000), linetype = "dotted") + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
# 

# Plot GC  ----------------------------------------------------------------


ggplot(dat.debins.lst %>% bind_rows(), aes(x = mark, y = gc)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point() + 
  xlab("") + 
  ylab("GC content") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


if (make.plots){
  dev.off()
}


