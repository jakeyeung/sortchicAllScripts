# Jake Yeung
# Date of Creation: 2021-03-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/dendrograms/5-downstream_make_plots_of_dynamic_bins.from_raw_CPM.check_nozscore_no_log2FC.R
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


GetPbulks <- function(dat.metas, jmarks, jsuffixmain, scale.factor = 1000, merge.ctypes.by.lineage = FALSE){
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
jscale <- 10^6
jymax <- 20000

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
cluster2col <- hash(dat.metas$H3K4me1$cluster, dat.metas$H3K4me1$clustercol)
cluster2col[["Erythroid"]] <- "#0072B2"
cluster2col[["Lymphoid"]] <- "#56B4E9"
cluster2col[["Myeloid"]] <- "#D55E00"


cname2color <- hash::hash(dat.metas$H3K4me1$cluster, dat.metas$H3K4me1$clustercol)
colpalette <- colorRampPalette(c("grey1", "grey35", "grey99"))(1024)


# Load pbulks  ------------------------------------------------------------

jsuffixmain <- "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table"
# jsuffixmain2 <- "celltype_specific_genes.TSS_10kb.txt"
# jsuffixmain3 <- "celltype_specific_genes.TSS_TES.txt"

out.lst <- GetPbulks(dat.metas, jmarks, jsuffixmain, merge.ctypes.by.lineage = FALSE, scale.factor = jscale)
# out.lst.tss <- GetPbulks(dat.metas, jmarks, jsuffixmain2, merge.ctypes.by.lineage = FALSE)
# out.lst.tes <- GetPbulks(dat.metas, jmarks[1:3], jsuffixmain3, merge.ctypes.by.lineage = FALSE)

# jmark <- "H3K4me1"
# jmat2 <- log2(out.lst.tes[[jmark]]$dat.impute.pbulk)
# heatmap3::heatmap3(jmat2, 
#                    Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, 
#                    main = paste("TSS", jmark, jmethod), margins = c(5, 8), cexRow = 0.5, method = jmethod, ColSideColors = cnames.color, col = colpalette)

# check K4me1
jinf.check <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps/table.H3K4me1.dynamic_bins.50kb.ward.D2.1000000.2021-02-19.txt")
jcheck <- read.table(jinf.check, header = TRUE, row.names = 1)
jref <- out.lst$H3K4me1$dat.raw.pbulk
jref2 <- GetPbulks(dat.metas, jmarks[["H3K4me1"]], jsuffixmain, merge.ctypes.by.lineage = FALSE, scale.factor = 1000000)[[1]]$dat.raw.pbulk
jref3 <- GetPbulks(dat.metas, jmarks[["H3K4me1"]], jsuffixmain, merge.ctypes.by.lineage = FALSE, scale.factor = 1000)[[1]]$dat.raw.pbulk

dim(jcheck)
dim(jref)

rnames.check <- rownames(jcheck)
rnames.ref <- rownames(jref)
rnames.common <- intersect(rnames.check, rnames.ref)


# Make heatmaps  ----------------------------------------------------------

#     heatmap3::heatmap3(log2(jmat2), Rowv = TRUE, Colv = TRUE, scale = "column",  revC = TRUE, main = paste0(jmark.tmp, " 50kb bins and TES. Complete"), margins = c(5, 8), RowSideColors = rnames.color, ColSideColors = cnames.color, col = colpalette, cexRow = 0.5, method = "complete")

jmethod <- "ward.D2"
# outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps/heatmaps_final.impute.", jmethod, ".pdf")
# pdf(outpdf, useDingbats = FALSE)
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps"

if (make.plots){
  pdf(file.path(outdir, paste0("primetime_top_6085_features_and_heatmaps.from_CPM.", Sys.Date(), ".noscaling_nolog2.pdf")), useDingbats = FALSE)
}

jmethod <- "complete"
for (jmark in jmarks){
  # jmat2 <- log2(out.lst[[jmark]]$dat.impute.pbulk)
  jmat2 <- log2(out.lst[[jmark]]$dat.raw.pbulk)
  cnames.color <- sapply(colnames(jmat2), function(x) AssignHash(x, cname2color, null.fill = x))
  heatmap3::heatmap3(jmat2, 
                     Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, 
                     main = paste("50kb", jmark, jmethod), margins = c(5, 8), cexRow = 0.5, method = jmethod, ColSideColors = cnames.color, col = colpalette)
}
# dev.off()

for (jmark in jmarks){
  # jmat2 <- log2(out.lst[[jmark]]$dat.impute.pbulk)
  jmat2 <- log2(out.lst[[jmark]]$dat.raw.pbulk)
  cnames.color <- sapply(colnames(jmat2), function(x) AssignHash(x, cname2color, null.fill = x))
  heatmap3::heatmap3(jmat2, 
                     Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, 
                     main = paste("50kb", jmark, jmethod), margins = c(5, 8), cexRow = 0.5, method = jmethod, ColSideColors = cnames.color, col = colpalette)
}


for (jmark in jmarks){
  # jmat2 <- log2(out.lst[[jmark]]$dat.impute.pbulk)
  jmat2 <- log2(out.lst[[jmark]]$dat.raw.pbulk)
  jmat2 <- sweep(jmat2, MARGIN = 1, STATS = rowMeans(jmat2), FUN = "-")
  cnames.color <- sapply(colnames(jmat2), function(x) AssignHash(x, cname2color, null.fill = x))
  heatmap3::heatmap3(jmat2, 
                     Rowv = TRUE, Colv = TRUE, scale = "none",  revC = TRUE, 
                     main = paste("50kb", jmark, jmethod, "log2 center no scale"), margins = c(5, 8), cexRow = 0.5, method = jmethod, ColSideColors = cnames.color, col = colpalette)
}

for (jmark in jmarks){
  # jmat2 <- log2(out.lst[[jmark]]$dat.impute.pbulk)
  jmat2 <- out.lst[[jmark]]$dat.raw.pbulk
  cnames.color <- sapply(colnames(jmat2), function(x) AssignHash(x, cname2color, null.fill = x))
  heatmap3::heatmap3(jmat2, 
                     Rowv = TRUE, Colv = TRUE, scale = "none",  revC = TRUE, 
                     main = paste("50kb", jmark, jmethod, "nolog2 noscale"), margins = c(5, 8), cexRow = 0.5, method = jmethod, ColSideColors = cnames.color, col = colpalette)
}




# Calculate variance: both pos and neg ------------------------------------

dat.vars.lst <- lapply(out.lst, function(jout){
  # dat <- jout$dat.impute.pbulk
  jdat <- jout$dat.raw.pbulk
  jdat.log2 <- log2(jdat)
  jcheck <- sweep(jdat.log2, MARGIN = 1, STATS = rowMeans(jdat.log2), FUN = "-") 
  jcheck <- jcheck ^ 2
})

dat.diffs.lst <- lapply(out.lst, function(jout){
  # jdat <- jout$dat.impute.pbulk
  jdat <- jout$dat.raw.pbulk
  jdat.log2 <- log2(jdat)
  jcheck <- sweep(jdat.log2, MARGIN = 1, STATS = rowMeans(jdat.log2), FUN = "-") 
  # jcheck <- jcheck ^ 2
})

# calculate variance, break down into pos or neg
dat.vars.posneg.long.lst <- lapply(jmarks, function(jmark){
  dat.diffs <- dat.diffs.lst[[jmark]]
  dat.diffs.long <- melt(dat.diffs)
  colnames(dat.diffs.long) <- c("rname", "celltype", "log2diff")
  dat.diffs.sum.long <- dat.diffs.long %>%
    rowwise() %>%
    mutate(is.pos = log2diff > 0) %>%
    group_by(celltype, is.pos) %>%
    summarise(vartotal = sum(log2diff ^ 2)) %>%
    rowwise() %>%
    mutate(lineage = ctypes[[celltype]],
           clustercol = cluster2col[[as.character(celltype)]])
  dat.diffs.sum.long$mark <- jmark
  return(dat.diffs.sum.long)
})

# check they're equal 
dat.vars.long.lst <- lapply(jmarks, function(jmark){
  dat.vars <- dat.vars.lst[[jmark]]
  melt(dat.vars) %>%
    dplyr::rename(bin = Var1, 
                  celltype = Var2) %>%
    group_by(celltype) %>%
    summarise(vartotal = sum(value)) %>%
    mutate(mark = jmark) %>%
    rowwise() %>%
    mutate(lineage = ctypes[[celltype]],
           clustercol = cluster2col[[as.character(celltype)]])
})
# dat.vars.long <- dat.vars.long.lst %>%
#   bind_rows()


# Plot variance, then split by posneg -------------------------------------


for (jmark in jmarks){
  m <- ggplot(dat.vars.long.lst[[jmark]], aes(x = forcats::fct_reorder(.f = celltype, .x = vartotal, .fun = median, .desc = TRUE), y = vartotal, fill = clustercol)) + 
    geom_col() + 
    ylim(c(0, jymax)) + 
    scale_fill_identity() + 
    xlab("") + 
    ylab("Total Variance") + 
    facet_wrap(~mark, scales = "free_x") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
}

for (jmark in jmarks){
  m <- ggplot(dat.vars.posneg.long.lst[[jmark]], aes(x = forcats::fct_reorder(.f = celltype, .x = vartotal, .fun = median, .desc = TRUE), y = vartotal, fill = is.pos)) + 
    # geom_col(position = "stack") + 
    geom_col() + 
    ylim(c(0, jymax)) + 
    facet_wrap(~mark, scales = "free_x") + 
    theme_bw() + 
    xlab("") + 
    ylab("Total Variance") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")
  print(m)
}

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

dat.bins.merge.long <- lapply(jmarks, function(jmark){
  rbind(subset(dat.highbins.lst[[jmark]], select = c(distanceToTSS, region_coord, mark)) %>% mutate(type = "High"), 
        subset(dat.debins.lst[[jmark]], select = c(distanceToTSS, region_coord, mark)) %>% mutate(type = "DE"))
}) %>%
  bind_rows()

ggplot(dat.bins.merge.long, aes(x = log10(distanceToTSS + 1), fill = type)) + 
  geom_density(alpha = 0.25) + 
  facet_wrap(~mark, nrow = 1) + 
  geom_vline(xintercept = log10(25000), linetype = "dotted") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


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


