# Jake Yeung
# Date of Creation: 2021-02-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/dendrograms/3-calculate_variance_on_pseudobulk.R
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

library(hash)

merge.ctypes.by.lineage <- FALSE
# merge.ctypes.by.lineage <- TRUE

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps"
outpdf <- file.path(outdir, paste0("variance_by_clusters_top_6085_bins.merge_ctypes.", merge.ctypes.by.lineage, ".", Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)

# merge some celltypes
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

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"


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


# Load LDA  ---------------------------------------------------------------


jmark <- "H3K4me3"
# jsuffixmain <- paste0("dynamic_bins.50kb")
jsuffixmain <- "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table"
# jsuffix <- "celltype_specific_genes.TSS_10kb.txt"
# jsuffix <- "celltype_specific_genes.TSS_TES.txt"
# dat.impute.pbulk.lst <- lapply(jmarks, function(jmark){

jsuffix.lst <- list("dynamic_bins.50kb", "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table")
names(jsuffix.lst) <- jsuffix.lst

dat.raw.pbulk.lst.lst <- lapply(jsuffix.lst, function(jsuffixmain){
  out.lst <- lapply(jmarks, function(jmark){
    # jsuffixmain <- "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table.H3K9me3.2021-02-15.txt"
    
    if (jsuffixmain == "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table"){
      jsuffix <- paste0(jsuffixmain, ".", jmark, ".2021-02-15.txt")
    } else {
      jsuffix <- paste0("dynamic_bins.50kb.", jmark, ".txt")
    }
    
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
    dat.raw.pbulk <- sweep(dat.raw.pbulk, MARGIN = 2, STATS = colSums(dat.raw.pbulk), FUN = "/", check.margin = TRUE) * 1000 + 1
    return(list(dat.impute.pbulk = dat.impute.pbulk, dat.raw.pbulk = dat.raw.pbulk))
  })
  
  dat.impute.pbulk.lst <- lapply(out.lst, function(x) x$dat.impute.pbulk)
  dat.raw.pbulk.lst <- lapply(out.lst, function(x) x$dat.raw.pbulk)
  return(dat.raw.pbulk.lst)
})

# 
# # compare H3K4me1
# dat.raw.pbulk1 <- dat.raw.pbulk.lst.lst$dynamic_bins.50kb$H3K9me3
# dat.raw.pbulk2 <- dat.raw.pbulk.lst.lst$DE_bins_all_marks_top_6085_dists_to_TSS.annot_table$H3K9me3
# 
# # check inputs to matrix of fig6
# 
# bsize <- 50000
# inf.k9 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/DE_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.H3K9me3.2021-01-30.txt"
# dat.k9 <- fread(inf.k9)
# dat.k9 <- dat.k9 %>%
#   rowwise() %>%
#   mutate(startExtend = start + 1 - bsize / 2,
#          endExtend = end - 1 + bsize / 2)
# dat.k9$region_coordExtend <- paste(dat.k9$seqnames, paste(dat.k9$startExtend, dat.k9$endExtend, sep = "-"), sep = ":")
# 
# inf.k9.dynamic.bins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/regions_H3K4me1_H3K9me3_dynamic_regions/H3K4me1_H3K9me3_celltype_specific_genes_and_bins.2021-01-29.txt"
# k9.dynamic.bins <- subset(fread(inf.k9.dynamic.bins), grepl("^chr", V4))$V4
# 
# 
# inf.k9.check <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_differential_expression_outputs/heatmaps_scchix/K9_regions_matrix_input.use_imputed.TRUE.mark.H3K9me3.log2scale.txt")
# dat.k9.check <- read.table(inf.k9.check)
# 
# rcommon <- intersect(rownames(dat.raw.pbulk1), rownames(dat.raw.pbulk2))
# rcommon <- intersect(rownames(dat.raw.pbulk1), rownames(dat.k9.check))
# rcommon <- intersect(rownames(dat.raw.pbulk2), rownames(dat.k9.check))
# rcommon <- intersect(dat.k9$region_coordExtend, sapply(rownames(dat.k9.check), function(x) strsplit(x, ";")[[1]][[2]]))
# rcommon <- intersect(k9.dynamic.bins, sapply(rownames(dat.k9.check), function(x) strsplit(x, ";")[[1]][[2]]))
# 
# length(rcommon)



# Calculate variance for each gene attribute to cluster -------------------

dat.raw.pbulk.lst <- dat.raw.pbulk.lst.lst$DE_bins_all_marks_top_6085_dists_to_TSS.annot_table

head(dat.raw.pbulk.lst[[1]])
head(log2(dat.raw.pbulk.lst[[1]]))

jmat <- log2(dat.raw.pbulk.lst[[1]])
jcheck <- sweep(jmat, MARGIN = 1, STATS = rowMeans(jmat), FUN = "-") 


dat.vars.lst <- lapply(dat.raw.pbulk.lst, function(jdat){
  jdat.log2 <- log2(jdat)
  jcheck <- sweep(jdat.log2, MARGIN = 1, STATS = rowMeans(jdat.log2), FUN = "-") 
  jcheck <- jcheck ^ 2
})


# Calculate variance of each gene  ----------------------------------------

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

dat.vars.long <- dat.vars.long.lst %>%
  bind_rows()

# normalize var
dat.vars.long <- dat.vars.long %>%
  group_by(mark) %>%
  mutate(varfrac = vartotal / sum(vartotal))

ggplot(dat.vars.long, aes(x = forcats::fct_reorder(.f = celltype, .x = vartotal, .fun = median, .desc = TRUE), y = vartotal, fill = clustercol)) + 
  geom_col() + 
  scale_fill_identity() + 
  facet_wrap(~mark, scales = "free_x") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


ggplot(dat.vars.long, aes(x = forcats::fct_reorder(.f = celltype, .x = varfrac, .fun = median, .desc = TRUE), y = varfrac, fill = clustercol)) + 
  geom_col() + 
  scale_fill_identity() + 
  facet_wrap(~mark, scales = "free_x") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


for (jmark in jmarks){
  m <- ggplot(dat.vars.long %>% filter(mark == jmark), aes(x = forcats::fct_reorder(.f = celltype, .x = vartotal, .fun = median, .desc = TRUE), y = vartotal, fill = clustercol)) + 
    geom_col() + 
    ylim(c(0, 250)) + 
    scale_fill_identity() + 
    facet_wrap(~mark, scales = "free_x") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
}

nbins.lst <- lapply(jmarks, function(jmark){
  nrow(dat.raw.pbulk.lst$H3K4me1)
})

for (jmark in jmarks){
  nbins <- nbins.lst[[jmark]]
  m <- ggplot(dat.vars.long %>% filter(mark == jmark), aes(x = forcats::fct_reorder(.f = celltype, .x = vartotal, .fun = median, .desc = TRUE), y = vartotal / nbins, fill = clustercol)) + 
    geom_col() + 
    scale_fill_identity() + 
    ylim(c(0, 0.05)) + 
    facet_wrap(~mark, scales = "free_x") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
}

for (jmark in jmarks){
  m <- ggplot(dat.vars.long %>% filter(mark == jmark), aes(x = forcats::fct_reorder(.f = celltype, .x = varfrac, .fun = median, .desc = TRUE), y = varfrac, fill = clustercol)) + 
    geom_col() + 
    ylim(c(0, 0.5)) + 
    scale_fill_identity() + 
    facet_wrap(~mark, scales = "free_x") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
}


dev.off()