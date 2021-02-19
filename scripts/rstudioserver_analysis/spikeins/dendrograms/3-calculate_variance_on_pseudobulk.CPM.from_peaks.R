# Jake Yeung
# Date of Creation: 2021-02-19
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/dendrograms/3-calculate_variance_on_pseudobulk.CPM.from_peaks.R
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

jymax <- 500

# merge.ctypes.by.lineage <- FALSE
# merge.ctypes.by.lineage <- FALSE

jscale <- 1000
log2filt <- -0.75

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps"
outpdf <- file.path(outdir, paste0("variance_by_clusters_by_peaks.", Sys.Date(), ".", jscale, ".log2filt_", log2filt, ".pdf"))
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

cname2color <- hash::hash(dat.metas$H3K4me1$cluster, dat.metas$H3K4me1$clustercol)
colpalette <- colorRampPalette(c("grey1", "grey35", "grey99"))(1024)

# Load mats ---------------------------------------------------------------


# jmark <- "H3K9me3"
niter <- "100"
inf.glmpca.lst <- lapply(jmarks, function(jmark){
  if (jmark != "H3K27me3"){
    inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/BM_from_matadj/glmpca.", jmark, ".from_matadj.platename_jrep.szname_none.niter_", niter, ".RData"))
    assertthat::assert_that(file.exists(inf.glmpca))
  } else {
    inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/H3K27me3_rep2rep3reseq.peaks.varfilt/glmpca.", jmark, ".bincutoff_0.binskeep_0.platename_plate.szname_none.niter_500.reorder_rownames.dupfilt.suffix_peaks.RData"))
  }
  return(inf.glmpca)
})


count.mat.peaks.raw.lst <- lapply(inf.glmpca.lst, function(jinf){
  load(jinf, v=T)
  return(glm.inits$Y.filt)
})

count.mat.peaks.raw.filt.lst <- lapply(count.mat.peaks.raw.lst, function(count.mat){
  rmeans <- log10(rowMeans(count.mat))
  bfilt.keep <- rmeans > log2filt
  count.mat[bfilt.keep, ]
})

lapply(count.mat.peaks.raw.lst, dim)
lapply(count.mat.peaks.raw.filt.lst, dim)


# Load LDA  ---------------------------------------------------------------

dat.raw.pbulk.lst <- lapply(jmarks, function(jmark){
  count.mat <- count.mat.peaks.raw.filt.lst[[jmark]]
  cnames.keep.lst <- split(x = dat.metas[[jmark]]$cell, f = dat.metas[[jmark]]$cluster)
  dat.raw.pbulk <- SumAcrossClusters(count.mat, cnames.keep.lst)
  dat.raw.pbulk <- do.call(cbind, dat.raw.pbulk)
  # normalize
  dat.raw.pbulk <- sweep(dat.raw.pbulk, MARGIN = 2, STATS = colSums(dat.raw.pbulk), FUN = "/", check.margin = TRUE) * jscale + 1
})


# Calculate variance for each gene attribute to cluster -------------------

# dat.raw.pbulk.lst <- dat.raw.pbulk.lst.lst$DE_bins_all_marks_top_6085_dists_to_TSS.annot_table

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
    ylim(c(0, jymax)) + 
    scale_fill_identity() + 
    facet_wrap(~mark, scales = "free_x") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
}

nbins.lst <- lapply(jmarks, function(jmark){
  nrow(dat.raw.pbulk.lst[[jmark]])
})

for (jmark in jmarks){
  nbins <- nbins.lst[[jmark]]
  m <- ggplot(dat.vars.long %>% filter(mark == jmark), aes(x = forcats::fct_reorder(.f = celltype, .x = vartotal, .fun = median, .desc = TRUE), y = vartotal / nbins, fill = clustercol)) + 
    geom_col() + 
    scale_fill_identity() + 
    ylim(c(0, 1.5)) + 
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

# split variance by pos or neg


dat.diffs.lst <- lapply(dat.raw.pbulk.lst, function(jdat){
  # jdat <- jout$dat.impute.pbulk
  # jdat <- jout$dat.raw.pbulk
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

# make heatmaps 

# jmark <- "H3K9me3"
jmethod <- "ward.D2"
# jmethod <- "complete"
for (jmark in jmarks){
  jmat2 <- log2(dat.raw.pbulk.lst[[jmark]])
  cnames.color <- sapply(colnames(jmat2), function(x) AssignHash(x, cname2color, null.fill = x))
  heatmap3::heatmap3(jmat2, 
                     Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, 
                     main = paste("peaks", jmark, jmethod), margins = c(5, 8), cexRow = 0.5, method = jmethod, ColSideColors = cnames.color, col = colpalette)
}


dev.off()
