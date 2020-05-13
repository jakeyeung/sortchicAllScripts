# Jake Yeung
# Date of Creation: 2020-05-10
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/pretty_plots_WKM_BM.R
# Load objects from previous analyses and make plots that are coherent 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(hash)
library(JFuncs)
library(scchicFuncs)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks


cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# rename some clusters
wkm.rename <- hash(c("eryth1", "eryth2", "HSC1", "HSC2", "monocyte"), c("eryth", "eryth", "HSPCs", "HSPCs", "granu"))
bm.rename <- as.list(hash(c("Bcells", "Eryth", "HSCs", "Neutrophils"), c("lymph", "eryth", "HSPCs", "granu")))

RenameClusterBM <- function(clstr.orig, bm.rename){
  # clstr.orig <- "Bcells-Cd83_topic10"
  clstr.new <- paste0("z", clstr.orig)
  for (cname in names(bm.rename)){
    if (startsWith(clstr.orig, prefix = cname)){
      clstr.new <- bm.rename[[cname]]
    } else{
    }
  }
  return(clstr.new)
}


pdfout <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/WKM_BM_together.", Sys.Date(), ".pdf")

pdf(pdfout, useDingbats = FALSE)

# Load UMAPs for BM and WKM  ----------------------------------------------


dat.annot.lst.WKM <- lapply(jmarks, function(jmark){
  print(jmark)
  # filter by previously defined celltypes? 
  inf.annot.louv <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/LDA_downstream/LDA_downstream_ZF.2020-04-23.imputevarfilt.lessstringent/ZF_LDA_output.", jmark, ".keepn_150.final.ClusterTables.txt")
  assertthat::assert_that(file.exists(inf.annot.louv))
  
  
  inf.annot.glmpca <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/celltyping/from_louvain/cell_to_cluster_table.", jmark, ".2020-04-13.txt")
  assertthat::assert_that(file.exists(inf.annot.glmpca))
  
  annot.louv <- fread(inf.annot.louv)
  annot.louv$clusterplate <- paste(annot.louv$cluster, annot.louv$plate, "_")
  
  annot.glmpca <- fread(inf.annot.glmpca)
  annot.glmpca.filt <- subset(annot.glmpca, cell %in% annot.louv$cell) %>%
    rowwise() %>%
    mutate(clusterplate = paste(cluster, plate, sep = "_")) %>%
    mutate(cluster = ifelse(cluster %in% c("lymph1", "lymph2"), "lymph", cluster)) %>%   # rename lymph1 and lymph2 into lymph
    ungroup() %>%
    filter(cluster != "Unknown") %>%
    left_join(., subset(annot.louv, select = c(cell, var.imputed)))
  # rename clusters
  annot.glmpca.filt$cluster <- sapply(annot.glmpca.filt$cluster, function(jclst) AssignHash(jclst, wkm.rename, null.fill = jclst))
  print("annot glmpca filt")
  print(annot.glmpca.filt)
  return(annot.glmpca.filt)
})


dat.annot.lst.BM <- lapply(jmarks, function(jmark){
  inf.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
  assertthat::assert_that(file.exists(inf.annot))
  load(inf.annot, v=T)
  dat.umap.glm.fillNAs <- subset(dat.umap.glm.fillNAs, !is.na(cluster))
  dat.umap.glm.fillNAs$cluster <- sapply(dat.umap.glm.fillNAs$cluster, RenameClusterBM, bm.rename)
  return(dat.umap.glm.fillNAs)
})


# Load pseudobulks  -------------------------------------------------------

# zebrafish WKM

jwinsize <- 10000L
downsample.reads <- TRUE
downsample.cells <- FALSE
jdate <- "2020-05-09"  # downsample.cells FALSE now does not down sample 
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/pseudobulk_analysis.integrated_analysis.TSS.readsDS.likeBM.OtherWinsizes"

# inprefix <- "integrated_analysis.2020-05-09.UseTSSfromH3K4me3.DSreads_TRUE.DScells_FALSE.likeBM.winsize_10000.RData"
inprefix <- paste0("integrated_analysis.", jdate, ".UseTSSfromH3K4me3.DSreads_", downsample.reads, ".DScells_", downsample.cells, ".likeBM.winsize_", jwinsize)
inrdata <- paste0(inprefix, ".RData")
infrdata <- file.path(indir, inrdata)

assertthat::assert_that(file.exists(infrdata))
load(infrdata, v=T)
jlong.diff.genesets.WKM <- jlong.diff.genesets
jlong.diff.genesets.WKM$cluster <- sapply(jlong.diff.genesets.WKM$cluster, function(x) gsub("monocyte", "granu", x = x))
jlong.diff.genesets.WKM$cluster <- sapply(jlong.diff.genesets.WKM$cluster, function(x) gsub("HSC", "HSPCs", x = x))
jlong.diff.genesets.WKM$cluster <- factor(jlong.diff.genesets.WKM$cluster, levels = c("HSPCs", "granu", "lymph", "eryth"))
jlong.diff.genesets.WKM$mark <- factor(jlong.diff.genesets.WKM$mark, jmarks)
jlong.diff.genesets.WKM$geneset <- factor(jlong.diff.genesets.WKM$geneset, c("HSPCs", "granulocytes", "lymphocytes", "erythrocytes", "HighExprs", "LowExprs", "zOther"))
gsetfilt.WKM <- c("HSPCs", "granulocytes", "lymphocytes", "erythrocytes")
  
# mouse BM

inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.stringentDE/integrated_analysis_3_marks.stringentDE.WithHighLowExprs.pseudobulk.fewerk27me3_TRUE.DSreads_TRUEDScells.FALSE.2020-05-04.RData"
load(inf, v=T)
jlong.diff.genesets.BM <- jlong.diff.genesets
print(unique(jlong.diff.genesets.BM$cluster))
jlong.diff.genesets.BM$cluster <- factor(jlong.diff.genesets.BM$cluster, c("HSPCs", "Granulocytes", "Bcells", "Erythroblasts"))
jlong.diff.genesets.BM$mark <- factor(jlong.diff.genesets.BM$mark, jmarks)
jlong.diff.genesets.BM$geneset <- factor(jlong.diff.genesets.BM$geneset, c("HSCs", "Neutrophil", "Bcell", "Erythroblast", "HighExprs", "LowExprs", "zOther"))
gsetfilt.BM <- c("HSCs", "Neutrophil", "Bcell", "Erythroblast")


# Plot UMAPs --------------------------------------------------------------

print(lapply(dat.annot.lst.BM, function(x) sort(unique(x$cluster))))


m.umaps.wkm <- lapply(jmarks, function(jmark){
  jdat <- dat.annot.lst.WKM[[jmark]]
  ggplot(jdat, aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_manual(values = cbPalette) + ggtitle(jmark)
})
multiplot(m.umaps.wkm[[1]], m.umaps.wkm[[2]], m.umaps.wkm[[3]], cols = 3)

m.umaps.bm <- lapply(jmarks, function(jmark){
  jdat <- dat.annot.lst.BM[[jmark]]
  ggplot(jdat, aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_manual(values = cbPalette) + ggtitle(jmark)
})
multiplot(m.umaps.bm[[1]], m.umaps.bm[[2]], m.umaps.bm[[3]], cols = 3)


# Plot boxplots -----------------------------------------------------------

m.box.zscore.wkm <- ggplot(jlong.diff.genesets.WKM, aes(x = mark, y = zscore, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette) + ggtitle("Zebrafish WKM")
print(m.box.zscore.wkm)

m.box.log2fc.wkm <- ggplot(jlong.diff.genesets.WKM, aes(x = mark, y = log2fc, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette) + ggtitle("Zebrafish WKM")
print(m.box.log2fc.wkm)

m.box.logcounts.wkm <- ggplot(jlong.diff.genesets.WKM, aes(x = mark, y = log2p1counts, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette) + ggtitle("Zebrafish WKM")
print(m.box.logcounts.wkm)

m.box.counts.wkm <- ggplot(jlong.diff.genesets.WKM, aes(x = mark, y = sqrt(counts), fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette) + ggtitle("Zebrafish WKM")
print(m.box.counts.wkm)


m.box.zscore.filt.wkm <- ggplot(subset(jlong.diff.genesets.WKM, geneset %in% gsetfilt.WKM), aes(x = mark, y = zscore, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette) + ggtitle("Zebrafish WKM")
print(m.box.zscore.filt.wkm)

m.box.log2fc.filt.wkm <- ggplot(subset(jlong.diff.genesets.WKM, geneset %in% gsetfilt.WKM), aes(x = mark, y = log2fc, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette) + ggtitle("Zebrafish WKM")
print(m.box.log2fc.filt.wkm)

m.box.logcounts.filt.wkm <- ggplot(subset(jlong.diff.genesets.WKM, geneset %in% gsetfilt.WKM), aes(x = mark, y = log2p1counts, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette) + ggtitle("Zebrafish WKM")
print(m.box.logcounts.filt.wkm)

m.box.counts.filt.wkm <- ggplot(subset(jlong.diff.genesets.WKM, geneset %in% gsetfilt.WKM), aes(x = mark, y = sqrt(counts), fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette) + ggtitle("Zebrafish WKM")
print(m.box.counts.filt.wkm)




m.box.zscore.bm <- ggplot(jlong.diff.genesets.BM, aes(x = mark, y = zscore, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette) + ggtitle("Mouse BM")
print(m.box.zscore.bm)

m.box.log2fc.bm <- ggplot(jlong.diff.genesets.BM, aes(x = mark, y = log2fc, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette) + ggtitle("Mouse BM")
print(m.box.log2fc.bm)

m.box.logcounts.bm <- ggplot(jlong.diff.genesets.BM, aes(x = mark, y = log2p1counts, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette) + ggtitle("Mouse BM")
print(m.box.logcounts.bm)

m.box.counts.bm <- ggplot(jlong.diff.genesets.BM, aes(x = mark, y = sqrt(counts), fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette) + ggtitle("Mouse BM")
print(m.box.counts.bm)







m.box.zscore.filt.bm <- ggplot(subset(jlong.diff.genesets.BM, geneset %in% gsetfilt.BM), aes(x = mark, y = zscore, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette) + ggtitle("Mouse BM")
print(m.box.zscore.filt.bm)

m.box.log2fc.filt.bm <- ggplot(subset(jlong.diff.genesets.BM, geneset %in% gsetfilt.BM), aes(x = mark, y = log2fc, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette) + ggtitle("Mouse BM")
print(m.box.log2fc.filt.bm)

m.box.logcounts.filt.bm <- ggplot(subset(jlong.diff.genesets.BM, geneset %in% gsetfilt.BM), aes(x = mark, y = log2p1counts, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette) + ggtitle("Mouse BM")
print(m.box.logcounts.filt.bm)

m.box.counts.filt.bm <- ggplot(subset(jlong.diff.genesets.BM, geneset %in% gsetfilt.BM), aes(x = mark, y = sqrt(counts), fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette) + ggtitle("Mouse BM")
print(m.box.counts.filt.bm)



dev.off()
