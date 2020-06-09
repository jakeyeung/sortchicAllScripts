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

WrangleRefCtype <- function(zscore.bm.arrows, jclst, jgset){
  zscore.bm.arrows.cspec <- zscore.bm.arrows %>%
    filter(cluster %in% c("HSPCs", jclst)) %>%
    rowwise() %>%
    mutate(geneset = ifelse(geneset == jgset, "CSpecific", "NotCtypeSpecific")) %>%
    group_by(geneset, cluster) %>%
    summarise(H3K4me1 = mean(H3K4me1),
              H3K4me3 = mean(H3K4me3),
              H3K27me3 = mean(H3K27me3))
  zscore.bm.arrows.cspec.wide <- zscore.bm.arrows.cspec %>%
    data.table::melt(., id.vars = c("geneset", "cluster"), measure.vars = c("H3K4me1", "H3K4me3", "H3K27me3"), variable.name = "mark", value.name = "zscore") %>% 
    rowwise() %>%
    mutate(cluster = ifelse(cluster == jclst, "Ctype", "HSPCs")) %>%
    data.table::dcast(., geneset ~ cluster + mark, value.var = "zscore") %>%
    mutate(ref.ctype = jclst)
  return(zscore.bm.arrows.cspec.wide)
}

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

CalculateGeometricMedian <- function(x1, x2, x3, cname.out, cnames = c("H3K4me1", "H3K4me3", "H3K27me3")){
  # create matrix
  X <- as.matrix(data.frame(x1 = x1, x2 = x2, x3 = x3, stringsAsFactors = FALSE))
  colnames(X) <- cnames
  gmed.out <- Gmedian::Gmedian(X)
  colnames(gmed.out) <- cnames
  gmed.out <- as.data.frame(gmed.out)
  return(gmed.out[[cname.out]])
}

ReannotateGeneset <- function(geneset, genesets.ctype, genesets.other, na.fill = NA){
  genesets.ctype.str <- paste(genesets.ctype, collapse = "&")
  genesets.other.str <- paste(genesets.other, collapse = "&")
  if (geneset %in% genesets.ctype){
    gset <- paste0(genesets.ctype.str, "-specificGenes")
  } else if (geneset %in% genesets.other){
    gset <- paste0("z", genesets.other.str, "-specificGenes")
  } else {
    gset <- na.fill
  }
  return(gset)
}


# Constants ---------------------------------------------------------------


make.plots <- FALSE


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks


cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
cbPalette2 <- c("#56B4E9", "#32CD32", "#FFB6C1", "#696969", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
cbPalette2.shift <- c("#32CD32", "#FFB6C1", "#696969", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# rename some clusters
wkm.rename <- hash(c("eryth1", "eryth2", "HSC1", "HSC2", "monocyte"), c("eryth", "eryth", "HSPCs", "HSPCs", "granu"))
bm.rename <- as.list(hash(c("Bcells", "Eryth", "HSCs", "Neutrophils"), c("lymph", "eryth", "HSPCs", "granu")))




pdfout <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/WKM_BM_together.", Sys.Date(), ".2Dcloud.pdf")

if (make.plots){
  pdf(pdfout, useDingbats = FALSE)
}

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

lapply(m.umaps.wkm, function(m){
  mnew <- m + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")
  return(mnew)
})

m.umaps.bm <- lapply(jmarks, function(jmark){
  jdat <- dat.annot.lst.BM[[jmark]]
  ggplot(jdat, aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_manual(values = cbPalette) + ggtitle(jmark)
})
multiplot(m.umaps.bm[[1]], m.umaps.bm[[2]], m.umaps.bm[[3]], cols = 3)

lapply(m.umaps.bm, function(m){
  mnew <- m + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")
  return(mnew)
})


# Plot boxplots -----------------------------------------------------------

m.box.zscore.wkm <- ggplot(jlong.diff.genesets.WKM, aes(x = mark, y = zscore, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette2) + ggtitle("Zebrafish WKM")
print(m.box.zscore.wkm)

m.box.log2fc.wkm <- ggplot(jlong.diff.genesets.WKM, aes(x = mark, y = log2fc, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette2) + ggtitle("Zebrafish WKM")
print(m.box.log2fc.wkm)

m.box.logcounts.wkm <- ggplot(jlong.diff.genesets.WKM, aes(x = mark, y = log2p1counts, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette2) + ggtitle("Zebrafish WKM")
print(m.box.logcounts.wkm)

m.box.counts.wkm <- ggplot(jlong.diff.genesets.WKM, aes(x = mark, y = sqrt(counts), fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette2) + ggtitle("Zebrafish WKM")
print(m.box.counts.wkm)


m.box.zscore.filt.wkm <- ggplot(subset(jlong.diff.genesets.WKM, geneset %in% gsetfilt.WKM), aes(x = mark, y = zscore, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette2) + ggtitle("Zebrafish WKM")
print(m.box.zscore.filt.wkm)

m.box.log2fc.filt.wkm <- ggplot(subset(jlong.diff.genesets.WKM, geneset %in% gsetfilt.WKM), aes(x = mark, y = log2fc, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette2) + ggtitle("Zebrafish WKM")
print(m.box.log2fc.filt.wkm)

m.box.logcounts.filt.wkm <- ggplot(subset(jlong.diff.genesets.WKM, geneset %in% gsetfilt.WKM), aes(x = mark, y = log2p1counts, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette2) + ggtitle("Zebrafish WKM")
print(m.box.logcounts.filt.wkm)

m.box.counts.filt.wkm <- ggplot(subset(jlong.diff.genesets.WKM, geneset %in% gsetfilt.WKM), aes(x = mark, y = sqrt(counts), fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette2) + ggtitle("Zebrafish WKM")
print(m.box.counts.filt.wkm)




m.box.zscore.bm <- ggplot(jlong.diff.genesets.BM, aes(x = mark, y = zscore, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette2) + ggtitle("Mouse BM")
print(m.box.zscore.bm)

m.box.log2fc.bm <- ggplot(jlong.diff.genesets.BM, aes(x = mark, y = log2fc, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette2) + ggtitle("Mouse BM")
print(m.box.log2fc.bm)

m.box.logcounts.bm <- ggplot(jlong.diff.genesets.BM, aes(x = mark, y = log2p1counts, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette2) + ggtitle("Mouse BM")
print(m.box.logcounts.bm)

m.box.counts.bm <- ggplot(jlong.diff.genesets.BM, aes(x = mark, y = sqrt(counts), fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette2) + ggtitle("Mouse BM")
print(m.box.counts.bm)







m.box.zscore.filt.bm <- ggplot(subset(jlong.diff.genesets.BM, geneset %in% gsetfilt.BM), aes(x = mark, y = zscore, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette2) + ggtitle("Mouse BM")
print(m.box.zscore.filt.bm)

m.box.log2fc.filt.bm <- ggplot(subset(jlong.diff.genesets.BM, geneset %in% gsetfilt.BM), aes(x = mark, y = log2fc, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette2) + ggtitle("Mouse BM")
print(m.box.log2fc.filt.bm)

m.box.logcounts.filt.bm <- ggplot(subset(jlong.diff.genesets.BM, geneset %in% gsetfilt.BM), aes(x = mark, y = log2p1counts, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette2) + ggtitle("Mouse BM")
print(m.box.logcounts.filt.bm)

m.box.counts.filt.bm <- ggplot(subset(jlong.diff.genesets.BM, geneset %in% gsetfilt.BM), aes(x = mark, y = sqrt(counts), fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~geneset, ncol = 4) + geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette2) + ggtitle("Mouse BM")
print(m.box.counts.filt.bm)



# Show the 2D cloud of zscores  -------------------------------------------

# create matrix of zscores for each mark, labeled by geneset

zscore.bm.mat <- reshape2::dcast(data = jlong.diff.genesets.BM, formula = bin + ens + geneset + cluster ~ mark, value.var = "zscore")
zscore.wkm.mat <- reshape2::dcast(data = jlong.diff.genesets.WKM, formula = bin + ens + geneset + cluster ~ mark, value.var = "zscore")

# zscore.bm.arrows <- zscore.bm.mat %>%
#   group_by(geneset, cluster) %>%
#   summarise(H3K4me3 = CalculateGeometricMedian(x1 = H3K4me1, x2 = H3K4me3, x3 = H3K27me3, cname.out = "H3K4me3", cnames = c("H3K4me1", "H3K4me3", "H3K27me3")),
#             H3K27me3 = CalculateGeometricMedian(x1 = H3K4me1, x2 = H3K4me3, x3 = H3K27me3, cname.out = "H3K27me3", cnames = c("H3K4me1", "H3K4me3", "H3K27me3")),
#             H3K4me1 = CalculateGeometricMedian(x1 = H3K4me1, x2 = H3K4me3, x3 = H3K27me3, cname.out = "H3K4me1", cnames = c("H3K4me1", "H3K4me3", "H3K27me3")),
#             H3K4me1.start = 0,
#             H3K4me3.start = 0,
#             H3K27me3.start = 0)
# zscore.wkm.arrows <- zscore.wkm.mat %>%
#   group_by(geneset, cluster) %>%
#   summarise(H3K4me3 = CalculateGeometricMedian(x1 = H3K4me1, x2 = H3K4me3, x3 = H3K27me3, cname.out = "H3K4me3", cnames = c("H3K4me1", "H3K4me3", "H3K27me3")),
#             H3K27me3 = CalculateGeometricMedian(x1 = H3K4me1, x2 = H3K4me3, x3 = H3K27me3, cname.out = "H3K27me3", cnames = c("H3K4me1", "H3K4me3", "H3K27me3")),
#             H3K4me1 = CalculateGeometricMedian(x1 = H3K4me1, x2 = H3K4me3, x3 = H3K27me3, cname.out = "H3K4me1", cnames = c("H3K4me1", "H3K4me3", "H3K27me3")),
#             H3K4me1.start = 0,
#             H3K4me3.start = 0,
#             H3K27me3.start = 0)

zscore.bm.arrows <- zscore.bm.mat %>%
  group_by(geneset, cluster) %>%
  summarise(H3K4me3 = mean(H3K4me3, na.rm = TRUE),
            H3K27me3 = mean(H3K27me3, na.rm = TRUE),
            H3K4me1 = mean(H3K4me1, na.rm = TRUE),
            H3K4me1.start = 0,
            H3K4me3.start = 0,
            H3K27me3.start = 0)

zscore.wkm.arrows <- zscore.wkm.mat %>%
  group_by(geneset, cluster) %>%
  summarise(H3K4me3 = mean(H3K4me3, na.rm = TRUE),
            H3K27me3 = mean(H3K27me3, na.rm = TRUE),
            H3K4me1 = mean(H3K4me1, na.rm = TRUE),
            H3K4me1.start = 0,
            H3K4me3.start = 0,
            H3K27me3.start = 0)
  

# calculate arrows in 3D?

m.bm.act <- ggplot(subset(zscore.bm.mat, geneset %in% gsetfilt.BM), aes(x = H3K4me3, y = H3K4me1, color = geneset)) + geom_point(alpha = 0.1) + 
  facet_grid(geneset ~ cluster) +
  # facet_grid(cluster ~ geneset) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  geom_segment(data = subset(zscore.bm.arrows, geneset %in% gsetfilt.BM), mapping = aes(xend = 0, yend = 0), arrow = arrow(length=unit(0.10,"cm"), ends = "first"), color = "black", alpha = 0.5) + 
  ggtitle("Mouse Bone Marrow: active vs active")
print(m.bm.act)

m.bm.rep <- ggplot(subset(zscore.bm.mat, geneset %in% gsetfilt.BM), aes(x = H3K4me3, y = H3K27me3, color = geneset)) + geom_point(alpha = 0.1) + 
  facet_grid(geneset ~ cluster) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  geom_segment(data = subset(zscore.bm.arrows, geneset %in% gsetfilt.BM), mapping = aes(xend = H3K4me3.start, yend = H3K27me3.start), arrow = arrow(length=unit(0.10,"cm"), ends = "first"), color = "black", alpha = 0.5) + 
  ggtitle("Mouse Bone Marrow: active vs repressive")
print(m.bm.rep)

m.bm.rep2 <- ggplot(subset(zscore.bm.mat, geneset %in% gsetfilt.BM), aes(x = H3K4me1, y = H3K27me3, color = geneset)) + geom_point(alpha = 0.1) + 
  facet_grid(geneset ~ cluster) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  geom_segment(data = subset(zscore.bm.arrows, geneset %in% gsetfilt.BM), mapping = aes(xend = H3K4me1.start, yend = H3K27me3.start), arrow = arrow(length=unit(0.10,"cm"), ends = "first"), color = "black", alpha = 0.5) + 
  ggtitle("Mouse Bone Marrow: active vs repressive 2")
print(m.bm.rep2)



m.wkm.act <- ggplot(subset(zscore.wkm.mat, geneset %in% gsetfilt.WKM), aes(x = H3K4me3, y = H3K4me1, color = geneset)) + geom_point(alpha = 0.1) + 
  facet_grid(geneset ~ cluster) +
  # facet_grid(cluster ~ geneset) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  geom_segment(data = subset(zscore.wkm.arrows, geneset %in% gsetfilt.WKM), mapping = aes(xend = 0, yend = 0), arrow = arrow(length=unit(0.10,"cm"), ends = "first"), color = "black", alpha = 0.5) + 
  ggtitle("WKM: active vs active")
print(m.wkm.act)

m.wkm.rep <- ggplot(subset(zscore.wkm.mat, geneset %in% gsetfilt.WKM), aes(x = H3K4me3, y = H3K27me3, color = geneset)) + geom_point(alpha = 0.1) + 
  facet_grid(geneset ~ cluster) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  geom_segment(data = subset(zscore.wkm.arrows, geneset %in% gsetfilt.WKM), mapping = aes(xend = 0, yend = 0), arrow = arrow(length=unit(0.10,"cm"), ends = "first"), color = "black", alpha = 0.5) + 
  ggtitle("WKM: active vs rep")
print(m.wkm.rep)

m.wkm.rep2 <- ggplot(subset(zscore.wkm.mat, geneset %in% gsetfilt.WKM), aes(x = H3K4me1, y = H3K27me3, color = geneset)) + geom_point(alpha = 0.1) + 
  facet_grid(geneset ~ cluster) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  geom_segment(data = subset(zscore.wkm.arrows, geneset %in% gsetfilt.WKM), mapping = aes(xend = 0, yend = 0), arrow = arrow(length=unit(0.10,"cm"), ends = "first"), color = "black", alpha = 0.5) + 
  ggtitle("WKM: active vs rep")
print(m.wkm.rep2)

# collapse these arrows into one plot to summarize everything 

# do this cellltype by celltype (tedious?)



# 3 possible celltypes: granu, lymph, eryth
print(gsetfilt.BM)
print(unique(as.character(zscore.bm.arrows$cluster)))
clstrfilt.BM <- c("HSPCs", "Granulocytes", "Bcells", "Erythroblasts"); names(clstrfilt.BM) <- clstrfilt.BM
jclsts <- clstrfilt.BM[which(clstrfilt.BM != "HSPCs")]

print(gsetfilt.WKM)
print(unique(as.character(zscore.wkm.arrows$cluster)))
clstrfilt.WKM <- c("HSPCs", "granu", "lymph", "eryth"); names(clstrfilt.BM) <- clstrfilt.BM
jclsts.WKM <- clstrfilt.WKM[which(clstrfilt.WKM != "HSPCs")]


zscore.cspec.BM.wide <- lapply(jclsts, function(jclst){
  print(jclst)
  jgset <- gsetfilt.BM[which(clstrfilt.BM  == jclst)]
  zscore.bm.arrows.cspec.wide <- WrangleRefCtype(zscore.bm.arrows, jclst, jgset)
}) %>%
  bind_rows()

zscore.cspec.WKM.wide <- lapply(jclsts.WKM, function(jclst){
  print(jclst)
  jgset <- gsetfilt.WKM[which(clstrfilt.WKM  == jclst)]
  zscore.bm.arrows.cspec.wide <- WrangleRefCtype(zscore.wkm.arrows, jclst, jgset)
}) %>%
  bind_rows() 
zscore.cspec.WKM.wide$ref.ctype <- factor(zscore.cspec.WKM.wide$ref.ctype, levels = c("granu", "lymph", "eryth"))

ggplot(zscore.cspec.BM.wide, aes(x = HSPCs_H3K4me3, xend = Ctype_H3K4me3, y = HSPCs_H3K27me3, yend = Ctype_H3K27me3, color = ref.ctype)) + 
  geom_segment(arrow = arrow(length=unit(0.10,"cm"), ends = "last")) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~geneset) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  xlab("Zscore H3K4me3") + ylab("Zscore H3K27me3")  + ggtitle("Mouse Bone Marrow") + 
  scale_color_manual(values = cbPalette2.shift)

ggplot(zscore.cspec.WKM.wide, aes(x = HSPCs_H3K4me3, xend = Ctype_H3K4me3, y = HSPCs_H3K27me3, yend = Ctype_H3K27me3, color = ref.ctype)) + 
  geom_segment(arrow = arrow(length=unit(0.10,"cm"), ends = "last")) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~geneset) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  xlab("Zscore H3K4me3") + ylab("Zscore H3K27me3")  + ggtitle("ZF WKM") + 
  scale_color_manual(values = cbPalette2.shift)
  

if (make.plots){
  dev.off()
}


# Show flows --------------------------------------------------------------

# show BM active-repressed plane 

logcounts.bm.mat <- reshape2::dcast(data = jlong.diff.genesets.BM, formula = bin + ens + geneset + cluster ~ mark, value.var = "log2p1counts")
logcounts.wkm.mat <- reshape2::dcast(data = jlong.diff.genesets.WKM, formula = bin + ens + geneset + cluster ~ mark, value.var = "log2p1counts")


counts.bm.mat <- reshape2::dcast(data = jlong.diff.genesets.BM, formula = bin + ens + geneset + cluster ~ mark, value.var = "counts")
counts.wkm.mat <- reshape2::dcast(data = jlong.diff.genesets.WKM, formula = bin + ens + geneset + cluster ~ mark, value.var = "counts")

ggplot(logcounts.bm.mat, aes(x = H3K4me3, y = H3K27me3)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(counts.bm.mat, aes(x = H3K4me3, y = H3K27me3)) + 
  geom_point(alpha = 0.1) + 
  geom_density_2d() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cluster)

# show the flows
cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
jscale <- 5
jctype.init <- "HSPCs"

# do it in a loop

print(unique(as.character(counts.bm.mat$cluster)))
print(unique(as.character(counts.bm.mat$geneset)))

jctype <- "Granulocytes"
genesets.ctype <- "Neutrophil"
genesets.other <- c("Bcell", "Erythroblast")

jctype.vec <- list("Granulocytes", "Bcells", "Erythroblasts")
genesets.ctype.vec <- list("Neutrophil", "Bcell", "Erythroblast")
genesets.other.vec <- list(c("Bcell", "Erythroblast"), c("Neutrophil", "Erythroblast"), c("Neutrophil", "Bcell"))
names(jctype.vec) <- jctype.vec
names(genesets.ctype.vec) <- jctype.vec
names(genesets.other.vec) <- jctype.vec


pdfout <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/gene_flows/celltype_specific_gene_flows_HSC_to_celltype.", Sys.Date(), ".with_boxplot.pdf")
pdf(pdfout, width = 1020/72, height = 815/72, useDingbats = FALSE)
for (jctype in jctype.vec){
  
  print(jctype)
  
  genesets.ctype <- genesets.ctype.vec[[jctype]]
  genesets.other <- genesets.other.vec[[jctype]]
  # genesets.other <- c("Bcell", "Erythroblast")
  
  
  jtitle <- paste0(jctype.init, '->', jctype)
  
  # make endspoints granulocytes
  
  init.mat <- subset(counts.bm.mat, cluster == jctype.init)
  final.mat <- subset(counts.bm.mat, cluster %in% c(jctype)) %>%
    dplyr::rename(H3K4me1_end = H3K4me1, H3K4me3_end = H3K4me3, H3K27me3_end = H3K27me3)
  final.mat$cluster <- NULL
  init_final.mat <- left_join(init.mat, final.mat, by = c("bin", "ens", "geneset"))
  
  # reannotate geneset
  # init_final.mat <- subset(init_final.mat, geneset %in% c(genesets.ctype, genesets.other))
  
  # check genesets
  # get magnitude in log2fc or zscore
  log2fc.hsc_to_ctype.long <- jlong.diff.genesets.BM %>%
    filter(cluster %in% c(jctype.init, jctype)) %>%
    group_by(bin, gene, ens, geneset, mark) %>%
    summarise(log2fc = log2p1counts[[1]] - log2p1counts[[2]])
  
  log2fc.hsc_to_ctype.wide <- log2fc.hsc_to_ctype.long %>%
    reshape2::dcast(data = ., formula = "bin + ens + geneset ~ mark", value.var = "log2fc") %>%
    dplyr::rename(H3K4me1_log2fc = H3K4me1, H3K4me3_log2fc = H3K4me3, H3K27me3_log2fc = H3K27me3)
  
  init_final.log2fc.mat <- left_join(init_final.mat, log2fc.hsc_to_ctype.wide, by = c("bin", "ens", "geneset")) %>%
    rowwise() %>%
    mutate(gset = ReannotateGeneset(geneset, genesets.ctype, genesets.other, na.fill = NA))
  
  init_final.log2fc.mat <- subset(init_final.log2fc.mat, !is.na(gset))
  
  
  m.density_log2fc.active <- ggplot(init_final.log2fc.mat, aes(x = H3K4me3_log2fc, fill = gset)) + 
    geom_density(alpha = 0.3) + 
    # facet_wrap(~gset) + 
    geom_vline(xintercept = 0) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle, "across different celltype-specific genes")
  
  m.density_log2fc.repress <- ggplot(init_final.log2fc.mat, aes(x = H3K27me3_log2fc, fill = gset)) + 
    geom_density(alpha = 0.3) + 
    # facet_wrap(~gset) + 
    geom_vline(xintercept = 0) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(jtitle, "across different celltype-specific genes")
  
  
  
  m.phase.log2fc.full <- ggplot(data = init_final.log2fc.mat, aes(x = H3K4me3, y = H3K27me3, color = gset)) + 
    geom_point(alpha = 0.1, size = 0.1, color = 'black') +  
    geom_segment(mapping = aes(xend = H3K4me3 + jscale*H3K4me3_log2fc, yend = H3K27me3 + jscale*H3K27me3_log2fc),
                 arrow = arrow(length=unit(0.25,"cm"), ends = "last"), alpha = 0.8, size = 0.1) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
    ggtitle(jtitle)
  
  
  m.phase.log2fc.zoom <- ggplot(data = init_final.log2fc.mat, aes(x = H3K4me3, y = H3K27me3, color = gset)) + 
    geom_point(alpha = 0.1, size = 0.1, color = 'black') +  
    geom_segment(mapping = aes(xend = H3K4me3 + jscale*H3K4me3_log2fc, yend = H3K27me3 + jscale*H3K27me3_log2fc),
                 arrow = arrow(length=unit(0.25,"cm"), ends = "last"), alpha = 0.8, size = 0.1) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
    coord_cartesian(ylim=c(0, 50), xlim = c(0, 100)) + facet_wrap(~gset) + 
    ggtitle(jtitle, "zoomed in")
  
  
  print(m.phase.log2fc.full)
  print(m.phase.log2fc.zoom)
  print(m.density_log2fc.active)
  print(m.density_log2fc.repress)
}
dev.off()









