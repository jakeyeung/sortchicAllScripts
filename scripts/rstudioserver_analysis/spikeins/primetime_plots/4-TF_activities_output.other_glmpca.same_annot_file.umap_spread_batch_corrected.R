# Jake Yeung
# Date of Creation: 2020-12-28
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/4-TF_activities_output.other_glmpca.same_annot_file.umap_spread_batch_corrected.R
# description

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

library(ggrepel)

library(heatmap3)


source("~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/heatmap3revR.R")


# bwpalette <- colorRampPalette(c("grey99", "grey65", "grey1"))(1024)
bwpalette <- colorRampPalette(c("grey1", "grey75", "grey99"))(1024)
ctypes.arranged <-  c("Eryths", "Bcells", "NKs", "pDCs", "Granulocytes", "DCs", "Basophils", "HSPCs")


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")


# Load UMAP  ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"
# jmark <- "H3K27me3"
jmark <- "H3K4me1"

indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28")

inf.meta <- file.path(indir.meta, paste0("metadata_batch_corrected.", jmark, ".2020-12-28.txt"))
print(inf.meta)
dat.umap <- fread(inf.meta)

# Load clusters -----------------------------------------------------------

# dat.annot.ref <- dat.umap
dat.merge <- as.data.frame(dat.umap)
rownames(dat.merge) <- dat.merge$cell


outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/batch_corrected"

outbase <- paste0("motif_activity_", jmark, ".", Sys.Date(), ".uniquecolor.same_annot_file.fix_colors.umap_spread_batch_corrected")
outname <- paste0(outbase, ".pdf")
outnametxt <- paste0(outbase, ".txt")
outtxt <- file.path(outdir, outnametxt)
outpdf <- file.path(outdir, outname)
pdf(file = outpdf, useDingbats = FALSE)

# cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
ggplot(dat.merge, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load MARA output --------------------------------------------------------

binskeep <- 0
niter <- "500"
jsuffix2 <- paste0("BM_", jmark, ".BM_AllMerged3.glmpca_plate.bincutoff_0.binskeep_", binskeep, ".byplate.szname_none.niter_", niter, ".reorder_rownames.dupfilt.cleanuprows.same_annot_file")

indir.mara <- file.path(hubprefix, paste0("jyeung/data/scChiC/mara_analysis_BM-AllMerged3_Peaks/", jmark, "/mara_output.lessmemory/", jsuffix2, "-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, "/", jsuffix2))
print(indir.mara)
assertthat::assert_that(dir.exists(indir.mara))


mara.out <- LoadMARA(indir.mara, make.cnames = FALSE)

act.mat.clean <- t(subset(mara.out$act.mat, select = -motif))
colnames(act.mat.clean) <- mara.out$act.mat$motif
act.mat.clean.dat <- data.frame(cell = rownames(act.mat.clean), act.mat.clean, stringsAsFactors = FALSE) %>% 
  ungroup() %>%
  mutate(cell = gsub("\\.", "-", cell))


zscores.cutoff <- 0.7
motifs.keep <- subset(mara.out$zscores, zscore > zscores.cutoff)$motif

m.zscores <- mara.out$zscores %>%
  ungroup() %>% 
  mutate(rnk = seq(length(motif)),
         motiflab = ifelse(zscore > zscores.cutoff * 1.5, toupper(motif), NA)) %>%
  ggplot(., aes(x = rnk, y = zscore, label = motiflab)) + 
    geom_point() + 
    geom_text_repel() + 
    theme_bw() + 
    geom_hline(yintercept = zscores.cutoff, linetype = "dotted") + 
    xlab("Rank") + 
    ylab("Motif zscore") + 
    theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.zscores)

m.zscores2 <- mara.out$zscores %>%
  ungroup() %>% 
  mutate(rnk = seq(length(motif))) %>%
  ggplot(., aes(x = zscore)) + 
  geom_density(fill = "red", alpha = 0.25) + 
  geom_vline(xintercept = zscores.cutoff, linetype = "dotted") + 
  theme_bw() + 
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.zscores2)


# Plot zscores ------------------------------------------------------------




# Heatmap -----------------------------------------------------------------



# Do heatmap --------------------------------------------------------------

K <- 6
jmeth <- "ward.D2"

jmat <- as.matrix(subset(act.mat.clean.dat, select = -cell))
rownames(jmat) <- act.mat.clean.dat$cell


# order cells by cluster
ctypes <- unique(dat.merge$cluster)
# collapse some clusters into one
hscs <- list(ctypes[grepl("^HSP", ctypes)])
names(hscs) <- "HSPCs"
ilcs <- list(ctypes[grepl("^NK", ctypes)])
names(ilcs) <- "NKcells"
bcells <- list(ctypes[grepl("^Bcell", ctypes)])
names(bcells) <- "Bcells"
granus <- list(ctypes[grepl("^Granu", ctypes)])
names(granus) <- "Granulocytes"
basos <- list(ctypes[grepl("^Basophils", ctypes)])
names(basos) <- "Basophils"
dends <- list(ctypes[grepl("^DCs", ctypes)])
names(dends) <- c("Dendritic")
pdends <- list(ctypes[grepl("^pDCs", ctypes)])
names(pdends) <- c("pDendritic")

clstr.hash <- hash()
for (lst in list(hscs, ilcs, bcells, granus, basos, dends, pdends)){
  # print(lst)
  jname <- names(lst)[[1]]
  # print(jname)
  # print(lst[[1]])
  for (ctype in lst[[1]]){
    print(paste(ctype, jname))
    clstr.hash[[ctype]] <- jname
  }
}

# dat.annots.mergeclst <- dat.merge %>%
#   rowwise() %>%
#   mutate(clst.merged = ifelse(!is.null(clstr.hash[[cluster]]), clstr.hash[[cluster]], cluster),
#          umapdist = sqrt(umap1^2 + umap2^2)) %>%
#   ungroup()

dat.annots.mergeclst <- dat.merge %>%
  rowwise() %>%
  mutate(clst.merged = cluster,
         umapdist = sqrt(umap1^2 + umap2^2)) %>%
  ungroup()

dat.annots.mergeclst$clst.merged <- factor(dat.annots.mergeclst$clst.merged, levels = ctypes.arranged)

dat.annots.mergeclst <- dat.annots.mergeclst %>%
  mutate(clst.merged = forcats::fct_relevel(clst.merged, ctypes.arranged)) %>%
  arrange(clst.merged)


cells.ordered <- dat.annots.mergeclst$cell
jsub <- jmat[cells.ordered, motifs.keep]


dat.col <- dat.merge[cells.ordered, ] %>%
  dplyr::select(c(cell, clustercol))
colvec <- dat.col$clustercol

# plot color

dat.merge.col <- left_join(dat.merge, dat.col)

ggplot(dat.merge.col, aes(x = umap1, y = umap2, color = clustercol)) + 
  geom_point() + 
  scale_color_identity() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


par(mfrow=c(1,1), mar=c(1,1,1,1), mgp=c(3, 1, 0), las=0)
hm.out <- heatmap3(jsub, margins = c(5, 8), cexCol = 0.35, Colv = TRUE, Rowv = NA, 
                   col = bwpalette, 
                   # ColSideColors = rep("blue", ncol(jsub)), 
                   # ColSideColors = FALSE,
                   RowSideColors = colvec, 
                   # RowSideColors = rep("red", nrow(jsub)), 
                   RowSideLabs = "celltype", 
                   labRow = FALSE, scale = "column", revC = TRUE,
                   distfun = dist, hclustfun = hclust, method = jmeth)

hm.out <- heatmap3(jsub, margins = c(5, 8), cexCol = 0.35, Colv = TRUE, Rowv = NA, 
                   col = bwpalette,
                   # ColSideColors = rep("blue", ncol(jsub)), 
                   # ColSideColors = FALSE,
                   RowSideColors = colvec, 
                   # RowSideColors = rep("red", nrow(jsub)), 
                   RowSideLabs = "celltype", 
                   labRow = FALSE, scale = "column", revC = FALSE,
                   distfun = dist, hclustfun = hclust, method = jmeth)


hm.out.transpose <- heatmap3(t(jsub), margins = c(5, 8), cexCol = 0.35, Colv = NA, Rowv = TRUE, 
                             col = bwpalette,
                   # ColSideColors = rep("blue", ncol(jsub)), 
                   # ColSideColors = FALSE,
                   ColSideColors = colvec, 
                   # RowSideColors = rep("red", nrow(jsub)), 
                   ColSideLabs = "celltype", 
                   labCol = FALSE, scale = "row", revC = FALSE, 
                   distfun = dist, hclustfun = hclust, method = jmeth)

hm.out.transpose <- heatmap3(t(jsub), margins = c(5, 8), cexCol = 0.35, Colv = NA, Rowv = TRUE, 
                             col = bwpalette,
                             # ColSideColors = rep("blue", ncol(jsub)), 
                             # ColSideColors = FALSE,
                             ColSideColors = colvec, 
                             # RowSideColors = rep("red", nrow(jsub)), 
                             ColSideLabs = "celltype", 
                             labCol = FALSE, scale = "row", revC = TRUE,
                             distfun = dist, hclustfun = hclust, method = jmeth)



hm.out.transpose <- heatmap.3(t(jsub), cexCol = 0.35, Colv = NA, Rowv = TRUE, 
                              # col = bwpalette,
                   # ColSideColors = rep("blue", ncol(jsub)), 
                   # ColSideColors = FALSE,
                   ColSideColors = colvec, 
                   # RowSideColors = rep("red", nrow(jsub)), 
                   ColSideLabs = "celltype", 
                   labCol = FALSE, scale = "row", revC = FALSE, revR = TRUE,
                   distfun = dist, hclustfun = hclust, method = jmeth)



# Motif by motif ----------------------------------------------------------



jmotif <- "Cebpb"
jmotif <- "Erg"
jmotif <- "Ebf1"
jmotif <- "Tal1"
jmotif <- "Gata1"
jmotif <- "Irf1"
jmotif <- "Hoxb5"
jmotif <- "Runx1"
jmotif <- "Hoxc6"
jmotif <- "Hlf"
jmotif <- "Bcl3"
jmotif <- "Yy1"
jmotif <- "Hoxa2"


for (jmotif in motifs.keep){
  jzscore <- signif(subset(mara.out$zscores, motif == jmotif)$zscore, digits = 2)
  (jtitle <- paste(jmotif, "Zscore:", jzscore))
  
  dat.merge.motifs <- left_join(dat.merge, act.mat.clean.dat, by = "cell")
  # pdf("/home/jyeung/hub_oudenaarden/jyeung/tmp/motiftest.pdf", useDingbats = FALSE)
  m <- ggplot(dat.merge.motifs, aes_string(x = "umap1", y = "umap2", color = jmotif)) + 
    geom_point(size = 0.75) + 
    ggtitle(jtitle) + 
    theme_bw() + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  # m <- PlotXYWithColor(dat.merge.motifs, xvar = "umap1", yvar = "umap2", cname = jmotif, jtitle = jtitle, cont.color = TRUE, jsize = 0.5) + scale_color_viridis_c()
  print(m)
  # dev.off()
}


# save output with color name
fwrite(dat.merge.col, file = outtxt, quote = FALSE, sep = "\t", col.names = TRUE)

dev.off()


