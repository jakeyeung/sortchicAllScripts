# Jake Yeung
# Date of Creation: 2020-12-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/4-TF_activities_output.other_glmpca.same_annot_file.umap_spread.R
# 


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

library(heatmap3)


source("~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/heatmap3revR.R")

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")


# Load UMAP  ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"
# jmark <- "H3K27me3"
jmark <- "H3K4me1"

indir <- "jyeung/data/scChiC/glmpca_outputs/same_annot_file_rerun"

# niter <- "1000"
binskeep <- 0
niter <- "500"
# binskeep <- 1000
jsuffix <- paste0("glmpca.", jmark, ".bincutoff_0.binskeep_", binskeep, ".byplate.szname_none.niter_", niter, ".reorder_rownames.dupfilt")
fname <- paste0(jsuffix, ".RData")

# jsuffix2 <- paste0("glmpca_plate.bincutoff_0.binskeep_0.byplate.szname_none.niter_1000.reorder_rownames.dupfilt.cleanuprows")
jsuffix2 <- paste0("BM_", jmark, ".BM_AllMerged3.glmpca_plate.bincutoff_0.binskeep_", binskeep, ".byplate.szname_none.niter_", niter, ".reorder_rownames.dupfilt.cleanuprows.same_annot_file")
# jsuffix2 <- paste0(jsuffix, ".cleanuprows")

inf <- file.path(hubprefix, indir, fname)

load(inf, v=T)

#inf.umap <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.", jmark, "_cleaned/cell_cluster_table_with_spikeins.", jmark, ".2020-12-23.umap_spread.final.txt")
inf.umap <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned/cell_cluster_table_with_spikeins.", jmark, ".2020-12-23.umap_spread.final.txt")
assertthat::assert_that(file.exists(inf.umap))
dat.umap <- fread(inf.umap)

# dat.umap <- DoUmapAndLouvain(glm.out$factors, jsettings = jsettings)


# Load clusters -----------------------------------------------------------

# dat.annot.ref <- dat.umap
dat.merge <- as.data.frame(dat.umap)
rownames(dat.merge) <- dat.merge$cell

# jmarkref <- "H3K4me1"  # has all clusters so colors will match
# inf.annot <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.", jmark, ".txt"))
# inf.annot.ref <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.", jmarkref, ".txt"))
# dat.annot <- fread(inf.annot) %>%
#   as.data.frame()
# rownames(dat.annot) <- dat.annot$cell
# 
# dat.annot.ref <- fread(inf.annot.ref)
# 
# cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
# clusters <- sort(unique(dat.annot$cluster))
# colsvec.byclusters <- hash(clusters, cbPalette[seq(length(clusters))])
# dat.annot$col <- sapply(dat.annot$cluster, function(x) AssignHash(x = x, jhash = colsvec.byclusters, null.fill = "grey95"))
# 
# 
# dat.merge <- left_join(dat.umap, subset(dat.annot, select = c(cell, cluster, batch)))

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2"

outbase <- paste0("motif_activity_", jmark, ".niter_", niter, ".binskeep_", binskeep, ".", Sys.Date(), ".uniquecolor.same_annot_file.fix_colors.umap_spread")
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

dat.annots.mergeclst <- dat.merge %>%
  rowwise() %>%
  mutate(clst.merged = ifelse(!is.null(clstr.hash[[cluster]]), clstr.hash[[cluster]], cluster),
         umapdist = sqrt(umap1^2 + umap2^2)) %>%
  ungroup()

dat.annots.mergeclst$clst.merged <- factor(dat.annots.mergeclst$clst.merged, levels = c("Eryths", "Bcells", "NKcells", "Granulocytes", "Basophils", "Dendritic", "pDendritic", "HSPCs"))

dat.annots.mergeclst <- dat.annots.mergeclst %>%
  mutate(clst.merged = forcats::fct_relevel(clst.merged, c("Eryths", "Bcells", "NKcells", "Granulocytes", "Basophils", "Dendritic", "pDendritic", "HSPCs"))) %>%
  arrange(clst.merged)


cells.ordered <- dat.annots.mergeclst$cell
jsub <- jmat[cells.ordered, motifs.keep]


# cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# colvec <- sapply(as.numeric(dat.annots.mergeclst$clst.merged), function(x) cbPalette[[x]])
# colvec <- sapply(as.(dat.annots.mergeclst$cluster), function(x) cbPalette[[x]])
# dat.col <- data.frame(cell = dat.annots.mergeclst$cell, col = colvec, stringsAsFactors = FALSE)

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
                   # ColSideColors = rep("blue", ncol(jsub)), 
                   # ColSideColors = FALSE,
                   RowSideColors = colvec, 
                   # RowSideColors = rep("red", nrow(jsub)), 
                   RowSideLabs = "celltype", 
                   labRow = FALSE, scale = "column", revC = TRUE,
                   distfun = dist, hclustfun = hclust, method = jmeth)

hm.out <- heatmap3(jsub, margins = c(5, 8), cexCol = 0.35, Colv = TRUE, Rowv = NA, 
                   # ColSideColors = rep("blue", ncol(jsub)), 
                   # ColSideColors = FALSE,
                   RowSideColors = colvec, 
                   # RowSideColors = rep("red", nrow(jsub)), 
                   RowSideLabs = "celltype", 
                   labRow = FALSE, scale = "column", revC = FALSE,
                   distfun = dist, hclustfun = hclust, method = jmeth)


hm.out.transpose <- heatmap3(t(jsub), margins = c(5, 8), cexCol = 0.35, Colv = NA, Rowv = TRUE, 
                   # ColSideColors = rep("blue", ncol(jsub)), 
                   # ColSideColors = FALSE,
                   ColSideColors = colvec, 
                   # RowSideColors = rep("red", nrow(jsub)), 
                   ColSideLabs = "celltype", 
                   labCol = FALSE, scale = "row", revC = FALSE, 
                   distfun = dist, hclustfun = hclust, method = jmeth)

hm.out.transpose <- heatmap3(t(jsub), margins = c(5, 8), cexCol = 0.35, Colv = NA, Rowv = TRUE, 
                             # ColSideColors = rep("blue", ncol(jsub)), 
                             # ColSideColors = FALSE,
                             ColSideColors = colvec, 
                             # RowSideColors = rep("red", nrow(jsub)), 
                             ColSideLabs = "celltype", 
                             labCol = FALSE, scale = "row", revC = TRUE,
                             distfun = dist, hclustfun = hclust, method = jmeth)



hm.out.transpose <- heatmap.3(t(jsub), cexCol = 0.35, Colv = NA, Rowv = TRUE, 
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


