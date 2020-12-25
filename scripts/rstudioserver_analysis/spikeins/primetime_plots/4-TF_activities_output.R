# Jake Yeung
# Date of Creation: 2020-11-19
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/4-TF_activities_output.R
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

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


# Load UMAP  ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"
niter <- 1000
# fname <- "jyeung/data/scChiC/glmpca_outputs/glmpca.H3K4me1.bincutoff_0.binskeep_0.byplate.szname_none.reorder_rownames.dupfilt.RData"

indir <- "jyeung/data/scChiC/glmpca_outputs"
# jsuffix <- paste0("glmpca_plate.bincutoff_0.binskeep_0.byplate.szname_none.niter_", niter, ".reorder_rownames.dupfilt")
jsuffix <- paste0("glmpca.H3K4me1.bincutoff_0.binskeep_0.byplate.szname_none.niter_", niter, ".reorder_rownames.dupfilt")
# jsuffix2 <- paste0("glmpca_plate.bincutoff_0.binskeep_0.byplate.szname_none.niter_1000.reorder_rownames.dupfilt.cleanuprows")
jsuffix2 <- paste0("BM_H3K4me1.BM_AllMerged3.glmpca_plate.bincutoff_0.binskeep_0.byplate.szname_none.niter_1000.reorder_rownames.dupfilt.cleanuprows")
# jsuffix2 <- paste0(jsuffix, ".cleanuprows")
fname <- paste0(jsuffix, ".RData")

inf <- file.path(hubprefix, indir, fname)

load(inf, v=T)

dat.umap <- DoUmapAndLouvain(glm.out$factors, jsettings = jsettings)


# Load clusters -----------------------------------------------------------

inf.annot <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins/cell_cluster_table_with_spikeins.H3K4me1.2020-11-18.dupfilt.txt")
dat.annot <- fread(inf.annot)



dat.merge <- left_join(dat.umap, subset(dat.annot, select = c(cell, cluster, batch)))

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2"

outbase <- paste0("motif_activity_H3K4me1.niter_", niter, ".", Sys.Date(), ".uniquecolor")
outname <- paste0(outbase, ".pdf")
outnametxt <- paste0(outbase, ".txt")
outtxt <- file.path(outdir, outnametxt)
outpdf <- file.path(outdir, outname)
pdf(file = outpdf, useDingbats = FALSE)

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
ggplot(dat.merge, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load MARA output --------------------------------------------------------

indir.mara <- file.path(hubprefix, paste0("jyeung/data/scChiC/mara_analysis_BM-AllMerged3_Peaks/H3K4me1/mara_output.lessmemory/", jsuffix2, "-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.H3K4me1/", jsuffix2))
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
colvec <- sapply(as.numeric(dat.annots.mergeclst$clst.merged), function(x) cbPalette[[x]])

dat.col <- data.frame(cell = dat.annots.mergeclst$cell, col = colvec, stringsAsFactors = FALSE)

# plot color

dat.merge.col <- left_join(dat.merge, dat.col)

ggplot(dat.merge.col, aes(x = umap1, y = umap2, color = col)) + 
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
  
  m <- PlotXYWithColor(dat.merge.motifs, xvar = "umap1", yvar = "umap2", cname = jmotif, jtitle = jtitle, cont.color = TRUE) + scale_color_viridis_c()
  print(m)
}


# save output with color name
fwrite(dat.merge.col, file = outtxt, quote = FALSE, sep = "\t", col.names = TRUE)

dev.off()


