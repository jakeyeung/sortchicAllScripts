# Jake Yeung
# Date of Creation: 2020-09-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged_for_MARA/3-make_heatmsps_for_MARA_downstream.H3K27me3.unified_colors.R
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

library(heatmap3)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123



# Load GLMPCA  ------------------------------------------------------------


# load GLMPCA from bins 
jmark <- "H3K27me3"

outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/MARA_output_BM/heatmap_MARA_outputs.", jmark, ".", Sys.Date(), ".unified_colors.pdf")

pdf(outpdf, useDingbats = FALSE)

# jexperi <- "AllMerged"
# mergesize <- "1000"
# nbins <- "1000"
# jcovar.cname <- "ncuts.var.log2.CenteredAndScaled"
# jpenalty <- 1
jsuffix <- "RemoveSmallPeaks"

inf.glmpca <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.Faster.pois.", jsuffix, "KeepBestPlates2.LDAfromPeaks/PZ_", 
                     jmark, 
                     ".AllMerged.KeepBestPlates2.LDAfromPeaks.GLMPCA_var_correction.mergebinsize_100.binskeep_250.covar_cell.var.within.sum.norm.log2.CenteredAndScaled.penalty_1.5.winsorize_FALSE.2020-08-26.RData")
assertthat::assert_that(file.exists(inf.glmpca))

inf.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables/BM_AllMerged.", jmark, ".cell_cluster_table.txt")
mdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged2_Peaks_1000/", jmark, "/mara_output/BM_", jmark, ".BM_AllMerged.VarCorrection.Poisson.GLMPCA.", jsuffix, "-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, "/BM_", jmark, ".BM_AllMerged.VarCorrection.Poisson.GLMPCA.", jsuffix)
assertthat::assert_that(dir.exists(mdir))


load(inf.glmpca, v=T)
dat.annot <- fread(inf.annot)

dat.glmpca.umap <- DoUmapAndLouvain(glm.out$factors, jsettings)

dat.glmpca.umap <- dat.glmpca.umap %>%
  rowwise() %>%
  mutate(plate = ClipLast(as.character(cell), "_"))

dat.glmpca.umap$cond <- sapply(dat.glmpca.umap$cell, GetCondFromSamp, mark = jmark)
dat.glmpca.umap$cond <- factor(dat.glmpca.umap$cond, levels = c("Unenriched", "Linneg", "StemCell"))

dat.glmpca.umap.annot <- left_join(dat.glmpca.umap, subset(dat.annot, select = c(cell, cluster)))

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

m1 <- ggplot(dat.glmpca.umap.annot, aes(x = umap1, y = umap2, color = cluster)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  facet_wrap(~cond)

m2 <- ggplot(dat.glmpca.umap.annot, aes(x = umap1, y = umap2, color = cluster)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
  scale_color_manual(values = cbPalette, na.value = "grey85")

print(m1)
print(m2)

# Summarize in a heatmap --------------------------------------------------


mara.out <- LoadMARA(mdir = mdir, make.cnames = FALSE)

act.mat.clean <- t(subset(mara.out$act.mat, select = -motif))
colnames(act.mat.clean) <- mara.out$act.mat$motif
act.mat.clean.dat <- data.frame(cell = rownames(act.mat.clean), act.mat.clean, stringsAsFactors = FALSE) %>% 
  ungroup() %>%
  mutate(cell = gsub("\\.", "-", cell))

dat.merge <- left_join(dat.glmpca.umap.annot, act.mat.clean.dat)

# Plot output -------------------------------------------------------------

# jmark <- "H3K4me1"

jmotifs <- mara.out$zscores[1:10, ]$motif
# jmotif <- "Erg"

jmotifs <- c("Gata1", "Gfi1", "Ets1", "Runx1", "Rorc", "Bptf", "Egr1", "Otx2", "Zbtb16", "Mbd2", "Mecp2", "Tfdp1", "Erg", "Cbfb", "Plag1", jmotifs)

for (jmotif in jmotifs){
  print(jmotif)
  
  jzscore <- signif(subset(mara.out$zscores, motif == jmotif)$zscore, digits = 2)
  jtitle <- paste(jmotif, "Zscore:", jzscore)
  
  m <- PlotXYWithColor(dat.merge, xvar = "umap1", yvar = "umap2", cname = jmotif, jtitle = jtitle, cont.color = TRUE) + scale_color_viridis_c()
  print(m)
  
}

zscores.sub <- subset(mara.out$zscores, zscore > 0.7)
motifs.keep <- zscores.sub$motif

# show heatmap of output
K <- 6
jmeth <- "ward.D2"

jmat <- as.matrix(subset(act.mat.clean.dat, select = -cell))
rownames(jmat) <- act.mat.clean.dat$cell




# order cells by cluster
ctypes <- unique(dat.annot$cluster)
# collapse some clusters into one
hscs <- list(ctypes[grepl("^HSCs", ctypes)])
names(hscs) <- "HSC"
ilcs <- list(ctypes[grepl("^InnateLyph", ctypes)])
names(ilcs) <- "ILC"
bcells <- list(ctypes[grepl("^Bcell", ctypes)])
names(bcells) <- "Bcells"
granus <- list(ctypes[grepl("^Neutro|^Basophils", ctypes)])
names(granus) <- "Granulocytes"
dends <- list(ctypes[grepl("Dendritic", ctypes)])
names(dends) <- c("DendriticCells")
eryths <- list(ctypes[grepl("Eryth", ctypes)])
names(eryths) <- c("Erythroblasts")

clstr.hash <- hash()
for (lst in list(hscs, ilcs, bcells, granus, dends, eryths)){
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
  filter(cluster != "") %>%
  mutate(clst.merged = ifelse(!is.null(clstr.hash[[cluster]]), clstr.hash[[cluster]], cluster)) %>%
  # umapdist = sqrt(umap1^2 + umap2^2)) %>%
  arrange(clst.merged) %>%
  ungroup() %>%
  mutate(clst.merged = as.factor(clst.merged)) %>%
  filter(!is.na(cluster))



inf.newannot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/GLMPCA_peaks_primetime/H3K4me1_and_H3K27me3_glmpca_peaks_primetime.2020-09-27.", jmark, ".txt")
dat.newannot <- as.data.frame(fread(inf.newannot))

dat.formerge <- subset(dat.newannot, select = c(cell, cluster.renamed)) %>%
  ungroup() %>%
  mutate(cluster.renamed = as.factor(cluster.renamed))

dat.annots.mergeclst <- left_join(dat.annots.mergeclst, dat.formerge)



cbPalette2 <- cbPalette
cbPalette2[[8]] <- cbPalette[[5]]
cbPalette2[[5]] <- cbPalette[[8]]




m1 <- ggplot(dat.annots.mergeclst, aes(x = umap1, y = umap2, color = cluster.renamed)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette2) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m1)

# arrange clusters
jclsts <- c("HSCs", "DC", "Bcells", "InnateLymph", "Neutrophils", "Basophils", "Erythroblasts")

dat.annots.mergeclst$cluster.renamed.fct <- factor(dat.annots.mergeclst$cluster.renamed, levels = jclsts)

dat.annots.mergeclst <- dat.annots.mergeclst %>%
  arrange(cluster.renamed.fct)

cells.ordered <- dat.annots.mergeclst$cell
jsub <- jmat[cells.ordered, motifs.keep]
colvec <- sapply(as.numeric(dat.annots.mergeclst$cluster.renamed), function(x) cbPalette2[[x]])

dat.annots.mergeclst$colvec <- colvec

jcex <- 0.25
# cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
par(mfrow=c(1,1), mar=c(1,1,1,1), mgp=c(3, 1, 0), las=0)
hm.out <- heatmap3(jsub, margins = c(5, 8), Colv = TRUE, Rowv = NA, 
                   cexCol = jcex, cexRow = jcex,
                   # ColSideColors = rep("blue", ncol(jsub)), 
                   RowSideColors = colvec, 
                   # RowSideColors = rep("red", nrow(jsub)), 
                   labRow = FALSE, scale = "column", revC = TRUE,
                   distfun = dist, hclustfun = hclust, method = jmeth, main = jmark)

hm.out <- heatmap3(t(jsub), margins = c(5, 8), ColSideLabs = "celltype", 
                   cexCol = jcex, cexRow = jcex,
                   Colv = NA, Rowv = TRUE, 
                   # ColSideColors = rep("blue", ncol(jsub)), 
                   ColSideColors = colvec, 
                   # RowSideColors = rep("red", nrow(jsub)), 
                   labCol = FALSE, scale = "row", revC = TRUE,
                   distfun = dist, hclustfun = hclust, method = jmeth, main = jmark)

dev.off()
