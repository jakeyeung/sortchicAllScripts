# Jake Yeung
# Date of Creation: 2020-03-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/MARA_analysis/make_heatmap_of_MARA_outputs.R
# Load up stuff

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
jmark <- "H3K4me1"
jexperi <- "AllMerged"
mergesize <- "1000"
nbins <- "1000"
jcovar.cname <- "ncuts.var.log2.CenteredAndScaled"
jpenalty <- 1

inf.glm <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.KeepBestPlates2/PZ_", jmark, ".", jexperi, ".KeepBestPlates2.GLMPCA_var_correction.mergebinsize_", mergesize, ".binskeep_", nbins, ".covar_", jcovar.cname, ".penalty_", jpenalty, ".winsorize_TRUE.2020-02-11.RData")
assertthat::assert_that(file.exists(inf.glm))
load(inf.glm, v=T)

inf.annots <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
assertthat::assert_that(file.exists(inf.annots))
load(inf.annots, v=T)

mdir <- paste0("/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/", jmark, "/mara_output/countmat_PZ_fromHiddenDomains_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_250.ntopics_30.2020-02-11-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, "/countmat_PZ_fromHiddenDomains_H3K4me1.AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_250.ntopics_30.2020-02-11")
assertthat::assert_that(dir.exists(mdir))


dat.glmpca.umap <- DoUmapAndLouvain(glm.out$factors, jsettings)

dat.glmpca.umap <- dat.glmpca.umap %>%
  rowwise() %>%
  mutate(plate = ClipLast(as.character(cell), "_"))

dat.glmpca.umap$cond <- sapply(dat.glmpca.umap$cell, GetCondFromSamp, mark = jmark)
dat.glmpca.umap$cond <- factor(dat.glmpca.umap$cond, levels = c("Unenriched", "Linneg", "StemCell"))

dat.glmpca.umap.annot <- left_join(dat.glmpca.umap, subset(dat.umap.glm.fillNAs, select = c(cell, cluster, topic.weight)))

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

jmark <- "H3K4me1"
jmotif <- "Gata1"
jzscore <- signif(subset(mara.out$zscores, motif == jmotif)$zscore, digits = 2)
jtitle <- paste(jmotif, "Zscore:", jzscore)

m <- PlotXYWithColor(dat.merge, xvar = "umap1", yvar = "umap2", cname = jmotif, jtitle = jtitle, cont.color = TRUE) + scale_color_viridis_c()
print(m)

zscores.sub <- subset(mara.out$zscores, zscore > 0.7)
motifs.keep <- zscores.sub$motif

# show heatmap of output
K <- 6
jmeth <- "ward.D2"

jmat <- as.matrix(subset(act.mat.clean.dat, select = -cell))
rownames(jmat) <- act.mat.clean.dat$cell




# order cells by cluster
ctypes <- unique(dat.umap.glm.fillNAs$cluster)
# collapse some clusters into one
hscs <- list(ctypes[grepl("^HSCs", ctypes)])
names(hscs) <- "HSC"
ilcs <- list(ctypes[grepl("^ILC", ctypes)])
names(ilcs) <- "ILC"
bcells <- list(ctypes[grepl("^Bcell", ctypes)])
names(bcells) <- "Bcells"
granus <- list(ctypes[grepl("^Neutro|^Basophils", ctypes)])
names(granus) <- "Granulocytes"
dends <- list(ctypes[grepl("Dendritic", ctypes)])
names(dends) <- c("DendriticCells")

clstr.hash <- hash()
for (lst in list(hscs, ilcs, bcells, granus, dends)){
  # print(lst)
  jname <- names(lst)[[1]]
  # print(jname)
  # print(lst[[1]])
  for (ctype in lst[[1]]){
    print(paste(ctype, jname))
    clstr.hash[[ctype]] <- jname
  }
}

dat.annots.mergeclst <- dat.umap.glm.fillNAs %>%
  rowwise() %>%
  mutate(clst.merged = ifelse(!is.null(clstr.hash[[cluster]]), clstr.hash[[cluster]], cluster),
         umapdist = sqrt(umap1^2 + umap2^2)) %>%
  arrange(clst.merged, umapdist) %>%
  ungroup() %>%
  mutate(clst.merged = as.factor(clst.merged)) %>%
  filter(!is.na(cluster))
cells.ordered <- dat.annots.mergeclst$cell
jsub <- jmat[cells.ordered, motifs.keep]



cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
colvec <- sapply(as.numeric(dat.annots.mergeclst$clst.merged), function(x) cbPalette[[x]])
par(mfrow=c(1,1), mar=c(1,1,1,1), mgp=c(3, 1, 0), las=0)
hm.out <- heatmap3(jsub, margins = c(5, 8), cexCol = 0.35, Colv = TRUE, Rowv = NA, 
                   ColSideColors = rep("blue", ncol(jsub)), 
                   # RowSideColors = colvec, 
                   RowSideColors = rep("red", nrow(jsub)), 
                   labRow = FALSE, scale = "column", revC = TRUE,
                   distfun = dist, hclustfun = hclust, method = jmeth)


# save to output try to run plots on macbook
save.image(file = "/home/jyeung/hpc/scChiC/from_rstudioserver/rsessions/make_TF_heatmap.2020-03-04.pdf")
