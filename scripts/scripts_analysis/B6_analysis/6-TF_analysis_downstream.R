# Jake Yeung
# Date of Creation: 2019-05-13
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/6-TF_analysis_downstream.R
# Plot activities in umap

library(dplyr)
library(ggplot2)
library(data.table)
library(heatmap3)
library(ggrepel)
library(topicmodels)
library(hash)
library(JFuncs)

# source("~/data/scchic/m")
# source("/Users/yeung/data/scchic/from_cluster/mara_outputs_B6/mara_analysis_cluster_build95_B6_CorrPeakFilt.cells_from_bin_analysis")
source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/PlotFunctions.R")

# Load data ---------------------------------------------------------------

inf <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6/dat_umap_long_with_louvain.H3K4me1.RData"

assertthat::assert_that(file.exists(inf))

load(inf, v=T)


# Load trajectories -------------------------------------------------------

inf.traj <- paste0("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata")
load(inf.traj, v=T)

# Load MARA  --------------------------------------------------------------

kchoose <- 50
maradir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_outputs_B6/mara_analysis_cluster_build95_B6_CorrPeakFilt.cells_from_bin_analysis/H3K4me1/mara_output/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-B6_H3K4me1.filt_0.99.center_TRUE_K50-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--K", kchoose, "/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-B6_H3K4me1.filt_0.99.center_TRUE_K", kchoose)
assertthat::assert_that(dir.exists(maradir))

mara.out <- LoadMARA(maradir, make.cnames = FALSE)
mara.out$act.long$cell <- gsub("\\.", "-", mara.out$act.long$cell)


# Merge and plot ----------------------------------------------------------

dat.merge <- left_join(dat.umap.long, mara.out$act.long)

jmotif <- "Irf4"
jmotif <- "Zeb1"
jmotif <- "Tcf3"
jmotif <- "Stat2"
jmotif <- "Yy1"
jmotif <- "Hoxa5"
jmotif <- "Rara"
jmotif <- "Sox6"
jmotif <- "Zeb1"


jmotif <- "Cebpb"

jmotif <- "Spib"

jmotif <- "Gata3"
jmotif <- "Zeb1"

jmotif <- "Tcf3"
jmotif <- "Foxo3"
jmotif <- "Meis1"

jmotif <- "Ebf1"
jmotif <- "Ikzf2"
jmotif <- "Spic"



# Zscore ------------------------------------------------------------------

jtrajs <- c("granu", "lymphoid", "eryth", "mega")

jmark <- "H3K4me1"
traj.annot <- lapply(jtrajs, function(jtraj) data.frame(cell = trajs[[jmark]][[jtraj]]$cell, traj = jtraj, lambda = trajs[[jmark]][[jtraj]]$lambda)) %>%
  bind_rows()

# heatmap
# reorder columns by lambda
cell.ordering.dat <- traj.annot %>%
  filter(grepl(jmark, cell)) %>%
  group_by(traj) %>%
  arrange(traj, lambda)
cells <- cell.ordering.dat$cell

motifs.keep <- subset(mara.out$zscore, zscore >= 2)$motif
jmat <- subset(mara.out$act.mat, motif %in% motifs.keep)
colnames(jmat) <- gsub("\\.", "-", colnames(jmat))  # B6.13W1.BM.H3K4me1.3_378 -> B6-13W1-BM-H3K4me1-3_378
jmotifs <- jmat$motif
jmat$motif <- NULL

jmat <- as.matrix(t(scale(t(jmat), center=TRUE, scale=TRUE)))
# jmat <- as.matrix(scale(jmat, center=FALSE, scale=FALSE))
rownames(jmat) <- jmotifs

cells.keep <- cells[which(cells %in% colnames(jmat))]


# cut tree
K <- 6
jmeth <- "ward.D2"
# jmeth <- "centroid"
clusters <- hclust(dist(jmat[, cells.keep]), method = jmeth)
# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73")
assertthat::assert_that(length(cbPalette) == K)
clst <- cutree(clusters, k = K)
clst.dat <- data.frame(motif = names(clst), clstr = clst)

pdf(paste0("~/data/scchic/pdfs/B6_figures/tf_activities/TF_activities.K_", kchoose, ".", Sys.Date(), ".pdf"), useDingbats = FALSE)
  hm.out <- heatmap3(t(jmat[, cells.keep]), margins = c(5, 8), cexCol = 0.35, Colv = TRUE, Rowv = NA, 
                     labRow = FALSE, scale = "column", revC = TRUE,
                     distfun = dist, hclustfun = hclust, method = jmeth)
  
  hm.out.t <- heatmap3(jmat[, cells.keep], margins = c(5, 8), cexRow = 0.35, Colv = NA, Rowv = TRUE, 
                     labCol = FALSE, scale = "row", revC = TRUE, 
                     distfun = dist, hclustfun = hclust, method = jmeth)
  myplclust(clusters, clusters$labels, cbPalette[clst], main = jmeth)
  
  jmotifs <- c("Cebpb", "Spic", "Tal1", "Ebf1", "Ebf3", "Hmbox1", "Chd1", "Cdx1", "Gata3", "Bptf")
  for (jmotif in jmotifs){
    m <- PlotXYWithColor(dat.merge %>% filter(motif == jmotif), xvar = "umap1", yvar = "umap2", cname = "activity", jtitle = jmotif, jcol = scales::muted('darkred'))
    print(m)
  }
dev.off()

# Correlate with expression  ----------------------------------------------

jmotifs <- mara.out$zscores$motif
jscale <- 10^6
jpseudo <- 0

from.tss <- FALSE

if (from.tss){
  # from TSS
  inf.tss <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs_TSS_B6/lda_out_meanfilt.B6-BM-", jmark, ".CountThres0.K-50.Robj")
  load(inf.tss, v=T)
  # imput mat from TSS
  imput.mat <- t(posterior(out.lda[[1]])$topics %*% posterior(out.lda[[1]])$terms)
  # rename colnames
  prefix <- paste0("B6-13W1-BM-", jmark)
  cnames.new <- names(colnames(imput.mat))
  repl <- sapply(cnames.new, function(x) strsplit(x, "\\.")[[1]][[14]], USE.NAMES = FALSE)
  cellindx <- sapply(cnames.new, function(x) strsplit(strsplit(x, "_")[[1]][[5]], "\\.")[[1]][[1]], USE.NAMES = FALSE)
  cnames.new2 <- paste0(prefix, "-", repl, "_", cellindx)
  colnames(imput.mat) <- cnames.new2
  rownames(imput.mat) <- sapply(rownames(imput.mat), function(x) strsplit(x, ";")[[1]][[2]])
} else {
  # from bins
  load("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata", v=T)
  load("/Users/yeung/data/scchic/robjs/B6_objs/terms_filt_H3K4me1_bin_TRUE_k_50.RData", v=T)  # terms.filt
  # inf.bin <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_bins_B6/ldaAnalysisBins_BinCellFilt2/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.TRUE.no_filt/lda_out_meanfilt.B6_H3K4me1_pcutoff_0.CountThres0.K-25_30_40_50.Robj"
  # kchoose <- 50
  # out.objs <- LoadLDABins(jmark, jbin = NA, top.thres = 0.995, inf = inf.bin, convert.chr20.21.to.X.Y = FALSE, add.chr.prefix = TRUE, choose.k = kchoose)
  load("/Users/yeung/data/scchic/robjs/B6_objs/LDA_objects_all_marks.Rdata", v=T)
  out.objs <- out.objs[[jmark]]
  
  add <- subset(out.objs$regions.annotated, SYMBOL == "Cebpb")
  imput.mat <- t(posterior(out.objs$out.lda)$topics %*% posterior(out.objs$out.lda)$terms)
  gene.hash <- hash(terms.filt$term, terms.filt$gene)
  gene.hash[[add$region_coord]] <- add$SYMBOL
  rnames.new <- sapply(rownames(imput.mat), function(x){
    return(ifelse(is.null(gene.hash[[x]]), x, gene.hash[[x]]))
  }, USE.NAMES = FALSE)
  rownames(imput.mat) <- rnames.new
}


rows.keep <- which(rownames(imput.mat) %in% jmotifs)

exprs.long <- data.frame(gene = rownames(imput.mat[rows.keep, ]), imput.mat[rows.keep, ]) %>%
  tidyr::gather(key = cell, value = exprs, -gene) %>%
  mutate(exprs.log2 = log2(exprs * jscale + jpseudo))
# change . to - before merging
exprs.long$cell <- gsub("\\.", "-", exprs.long$cell)
# add umap
exprs.long <- left_join(exprs.long, dat.umap.long.trajs[[jmark]] %>% dplyr::select(umap1, umap2, cell))
# add activity
exprs.long <- left_join(exprs.long, mara.out$act.long %>% dplyr::rename(gene = motif))
# add trajectory

ctypes <- c("eryth", "granu", "lymphoid", "mega")
traj.annot <- lapply(ctypes, function(jtraj) data.frame(cell = trajs[[jmark]][[jtraj]]$cell, traj = jtraj, lambda = trajs[[jmark]][[jtraj]]$lambda)) %>%
  bind_rows()
# traj.annot <- left_join(traj.annot, dat.umap.long.trajs[[jmark]] %>% dplyr::select(cell))

exprs.long <- left_join(exprs.long, traj.annot)

m1 <- PlotXYWithColor(exprs.long %>% filter(gene == "Tal1"), xvar = "umap1", yvar = "umap2", cname = "exprs.log2")
m2 <- PlotXYWithColor(exprs.long %>% filter(gene == "Tal1" & !is.na(activity)), xvar = "umap1", yvar = "umap2", cname = "activity", jcol = scales::muted("darkred"))

multiplot(m1, m2, cols = 2)

# cbPalette2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "D3D3D3")
cbPalette2 <- c("#999999", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "D3D3D3", "#009E73")



# Do genome-wide correlations and z-score sizes ---------------------------

cor.sum <- exprs.long %>% 
  filter(!is.na(activity)) %>%
  group_by(gene) %>%
  summarise(cor.out = cor(x = exprs.log2, y = activity),
            range.sqr = abs(diff(range(activity))) * abs(diff(range(exprs.log2))),
            exprs.mean = mean(exprs.log2),
            act.mean = mean(activity),
            range.act = abs(diff(range(activity))),
            range.exprs = abs(diff(range(exprs.log2)))) %>%
  arrange(desc(abs(cor.out)))
cor.sum <- left_join(cor.sum, mara.out$zscores %>% dplyr::rename(gene = motif))

cor.sum <- cor.sum %>%
  rowwise() %>%
  mutate(motif.lab = ifelse(zscore > 2, gene, NA))

# plot by zscore first
m <- ggplot(cor.sum, aes(x = cor.out, y = zscore, label = motif.lab)) + geom_point()  + geom_text_repel() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m)

cor.sum.bytraj <- exprs.long %>%
  filter(!is.na(activity) & !is.na(traj)) %>%
  group_by(gene, traj) %>%
  summarise(cor.out = cor(x = exprs.log2, y = activity),
            range.sqr = abs(diff(range(activity))) * abs(diff(range(exprs.log2))),
            exprs.mean = mean(exprs.log2),
            act.mean = mean(activity),
            range.act = abs(diff(range(activity))),
            range.exprs = abs(diff(range(exprs.log2)))) %>%
  arrange(desc(abs(cor.out)))
cor.sum.bytraj <- left_join(cor.sum.bytraj, mara.out$zscores %>% dplyr::rename(gene = motif)) %>%
  rowwise() %>%
  mutate(motif.lab = ifelse(zscore > 2, gene, NA))

cor.sum.bytraj.filt <- cor.sum.bytraj %>%
  group_by(gene) %>%
  filter(abs(cor.out) == max(abs(cor.out)))

m <- ggplot(cor.sum.bytraj.filt, aes(x = cor.out, y = zscore, label = motif.lab, color = traj)) + geom_point()  + geom_text_repel() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)
print(m)

# Range sqr  --------------------------------------------------------------

# color by cluster 

zscore.cutoff <- 2
labsize <- 5
clusters$labels
myplclust(clusters, clusters$labels, cbPalette[clst], main = jmeth)
jhits <- c("Tal1", "Ets1", "Cebpa", "Foxc1", "Cebpb", "Spic", "Nr4a1", "Ebf3", "Ebf1", "Irf4", "Pax6", "Mafb", "Zeb1", "Meis1", "Bcl6", "Bptf")


pdf(paste0("~/data/scchic/pdfs/B6_figures/tf_activities/TF_activity_with_gene_exprs.", Sys.Date(), ".pdf"), useDingbats = FALSE)
m <- ggplot(left_join(cor.sum %>% dplyr::rename(motif = gene), clst.dat) %>% mutate(motif.lab = ifelse(zscore > zscore.cutoff, motif, NA)),
            aes(x = cor.out, y = range.sqr, size = zscore, label = motif.lab, color = as.character(clstr))) +
  geom_point(alpha = 1) + geom_text_repel(size = labsize) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Pearson Correlation") + ylab("Range of Activity * Range of Expression") +
  scale_color_manual(values = cbPalette, na.value = "lightgrey", name = "Cluster")
print(m)
for (jhit in jhits){
  m <- ggplot(exprs.long %>% filter(gene == jhit & !is.na(traj)), aes(x = exprs.log2, y = activity, color = traj)) + 
    geom_point(alpha = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_smooth(method = "lm", se = FALSE, color = "gray50") + 
    scale_color_manual(values = cbPalette2) + ggtitle(jhit) + 
    xlab("log2(scChIC Signal)") + ylab("Activity")
  print(m)
}
dev.off()

# m <- ggplot(cor.sum, aes(x = cor.out, y = range.sqr, label = motif.lab, size = zscore)) + geom_point()  + geom_text_repel(size = 5) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# print(m)

# 
# ggplot(exprs.long %>% filter(gene == "Ets1" & !is.na(traj)), aes(x = exprs.log2, y = activity, color = traj)) + 
#   geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   geom_smooth(method = "lm", se = FALSE, color = "gray50") + 
#   scale_color_manual(values = cbPalette2)
# 
# ggplot(exprs.long %>% filter(gene == "Cebpa" & !is.na(traj)), aes(x = exprs.log2, y = activity, color = traj)) + 
#   geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   geom_smooth(method = "lm", se = FALSE, color = "gray50") + 
#   scale_color_manual(values = cbPalette2)
# 
# ggplot(exprs.long %>% filter(gene == "Foxc1" & !is.na(traj)), aes(x = exprs.log2, y = activity, color = traj)) + 
#   geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   geom_smooth(method = "lm", se = FALSE, color = "gray50") + 
#   scale_color_manual(values = cbPalette2)
# 
# ggplot(exprs.long %>% filter(gene == "Cebpb" & !is.na(traj)), aes(x = exprs.log2, y = activity, color = traj)) + 
#   geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   geom_smooth(method = "lm", se = FALSE, color = "gray50") + 
#   scale_color_manual(values = cbPalette2)
# 
# ggplot(exprs.long %>% filter(gene == "Spic" & !is.na(traj)), aes(x = exprs.log2, y = activity, color = traj)) + 
#   geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   geom_smooth(method = "lm", se = FALSE, color = "gray50") + 
#   scale_color_manual(values = cbPalette2)
# 
# ggplot(exprs.long %>% filter(gene == "Nr4a1" & !is.na(traj)), aes(x = exprs.log2, y = activity, color = traj)) + 
#   geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   geom_smooth(method = "lm", se = FALSE, color = "gray50") + 
#   scale_color_manual(values = cbPalette2)

