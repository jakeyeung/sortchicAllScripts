# Jake Yeung
# Date of Creation: 2019-05-26
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/make_figs/make_fig3.R
# Make figure 3


library(dplyr)
library(ggplot2)
library(data.table)
library(heatmap3)
library(ggrepel)
library(topicmodels)
library(hash)
library(JFuncs)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/PlotFunctions.R")

# Constants ---------------------------------------------------------------

kchoose <- 50

jtrajs <- c("granu", "lymphoid", "eryth")

jmark <- "H3K4me1"

# cut tree
K <- 6
jmeth <- "ward.D2"

jscale <- 10^6
jpseudo <- 0

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")  # for TF clusters
cbPalette2 <- c("#999999", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "D3D3D3", "#009E73")  # for Louvain

# significant TFs and label size for dot plot
zscore.cutoff <- 2
labsize <- 5

# Check paths -------------------------------------------------------------

indir.robjs <- paste0("/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/B6_objs")
indir.louvain <- paste0("/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/bamlist_for_merging_build95_B6")
inf <- file.path(indir.louvain, "dat_umap_long_with_louvain.H3K4me1.RData")
assertthat::assert_that(file.exists(inf))
inf.traj <- file.path(indir.robjs, "traj_objs_all_marks.Rdata")
assertthat::assert_that(file.exists(inf.traj))
maramain <- "/hpc/hub_oudenaarden/jyeung/data/scChiC"
maradir <- file.path(maramain, paste0("mara_analysis_cluster_build95_B6_CorrPeakFilt.cells_from_bin_analysis/H3K4me1/mara_output/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-B6_H3K4me1.filt_0.99.center_TRUE_K50-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--K", kchoose, "/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-B6_H3K4me1.filt_0.99.center_TRUE_K", kchoose))
assertthat::assert_that(dir.exists(maradir))
inf.h3k4me1 <- file.path(indir.robjs, "terms_filt_H3K4me1_bin_TRUE_k_50.RData")  # only top 1000 regions considered. Expand this later TODO. Also other bins TODO
assertthat::assert_that(file.exists(inf.h3k4me1))
inf.lda <- file.path(indir.robjs, "LDA_objects_all_marks.Rdata")
assertthat::assert_that(file.exists(inf.lda))

# pdfoutdir <- paste0("/tmp")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1){
  pdfoutdir <- args[[1]]
} else {
  pdfoutdir <- "/tmp"
}
# pdfout <- paste0("~/data/scchic/pdfs/B6_figures/tf_activities/TF_activities.K_", kchoose, ".", Sys.Date(), ".fewer_trajs.pdf")
# pdfout2 <- paste0("~/data/scchic/pdfs/B6_figures/tf_activities/TF_activity_with_gene_exprs.", Sys.Date(), ".fewer_trajs.pdf")

pdfout <- file.path(pdfoutdir, paste0("TF_activities.K_", kchoose, ".", Sys.Date(), ".fewer_trajs.pdf")) 
pdfout2 <- file.path(pdfoutdir, paste0("TF_activity_with_gene_exprs.", Sys.Date(), ".fewer_trajs.pdf"))
# assertthat::assert_that(file.exists(pdfout))
# assertthat::assert_that(file.exists(pdfout2))


# Load data ---------------------------------------------------------------

load(inf, v=T)

# Load trajectories -------------------------------------------------------

load(inf.traj, v=T)

# Load MARA  --------------------------------------------------------------



mara.out <- LoadMARA(maradir, make.cnames = FALSE)
mara.out$act.long$cell <- gsub("\\.", "-", mara.out$act.long$cell)


# Load gene exprs data ----------------------------------------------------

load(inf.h3k4me1, v=T)
load(inf.lda, v=T)

# from bins
out.objs <- out.objs[[jmark]]
add <- subset(out.objs$regions.annotated, SYMBOL == "Cebpb")  # Cebpb is too lowly ranked to be in the terms_filt, add it anyways 
imput.mat <- t(posterior(out.objs$out.lda)$topics %*% posterior(out.objs$out.lda)$terms)
gene.hash <- hash(terms.filt$term, terms.filt$gene)
gene.hash[[add$region_coord]] <- add$SYMBOL
rnames.new <- sapply(rownames(imput.mat), function(x){
  return(ifelse(is.null(gene.hash[[x]]), x, gene.hash[[x]]))
}, USE.NAMES = FALSE)
rownames(imput.mat) <- rnames.new

# Merge and plot ----------------------------------------------------------

dat.merge <- left_join(dat.umap.long, mara.out$act.long)

# jmotif <- "Irf4"
# jmotif <- "Zeb1"
# jmotif <- "Tcf3"
# jmotif <- "Stat2"
# jmotif <- "Yy1"
# jmotif <- "Hoxa5"
# jmotif <- "Rara"
# jmotif <- "Sox6"
# jmotif <- "Zeb1"
# 
# 
# jmotif <- "Cebpb"
# 
# jmotif <- "Spib"
# 
# jmotif <- "Gata3"
# jmotif <- "Zeb1"
# 
# jmotif <- "Tcf3"
# jmotif <- "Foxo3"
# jmotif <- "Meis1"
# 
# jmotif <- "Ebf1"
# jmotif <- "Ikzf2"
# jmotif <- "Spic"



# Zscore ------------------------------------------------------------------


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



# jmeth <- "centroid"
clusters <- hclust(dist(jmat[, cells.keep]), method = jmeth)
# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73")
assertthat::assert_that(length(cbPalette) == K)
clst <- cutree(clusters, k = K)
clst.dat <- data.frame(motif = names(clst), clstr = clst)


pdf(pdfout, useDingbats = FALSE)
# pdf(, useDingbats = FALSE)
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


clusters$labels
myplclust(clusters, clusters$labels, cbPalette[clst], main = jmeth)
jhits <- c("Tal1", "Ets1", "Cebpa", "Foxc1", "Cebpb", "Spic", "Nr4a1", "Ebf3", "Ebf1", "Irf4", "Pax6", "Mafb", "Zeb1", "Meis1", "Bcl6", "Bptf")


pdf(pdfout2, useDingbats = FALSE)
m <- ggplot(left_join(cor.sum %>% dplyr::rename(motif = gene), clst.dat) %>% mutate(motif.lab = ifelse(zscore > zscore.cutoff, motif, NA)),
            aes(x = cor.out, y = range.sqr, size = zscore, label = motif.lab, color = as.character(clstr))) +
  geom_point(alpha = 1) + geom_text_repel(size = labsize) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Pearson Correlation") + ylab("Range of Activity * Range of Expression") +
  scale_color_manual(values = cbPalette, na.value = "lightgrey", name = "Cluster")
print(m)
for (jhit in jhits){
  m <- ggplot(exprs.long %>% filter(gene == jhit & traj %in% jtrajs), aes(x = exprs.log2, y = activity, color = traj)) + 
    geom_point(alpha = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_smooth(method = "lm", se = FALSE, color = "gray50") + 
    scale_color_manual(values = cbPalette2) + ggtitle(jhit) + 
    xlab("log2(scChIC Signal)") + ylab("Activity")
  print(m)
}
dev.off()
