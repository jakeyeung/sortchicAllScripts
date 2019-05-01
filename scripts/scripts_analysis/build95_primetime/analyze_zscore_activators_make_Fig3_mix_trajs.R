# Jake Yeung
# Date of Creation: 2019-04-23
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/analyze_zscore_activators_make_Fig3_mix_trajs.R
# Mix trajs

rm(list=ls())

library(dplyr)
library(ggplot2)
library(ggrepel)
library(hash)
library(data.table)
library(JFuncs)

# plot heatmap
library(heatmap3)
library(gplots)
library(made4)


source("scripts/Rfunctions/GetMetaData.R")
source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/MatchCellNameToSample.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/TrajFunctions.R")
source("scripts/Rfunctions/AnnotationFunctions.R")

myplclust <- function(hclust, lab = hclust$labels, lab.col = rep(1, length(hclust$labels)), 
                      hang = 0.1, ...) {
  ## modifiction of plclust for plotting hclust objects *in colour*!  Copyright
  ## Eva KF Chan 2009 Arguments: hclust: hclust object lab: a character vector
  ## of labels of the leaves of the tree lab.col: colour for the labels;
  ## NA=default device foreground colour hang: as in hclust & plclust Side
  ## effect: A display of hierarchical cluster with coloured leaf labels.
  y <- rep(hclust$height, 2)
  x <- as.numeric(hclust$merge)
  y <- y[which(x < 0)]
  x <- x[which(x < 0)]
  x <- abs(x)
  y <- y[order(x)]
  x <- x[order(x)]
  plot(hclust, labels = FALSE, hang = hang, ...)
  text(x = x, y = y[hclust$order] - (max(hclust$height) * hang), labels = lab[hclust$order], 
       col = lab.col[hclust$order], srt = 90, adj = c(1, 0.5), xpd = NA, ...)
}


# Constants ---------------------------------------------------------------

zscore.cutoff <- 2
# jcols <- c("#A29F9D", "#E7A100", "#6AB7E6")
jcols <- GetTrajColors(as.hash = FALSE)
jscale.fac <- 10^6
jpseudo <- 0
jtrajs <- c("eryth", "granu", "lymphoid")

inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"
load(inf, v=T)

traj.mix.lst <- GetTrajMixed()

trajs <- traj.mix.lst$trajs.mixed
dat.umap.long.trajs <- traj.mix.lst$dat.umap.mixed

# need X1, X2, mark, cell, louvain)
dat.umap.long.trajs$H3K4me1 <- dat.umap.long.trajs$H3K4me1 %>%
  dplyr::rename(X1 = umap1, X2 = umap2) %>%
  dplyr::select(X1, X2, mark, cell, louvain)
dat.umap.long.trajs$H3K4me3 <- dat.umap.long.trajs$H3K4me3 %>%
  dplyr::select(X1, X2, mark, cell, louvain)
dat.umap.long.trajs$H3K27me3 <- dat.umap.long.trajs$H3K27me3 %>%
  dplyr::select(X1, X2, mark, cell, louvain)
dat.umap.long.trajs$H3K9me3 <- dat.umap.long.trajs$H3K9me3 %>%
  dplyr::select(X1, X2, mark, cell, louvain)

dat.umap.long.trajs <- dat.umap.long.trajs %>%
  bind_rows()

# 
# dat.umap.long.trajs <- dat.umap.long.trajs %>%
#   bind_rows() %>%
#   dplyr::rename(X1 = umap1, X2 = umap2)
# jtrajs.all <- names(trajs[[jmark]])

jmark <- "H3K4me1"
traj.annot <- lapply(jtrajs, function(jtraj) data.frame(cell = trajs[[jmark]][[jtraj]]$cell, traj = jtraj, lambda = trajs[[jmark]][[jtraj]]$lambda)) %>%
  bind_rows()
# 
# traj.renamed <- SwapTrajName(traj.annot$traj)
# 
# traj.annot$traj <- traj.renamed
# 
# traj.annot <- subset(traj.annot, traj %in% jtrajs)

# Load zscores ------------------------------------------------------------

experis <- read.table("data/outputs/experinames.txt", stringsAsFactors = FALSE)
experihash <- hash(sapply(experis$V1, function(x) strsplit(x, "_")[[1]][[1]]), sapply(experis$V1, function(x) paste(strsplit(x, "_")[[1]][2:3], collapse="_")))
switch.rep.hash <- GetRepSwitchHash(experihash)
# add m1_S1 which is a bug from m1_S10
# fix a bug with GetTechRep which causes S1 to show up 
switch.rep.hash[["BM_H3K4me3_m1_S1"]] <- "BM_H3K4me3_m1_rep1"


indir.mara1 <- "/Users/yeung/data/scchic/from_cluster/mara_analysis_cluster_build95_CorrPeakFilt.withchr.cells_from_bin_analysis/H3K4me1/mara_output/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_H3K4me1.filt_0.99.center_TRUE_K50-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--K50/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_H3K4me1.filt_0.99.center_TRUE_K50"
indir.mara2 <- "/Users/yeung/data/scchic/from_cluster/mara_analysis_cluster_build95_CorrPeakFilt.withchr.cells_from_bin_analysis/H3K4me3/mara_output/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_H3K4me3.filt_0.99.center_TRUE_K50-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--K50/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_H3K4me3.filt_0.99.center_TRUE_K50"

assertthat::assert_that(dir.exists(indir.mara1))
assertthat::assert_that(dir.exists(indir.mara2))

mara.out1 <- LoadMARA(indir.mara1, rep.prefix = "", swap.tech.rep = switch.rep.hash, filesuffix = "")
mara.out2 <- LoadMARA(indir.mara2, rep.prefix = "", swap.tech.rep = switch.rep.hash, filesuffix = "")

act.long <- bind_rows(mara.out1$act.long, mara.out2$act.long)

# Plot zscores ------------------------------------------------------------

zscores.merged <- inner_join(mara.out1$zscores, mara.out2$zscores, by = "motif") %>%
  mutate(motif.lab = ifelse(zscore.x > 3 | zscore.y > 3, motif, NA))

m <- ggplot(zscores.merged, aes(x = zscore.x, zscore.y, label = motif.lab)) + geom_point() + geom_text_repel() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Zscore H3K4me1") + ylab("Zscore H3K4me3") 
print(m)


# Plot activity -----------------------------------------------------------

# check number of cells
mreps <- dat.umap.long.trajs %>% filter(mark == "H3K4me3")
mreps.uniq <- unique(paste(sapply(mreps$cell, function(x) strsplit(x, "_")[[1]][[3]]), sapply(mreps$cell, function(x) strsplit(x, "_")[[1]][[4]]), sep = "_"))

act.uniq <- unique(paste(sapply(mara.out2$act.long$cell, function(x) strsplit(x, "_")[[1]][[3]]), sapply(mara.out2$act.long$cell, function(x) strsplit(x, "_")[[1]][[4]]), sep = "_"))

jmarks <- c("H3K4me1", "H3K4me3")
dat.sub <- dat.umap.long.trajs %>% filter(mark %in% jmarks)

act.exprs.umap <- left_join(act.long, dat.sub %>% dplyr::select(X1, X2, mark, cell, louvain))

jmark <- "H3K4me3"
jmotif <- "Cebpb"
jmotif <- "Mecp2"

jmotif <- "Nrf1"
jmotif <- "Tal1"

jmotif <- "Runx1"
jmotif <- "Spic"
jmotif <- "Spi1"
jmotif <- "Cebpb"
jmotif <- "Cebpa"
jmotif <- "Etv4"
jmotif <- "Spib"
jmotif <- "Pax6"
jmotif <- "Gata3"
mlst <- lapply(jmarks, function(jmark){
  m1 <- PlotXYWithColor(act.exprs.umap %>% filter(motif == jmotif & mark == jmark), xvar = "X1", yvar = "X2", cname = "activity", jcol = scales::muted("darkred")) + ggtitle(paste(jmark, jmotif))
  return(m1)
})
multiplot(mlst[[1]], mlst[[2]], cols = 2)


# Do clustering of cell distributions for each mark -----------------------

jmark <- "H3K4me1"
(motifs.keep <- subset(zscores.merged, zscore.x > zscore.cutoff)$motif)
jmat <- subset(mara.out1$act.mat, motif %in% motifs.keep)
jmotifs <- jmat$motif
jmat$motif <- NULL
jmat <- as.matrix(scale(jmat, center=FALSE, scale=TRUE))
rownames(jmat) <- jmotifs

# clusters <- hclust(dist(scale(jmat)), method = "ward.D")
# jmeth <- "ward.D"
jmeth <- "ward.D2"

# clusters <- hclust(dist(scale(jmat)), method = jmeth)
clusters <- hclust(dist(jmat[, cells.keep]), method = jmeth)

# heatmap

# reorder columns by lambda 
cell.ordering.dat <- traj.annot %>%
  filter(grepl(jmark, cell)) %>%
  group_by(traj) %>%
  arrange(traj, lambda)
cells <- cell.ordering.dat$cell

cells.keep <- cells[which(cells %in% colnames(jmat))]

K <- 6
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0")
assertthat::assert_that(length(cbPalette) == K)
clst <- cutree(clusters, k = K)
clst.dat <- data.frame(motif = names(clst), clstr = clst)

pdf(paste0("~/data/scchic/pdfs/heatmap_global_analysis_motifs_H3K4me1.", Sys.Date(), ".ClstMeth.", jmeth, ".pdf"), useDingbats = FALSE)
  hm.out <- heatmap3(t(jmat[, cells.keep]), margins = c(5, 8), cexRow = 0.25, Colv = TRUE, Rowv = NA,
                     labRow = FALSE, scale = "column", revC = TRUE, 
                     distfun = dist, hclustfun = hclust, method = jmeth)
  myplclust(clusters, clusters$labels, cbPalette[clst])
dev.off()

# Add H3K4me1 expression --------------------------------------------------

tssdist <- 50000; jdate <- "2019-04-22"; Kvec <- "50"
inf.lda <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs_peaks_exprs_merge_clusters_explore/lda_out_meanfilt.PZ-BM-H3K4me1.CountThres0.K-", Kvec, "_GeneTSS.Dedup.", jdate, ".", tssdist, ".Robj")

assertthat::assert_that(file.exists(inf.lda))
load(inf.lda, v=T)

out.lda <- out.lda[[length(out.lda)]]
out.lda@documents <- SwitchColnames(unname(out.lda@documents), jsplit = "-")
tm.result <- posterior(out.lda)

mat.impute <- t(tm.result$topics %*% tm.result$terms)

genes <- sapply(rownames(mat.impute), function(x){
  g <- tryCatch({
    return(strsplit(x, ";")[[1]][[2]])
  }, error = function(e){
    return("Peak")
  })
}, USE.NAMES = FALSE)

genes.keep <- mara.out1$zscores$motif

genes.keep.i <- which(genes %in% genes.keep)

mat.impute.sub <- mat.impute[genes.keep.i, ]

exprs.long <- data.frame(peak = rownames(mat.impute.sub), gene = genes[genes.keep.i], as.data.frame(mat.impute.sub)) %>%
  tidyr::gather(key = "cell", value = "exprs", c(-peak, -gene))


# Plot 3 major hits -------------------------------------------------------

# test
act.exprs.umap2 <- inner_join(act.exprs.umap, exprs.long %>% dplyr::rename(motif = gene))
jsize <- 1

# jmotifs <- c("Cebpb", "Tal1", "Ebf1")
# jpaths <- c("granu", "eryth", "lymphoid")

# jmotifs <- motifs.
act.exprs.umap2 <- left_join(act.exprs.umap2, traj.annot)
act.exprs.umap2$exprs.log <- log2(act.exprs.umap2$exprs * jscale.fac + jpseudo)

pdf(paste0("~/data/scchic/pdfs/TF_analysis_activity_umap_Figure3MixedMixed.", Sys.Date(), ".pdf"), useDingbats = FALSE)
for (i in seq(length(motifs.keep))){
  jmotif <- motifs.keep[[i]]
  print(jmotif)
  # jpath <- jpaths[[i]]
  # infer path from motif
  act.highest <- act.exprs.umap2 %>% filter(motif == jmotif) %>%
    filter(!is.na(traj)) %>%
    group_by(traj) %>%
    summarise(activity = mean(activity)) %>%
    filter(activity == max(activity))
  if (nrow(act.highest) == 0){
    print(paste("Skip", jmotif))
    next
  }
  jpath <- act.highest$traj
  assertthat::assert_that(length(jpath) == 1)
  print(jpath)
  m1 <- PlotXYWithColor(act.exprs.umap2 %>% filter(motif == jmotif), xvar = "X1", yvar = "X2", cname = "activity", jcol = scales::muted("darkred")) + ggtitle(jmotif) + 
    # geom_path(data = trajs.spring[[jmark]][[jpath]], inherit.aes = FALSE, mapping = aes(x = umap1, y = umap2), size = jsize, color = "black")
    geom_path(data = trajs[[jmark]][[jpath]], inherit.aes = FALSE, mapping = aes(x = umap1, y = umap2), size = jsize, color = "black")
  print(m1)
}
dev.off()

# Show correlation with pseudotime  ---------------------------------------



# show positive correlations and negative ones 

jmotif <- "Foxc1"
jmotif <- "Tal1"
jmotif <- "Cebpb"
jmotif <- "Erbf1"
jmotif <- "Spib"

# plot globa
m.global <- ggplot(mara.out1$zscores %>% mutate(motif.lab = ifelse(zscore > zscore.cutoff, motif, NA), motif = factor(motif, levels = motif)), 
       aes(x = motif, y = zscore, label = motif.lab)) + geom_point(alpha = 0.2, size = 0.2)  +  geom_text_repel() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  xlab("")

jmotifs <- motifs.keep
pdf(paste0("~/data/scchic/pdfs/TF_analysis_gene_correlations_pseudotime_Figure3Mixed.", Sys.Date(), ".pdf"), useDingbats = FALSE)
  print(m.global)
for (jmotif in jmotifs){
  print(jmotif)
  merged.sub <- subset(act.exprs.umap2, motif == jmotif) %>%
    tidyr::gather(., key = "type", value = "exprsORactivity", activity, exprs.log) %>%
    filter(!is.na(traj)) %>%
    group_by(type, motif) %>%
    mutate(exprsORactivity = scale(exprsORactivity, center = TRUE, scale = TRUE))
  if (nrow(merged.sub) == 0){
    print(paste("Empty dataframe, skipping", jmotif))
    next
  }
  merged.sub.wide <- subset(act.exprs.umap2, motif == jmotif) %>% filter(!is.na(traj))
  m.pseudo <- ggplot(merged.sub, aes(x = lambda, y = exprsORactivity, group = type, color = type)) + geom_point(alpha = 0.2) + facet_wrap(~traj) + geom_smooth(method = "lm", se = FALSE) +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(jmotif) + 
    xlab("Pseudotime") + ylab("Activity or scChiC signal (scaled)")
  m.pseudo.xy <- ggplot(merged.sub.wide, aes(y = activity, x = exprs.log, color = traj)) + geom_point(alpha = 0.5) + 
    geom_smooth(method = "lm", se = FALSE, mapping = aes(x = exprs.log, y = activity), inherit.aes = FALSE, color = "black") +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    ggtitle(jmotif) + 
    xlab("log2(scChiC signal)") + ylab("Activity") + scale_color_manual(values = jcols)
  print(m.pseudo)
  print(m.pseudo.xy)
}
dev.off()


# Do global correlations --------------------------------------------------


# maybe correlation versus rectangle 
cor.sum <- act.exprs.umap2 %>%
  group_by(motif) %>%
  summarise(cor.out = cor(x = exprs.log, y = activity),
            range.sqr = abs(diff(range(activity))) * abs(diff(range(exprs.log))),
            exprs.mean = mean(exprs.log),
            act.mean = mean(activity),
            range.act = abs(diff(range(activity))),
            range.exprs = abs(diff(range(exprs.log))))
cor.sum <- left_join(cor.sum, mara.out1$zscores)


labsize <- 5
pdf(paste0("~/data/scchic/pdfs/TF_analysis_correlation_global_Figure3.", Sys.Date(), ".pdf"), useDingbats = FALSE)
  m <- ggplot(cor.sum %>% mutate(motif.lab = ifelse(zscore > zscore.cutoff, motif, NA)), 
              aes(x = cor.out, y = range.sqr, size = zscore, label = motif.lab)) + 
    geom_point(alpha = 0.3) + geom_text_repel(size = labsize) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    xlab("Pearson Correlation") + ylab("Range of Activity * Range of Expression")
  print(m)
  
  m <- ggplot(cor.sum %>% mutate(motif.lab = ifelse(zscore > zscore.cutoff, motif, NA)), 
              aes(x = cor.out, y = range.act, size = zscore, label = motif.lab)) + 
    geom_point(alpha = 0.3) + geom_text_repel(size = labsize) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  print(m)
  
  m <- ggplot(cor.sum %>% mutate(motif.lab = ifelse(zscore > zscore.cutoff, motif, NA)), 
              aes(x = cor.out, y = range.exprs, size = zscore, label = motif.lab)) + 
    geom_point(alpha = 0.3) + geom_text_repel(size = labsize) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  print(m)
  
  m <- ggplot(cor.sum %>% mutate(motif.lab = ifelse(zscore > zscore.cutoff, motif, NA)), 
              aes(x = cor.out, y = zscore, label = motif.lab)) + 
    geom_point(alpha = 0.3) + geom_text_repel(size = labsize) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  print(m)
  
dev.off()


# Do further downstream of th edendrogram  --------------------------------





# for each clstr what are the correlations?

pdf(paste0("/Users/yeung/data/scchic/pdfs/Fig3_TF_analysis/tf_cluster_act_and_exprs.", Sys.Date(), ".ClstMeth.", jmeth, ".pdf"), useDingbats = FALSE)

myplclust(clusters, clusters$labels, cbPalette[clst])
m <- ggplot(left_join(cor.sum, clst.dat) %>% mutate(motif.lab = ifelse(zscore > zscore.cutoff, motif, NA)), 
            aes(x = cor.out, y = range.sqr, size = zscore, label = motif.lab, color = as.character(clstr))) + 
  geom_point(alpha = 1) + geom_text_repel(size = labsize) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Pearson Correlation") + ylab("Range of Activity * Range of Expression") + 
  scale_color_manual(values = cbPalette, na.value = "lightgrey", name = "Cluster")
print(m)

for (k in seq(K)){
  print(k)
  for (jmotif in names(clst[which(clst == k)])){
    print(jmotif)
    jjtitle <- paste(jmotif, "K:", k)
    m.act <- PlotXYWithColor(act.exprs.umap2 %>% filter(motif == jmotif), xvar = "X1", yvar = "X2", cname = "activity", jtitle = jjtitle, jsize = 2, jcol = "darkred")
    m.exprs <- PlotXYWithColor(act.exprs.umap2 %>% filter(motif == jmotif), xvar = "X1", yvar = "X2", cname = "exprs.log", jtitle = jjtitle, jsize = 2, jcol = "darkblue")
    m.cor <- ggplot(act.exprs.umap2 %>% filter(motif == jmotif & !is.na(traj)), aes(x = exprs.log, y = activity, color = traj)) + geom_point() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle(jjtitle) + geom_smooth(method = "lm", se = FALSE, inherit.aes = FALSE, mapping = aes(x = exprs.log, y = activity), color = "black")
    m.cor2 <- ggplot(act.exprs.umap2 %>% filter(motif == jmotif & !is.na(traj)), aes(x = exprs.log, y = activity, color = traj)) + geom_point() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle(jjtitle) + geom_smooth(method = "lm", se = FALSE, inherit.aes = TRUE)
    multiplot(m.act, m.exprs, cols = 2)
    print(m.cor)
    print(m.cor2)
  }
}
dev.off()
