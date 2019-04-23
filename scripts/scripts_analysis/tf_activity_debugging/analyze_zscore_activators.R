# Jake Yeung
# Date of Creation: 2019-04-22
# File: ~/projects/scchic/scripts/scripts_analysis/tf_activity_debugging/analyze_zscore_activators.R
# look at Zscores of activators

rm(list=ls())

library(dplyr)
library(ggplot2)
library(ggrepel)

source("scripts/Rfunctions/GetMetaData.R")
source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/MatchCellNameToSample.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/TrajFunctions.R")
source("scripts/Rfunctions/AnnotationFunctions.R")


# Constants ---------------------------------------------------------------

jscale.fac <- 10^6
jpseudo <- 0

inf.trajs <- "/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-11.RData"
inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"

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
mreps <- dat.trajs.long %>% filter(mark == "H3K4me3")
mreps.uniq <- unique(paste(sapply(mreps$cell, function(x) strsplit(x, "_")[[1]][[3]]), sapply(mreps$cell, function(x) strsplit(x, "_")[[1]][[4]]), sep = "_"))

act.uniq <- unique(paste(sapply(mara.out2$act.long$cell, function(x) strsplit(x, "_")[[1]][[3]]), sapply(mara.out2$act.long$cell, function(x) strsplit(x, "_")[[1]][[4]]), sep = "_"))

jmarks <- c("H3K4me1", "H3K4me3")
dat.sub <- dat.trajs.long %>% filter(mark %in% jmarks)

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
(motifs.keep <- subset(zscores.merged, zscore.x > 2)$motif)
jmat <- subset(mara.out1$act.mat, motif %in% motifs.keep)
jmotifs <- jmat$motif
jmat$motif <- NULL
jmat <- as.matrix(scale(jmat, center=FALSE, scale=TRUE))
rownames(jmat) <- jmotifs

clusters <- hclust(dist(scale(jmat)), method = "ward.D")
plot(clusters)

# heatmap
# plot heatmap
library(heatmap3)
library(gplots)
library(made4)

# reorder columns by lambda 
cell.ordering.dat <- traj.annot %>%
  filter(grepl(jmark, cell)) %>%
  group_by(traj) %>%
  arrange(traj, lambda)
cells <- cell.ordering.dat$cell

pdf(paste0("~/data/scchic/pdfs/heatmap_global_analysis_motifs_H3K4me1.", Sys.Date(), ".pdf"), useDingbats = FALSE)
  heatmap3(t(jmat[, cells]), margins = c(5, 8), cexRow = 0.25, Colv = TRUE, Rowv = NA, labRow = FALSE, scale = "column", revC = TRUE)
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

genes <- sapply(rownames(mat.impute.sub), function(x){
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

jmotifs <- c("Cebpb", "Tal1", "Ebf1")
jpaths <- c("granu", "eryth", "lymphoid")



pdf(paste0("~/data/scchic/pdfs/TF_analysis_activity_umap_Figure3.", Sys.Date(), ".pdf"), useDingbats = FALSE)
for (i in seq(length(jmotifs))){
  jmotif <- jmotifs[[i]]
  jpath <- jpaths[[i]]
  print(jmotif)
  m1 <- PlotXYWithColor(act.exprs.umap2 %>% filter(motif == jmotif), xvar = "X1", yvar = "X2", cname = "activity", jcol = scales::muted("darkred")) + ggtitle(jmotif) + 
    geom_path(data = trajs.spring[[jmark]][[jpath]], inherit.aes = FALSE, mapping = aes(x = umap1, y = umap2), size = jsize, color = "black")
  print(m1)
}
dev.off()

# Show correlation with pseudotime  ---------------------------------------

act.exprs.umap2 <- left_join(act.exprs.umap2, traj.annot)
act.exprs.umap2$exprs.log <- log2(act.exprs.umap2$exprs * jscale.fac + jpseudo)

# show positive correlations and negative ones 

jmotif <- "Foxc1"
jmotif <- "Tal1"
jmotif <- "Cebpb"
jmotif <- "Erbf1"
jmotif <- "Spib"

# plot globa
m.global <- ggplot(mara.out1$zscores %>% mutate(motif.lab = ifelse(zscore > 2, motif, NA), motif = factor(motif, levels = motif)), 
       aes(x = motif, y = zscore, label = motif.lab)) + geom_point(alpha = 0.2, size = 0.2)  +  geom_text_repel() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  xlab("")

jcols <- c("#A29F9D", "#6AB7E6", "#E7A100")
jmotifs <- motifs.keep
pdf(paste0("~/data/scchic/pdfs/TF_analysis_gene_correlations_pseudotime_Figure3.", Sys.Date(), ".pdf"), useDingbats = FALSE)
for (jmotif in jmotifs){
  print(m.global)
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

