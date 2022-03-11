# Jake Yeung
# Date of Creation: 2022-01-04
# File: ~/projects/scchic/scripts/macbook_analysis_2021/double_incubated_analysis/4-check_scchix_from_celltype_annots.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)
library(ggforce)
library(JFuncs)
library(scchicFuncs)


outpdf <- paste0("/Users/yeung/data/scchic/from_cluster_2021/primetime_plots/coords_dbl_assignment_knn_impute_with_UMAP.", Sys.Date(), ".pdf")
outrds <- paste0("/Users/yeung/data/scchic/from_cluster_2021/primetime_plots/coords_dbl_assignment_knn_impute.", Sys.Date(), ".rds")
pdf(outpdf, useDingbats = FALSE)

# load metas --------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K9me3"); names(jmarks) <- jmarks
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

inf.meta <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/from_cluster_2021/double_stain_outputs/scchix_outputs.H3K4me1-H3K9me3.2021-02-11.setseed.RData"
load(inf.meta, v=T)

cells.keep <- dat.umap.long.annot$cell


m.grid <- ggplot(dat.umap.long.annot, aes(x = louv.act.impute, y = louv.repress.impute, color = louv.act.impute)) +
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  theme(aspect.ratio=0.6) +
  scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) +
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid)

# Load outputs ------------------------------------------------------------


wmax <- c("0.5", "0.9", "0.99"); names(wmax) <- wmax

fits.out.lst <- lapply(wmax, function(w){
  inf <- paste0("/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/snakemake_outputs/snakemake_wlimits.0.01_", w, "/snakemake_outputs/scchix_outputs_objs/scchix_inputs_clstr_by_celltype_H3K4me1-H3K9me3.RData")
  if (w == "0.99"){
    inf <- "/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/count_tables.BM.k4_k9_dynamic_regions/scchix_outputs.50kb_genomewide/unmix_scchix_inputs_clstr_by_celltype_H3K4me1xH3K9me3.removeNA_FALSE.RData"
  }
  load(inf, v=T)
  return(act.repress.coord.lst)
})


coords.dbl.lst <- lapply(fits.out.lst, function(fits.out){
  # if louvains are now from clusters need eto rethink jcoord
  cell.vec <- names(fits.out)
  names(cell.vec) <- cell.vec
  coords.dbl <- lapply(cell.vec, function(jcell){
    jfit <- fits.out[[jcell]]
    jweight <- fits.out[[jcell]]$w
    p.mat <- SoftMax(jfit$ll.mat)
    jcoord <- which(jfit$ll.mat == max(jfit$ll.mat), arr.ind = TRUE)
    jmax <- max(p.mat)
    
    # rows are active, columns are repress I THINK?
    # TODO: assumes underscores be careful!
    jlouv.act <- rownames(p.mat)[[jcoord[[1]]]]
    jlouv.repress <- colnames(p.mat)[[jcoord[[2]]]]
    
    if (grepl("_", jlouv.act)){
      jlouv.act <- strsplit(jlouv.act, split = "_")[[1]][[2]]
    }
    if (grepl("_", jlouv.repress)){
      jlouv.repress <- strsplit(jlouv.repress, split = "_")[[1]][[2]]
    }
    
    out.dat <- data.frame(cell = jcell, louv.act = jlouv.act, louv.repress = jlouv.repress, lnprob = jmax, w = jweight, stringsAsFactors = FALSE)
    return(out.dat)
  }) %>%
    bind_rows()
  return(coords.dbl)
})

# coords.dbl <- rbind(subset(coords.dbl.lst$`0.9` %>% filter(louv.act == "Eryths")), subset(coords.dbl.lst$`0.9` %>% filter(louv.act != "Eryths")))

coords.dbl <- coords.dbl.lst$`0.99`

# plot distribution of probabilities
plot(density(exp(coords.dbl$lnprob)))

ggplot(coords.dbl, aes(x = exp(lnprob))) + 
  geom_histogram(bins = 100) + 
  xlab("Probability of Assigning to Best Model") + 
  theme_bw(24) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m.grid <- ggplot(coords.dbl, aes(x = louv.act, y = louv.repress, color = louv.act)) +
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  theme(aspect.ratio=0.6) +
  scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) +
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid)


# Clean up with imputes ---------------------------------------------------

inf.lda <- "/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/count_tables.BM.k4_k9_dynamic_regions/lda_and_clusters.50kb_genomewide/ClusterAnnot.lda_and_datmerged.k4_k9_50kb_genomewide.H3K4me1xH3K9me3.2021-01-31.RData"
load(inf.lda, v=T)


# Get nearest neighbors ---------------------------------------------------

library(hash)
library(igraph)
library(umap)
library(topicmodels)

tm.result <- posterior(out.lda)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
jsettings$spread <- 8

topics.mat <- tm.result$topics
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)




cell2repress <- hash::hash(dat.umap.long.annot$cell, dat.umap.long.annot$louv.repress)
cell2act <- hash::hash(dat.umap.long.annot$cell, dat.umap.long.annot$louv.act)

knn.out.repress <- umap.out$knn$indexes
rnames.orig.repress <- rownames(knn.out.repress); names(rnames.orig.repress) <- rnames.orig.repress
rownames(knn.out.repress) <- sapply(rownames(knn.out.repress), function(x) AssignHash(x, cell2repress, null.fill = x))

knn.out.act <- umap.out$knn$indexes
rnames.orig.act <- rownames(knn.out.act); names(rnames.orig.act) <- rnames.orig.act
rownames(knn.out.act) <- sapply(rownames(knn.out.act), function(x) AssignHash(x, cell2act, null.fill = x))


ctype.neighbors.counts.repress <- sapply(rnames.orig.repress, function(jcell){
  indx <- which(rnames.orig.repress == jcell)
  xsummary <- sort(table(rownames(knn.out.repress)[knn.out.repress[indx, ]]), decreasing = TRUE)
  return(names(xsummary)[[1]])
})

ctype.neighbors.counts.act <- sapply(rnames.orig.act, function(jcell){
  indx <- which(rnames.orig.act == jcell)
  xsummary <- sort(table(rownames(knn.out.act)[knn.out.act[indx, ]]), decreasing = TRUE)
  return(names(xsummary)[[1]])
})

coords.dbl.impute <- coords.dbl
coords.dbl.impute$louv.repress.impute <- sapply(coords.dbl.impute$cell, function(x) ctype.neighbors.counts.repress[[x]])
coords.dbl.impute$louv.act.impute <- sapply(coords.dbl.impute$cell, function(x) ctype.neighbors.counts.act[[x]])

m.grid <- ggplot(coords.dbl.impute, aes(x = louv.act.impute, y = louv.repress.impute, color = louv.act.impute)) +
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  theme(aspect.ratio=0.6) +
  scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Labels")) + ylab(paste0(jmarks[[2]], " Labels")) +
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid)

# rearrange clusters and write

# act.lvls <- levels(dat.umap.long.annot$louv.act.impute)
# repress.lvls <- levels(dat.umap.long.annot$louv.repress.impute)

act.lvls <- c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "DCs", "pDCs", "HSPCs"); names(act.lvls) <- act.lvls
repress.lvls <- c("Eryths", "Bcells", "Granulocytes", "HSPCs"); names(repress.lvls) <- repress.lvls

cluster2col.act <- hash::hash(dat.umap.long.annot$louv.act.impute, dat.umap.long.annot$clustercol.act)
cluster2col.repress <- hash::hash(dat.umap.long.annot$louv.repress.impute, dat.umap.long.annot$clustercol.repress)

coords.dbl.impute$louv.act.impute <- factor(coords.dbl.impute$louv.act.impute, levels = act.lvls)
coords.dbl.impute$louv.repress.impute <- factor(coords.dbl.impute$louv.repress.impute, levels = repress.lvls)

coords.dbl.impute$clustercol.act <- sapply(coords.dbl.impute$louv.act.impute, AssignHash, cluster2col.act)
coords.dbl.impute$clustercol.repress <- sapply(coords.dbl.impute$louv.repress.impute, AssignHash, cluster2col.repress)


cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# colorcode by H3K9me3 clusters
m.grid <- ggplot(coords.dbl.impute, aes(x = louv.act.impute, y = louv.repress.impute, color = clustercol.repress)) +
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  theme(aspect.ratio=0.6) +
  scale_color_identity() + 
  xlab(paste0(jmarks[[1]], " Clusters")) + ylab(paste0(jmarks[[2]], " Clusters")) +
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid)

# show UMAP 
dat.umap.long.annot.genomewide <- left_join(coords.dbl.impute, dat.umap.long, by = "cell")

m.grid <- ggplot(dat.umap.long.annot.genomewide, aes(x = louv.act.impute, y = louv.repress.impute, color = clustercol.repress)) +
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  theme(aspect.ratio=0.6) +
  scale_color_identity() + 
  xlab(paste0(jmarks[[1]], " Clusters")) + ylab(paste0(jmarks[[2]], " Clusters")) +
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid)


m.umap.act <- ggplot(dat.umap.long.annot.genomewide, aes(x = umap1, y = umap2, color = clustercol.act)) + 
  geom_point() + 
  theme_minimal() + 
  scale_color_identity() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m.umap.repress <- ggplot(dat.umap.long.annot.genomewide, aes(x = umap1, y = umap2, color = clustercol.repress)) + 
  geom_point() + 
  theme_minimal() + 
  scale_color_identity() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(m.umap.act)
print(m.umap.repress)

dev.off()



saveRDS(dat.umap.long.annot.genomewide, file = outrds)




