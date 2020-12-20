# Jake Yeung
# Date of Creation: 2020-12-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/2-find_K4me1_and_K9me3_bins.more_bins.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


jmarks <- c("H3K4me1", "H3K9me3"); names(jmarks) <- jmarks

# Load K4me1 and K9me3 cuts ------------------------------------------------

inf.out <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.H3K4me1_H3K9me3.reseq/K4me1_K9me3_objs.more_bins.twofilt.RData")
inf.out.count.mat1 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.H3K4me1_H3K9me3.reseq/H3K4me1_count_mat_all_bins.rds"
inf.out.count.mat2 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.H3K4me1_H3K9me3.reseq/H3K9me3_count_mat_all_bins.rds"

pdf.out <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.H3K4me1_H3K9me3.reseq/K4me1_K9me3_objs.more_bins.", Sys.Date(), ".pdf")

pdf(pdf.out, useDingbats = FALSE)



# 50 kb
hubprefix <- "/home/jyeung/hub_oudenaarden"
inf.lda.k4me1 <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmarks[[1]], ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmarks[[1]], ".remove_bad_clusters.2020-11-04.K-30.Robj"))
load(inf.lda.k4me1, v=T)

out.lda.k4me1 <- out.lda
count.mat.k4me1 <- count.mat

inf.lda.k9me3 <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmarks[[2]], ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmarks[[2]], ".remove_bad_clusters.2020-11-04.K-30.Robj"))
load(inf.lda.k9me3, v=T)
out.lda.k9me3 <- out.lda
count.mat.k9me3 <- count.mat


# Take common rows --------------------------------------------------------

rnames.all <- unique(unlist(lapply(list(count.mat.k4me1, count.mat.k9me3), FUN = function(x) rownames(x))))
# rnames.all <- unlist(unique(c(rownames(count.mat.k4me1))))
count.mat.merge <- cbind.fill.lst(mats.lst = list(count.mat.k4me1, count.mat.k9me3), all.rnames = rnames.all)

count.mat.k4me1.split <- count.mat.merge[, colnames(count.mat.k4me1)]
count.mat.k9me3.split <- count.mat.merge[, colnames(count.mat.k9me3)]


# Load meta dat -----------------------------------------------------------

inf.k4me1 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.H3K4me1.txt")
dat.k4me1 <- fread(inf.k4me1)

inf.k9me3 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/BM_celltypes.bincutoff_0.binskeep_1000.byplate.szname_none.niter_500.reorder_rownames.dupfilt.2020-11-23.H3K9me3.txt")
dat.k9me3 <- fread(inf.k9me3)

# cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.k9me3, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Get pseudobulk ----------------------------------------------------------

cell.clstr.k4me1 <- split(dat.k4me1$cell, dat.k4me1$cluster)
cell.clstr.k9me3 <- split(dat.k9me3$cell, dat.k9me3$cluster)

pseudobulk.k4me1 <- SumAcrossClusters(count.mat.k4me1.split, cnames.keep.lst = cell.clstr.k4me1)
pseudobulk.k4me1 <- as.data.frame(do.call(cbind, pseudobulk.k4me1))


pseudobulk.k9me3 <- SumAcrossClusters(count.mat.k9me3.split, cnames.keep.lst = cell.clstr.k9me3)
pseudobulk.k9me3 <- as.data.frame(do.call(cbind, pseudobulk.k9me3))

plot(pseudobulk.k4me1$Eryths, pseudobulk.k9me3$Eryth, pch = 20)

plot(pseudobulk.k4me1$Bcells, pseudobulk.k9me3$Lymphoid, pch = 20)

plot(pseudobulk.k4me1$Granulocytes, pseudobulk.k9me3$Granulocytes, pch = 20)

# plot(pseudobulk.k4me1$Bcells, pseudobulk.k9me3$Eryth, pch = 20)

plot(density(unlist(log2(pseudobulk.k9me3$Eryth))))
plot(density(unlist(log2(pseudobulk.k9me3$Granulocytes))))
plot(density(unlist(log2(pseudobulk.k9me3$Lymphoid))))
plot(density(unlist(log2(pseudobulk.k9me3$HSPCs))))

plot(density(unlist(log2(pseudobulk.k4me1$Eryths))))
plot(density(unlist(log2(pseudobulk.k4me1$Granulocytes))))
plot(density(unlist(log2(pseudobulk.k4me1$Bcells))))
plot(density(unlist(log2(pseudobulk.k4me1$HSPCs))))

colnames(pseudobulk.k9me3) <- paste("H3K9me3", colnames(pseudobulk.k9me3), sep = "_")
colnames(pseudobulk.k4me1) <- paste("H3K4me1", colnames(pseudobulk.k4me1), sep = "_")

pseudobulk <- cbind(pseudobulk.k4me1, pseudobulk.k9me3)
pseudobulk2 <- as.data.frame(cbind(k4me1 = rowSums(pseudobulk.k4me1), k9me3 = rowSums(pseudobulk.k9me3)))
pseudobulk2$rnames <- rownames(pseudobulk2)



plot(log2(pseudobulk2$k4me1), log2(pseudobulk2$k9me3), pch = 20)
abline(v = 8, col = 'blue')
abline(h = 8, col = 'red')

plot(density(log2(pseudobulk2$k4me1)))
abline(v = 8, col = 'red')
plot(density(log2(pseudobulk2$k9me3)))
abline(v = 8, col = 'blue')

pseudobulk.norm <- sweep(pseudobulk, MARGIN = 2, STATS = colSums(pseudobulk), FUN = "/")


# Define heterochromatin and active regions -------------------------------

# in log2
min.k9me3 <- 8
# max.k9me3 <- 8

min.k4me1 <- 8
# max.k4me1 <- 11

bins.k4 <- subset(pseudobulk2, log2(k9me3) < min.k9me3 & log2(k4me1) > min.k4me1)
bins.k9 <- subset(pseudobulk2, log2(k4me1) < min.k4me1 & log2(k9me3) > min.k9me3)


# Define bins  ------------------------------------------------------------

# saveRDS(object = count.mat.k4me1.split, file = inf.out.count.mat1)
# saveRDS(object = count.mat.k9me3.split, file = inf.out.count.mat2)
# save(pseudobulk, pseudobulk2, min.k9me3, min.k4me1, bins.k4, bins.k9, dat.k4me1, dat.k9me3, file = inf.out)



dev.off()


# Make heatmap  -----------------------------------------------------------


