# Jake Yeung
# Date of Creation: 2022-01-10
# File: ~/projects/scchic/scripts/macbook_analysis_2021/mara_downstream/check_sitecounts_H3K27me3_H3K9me3.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)


# Functions ---------------------------------------------------------------

MeanAcrossClusters <- function(count.mat, cnames.keep.lst){
  count.mat <- as.matrix(count.mat)
  count.vecs <- lapply(cnames.keep.lst, function(cnames.keep){
    cnames.keep.i <- which(colnames(count.mat) %in% cnames.keep)
    assertthat::assert_that(length(cnames.keep.i) > 0)
    rowMeans(count.mat[, cnames.keep.i])
  })
  return(count.vecs)
}




# Load K27me3 -------------------------------------------------------------

infN.k27 <- "/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/mara_objs/sitecount_mats/dynamic_bins_50kb__motevo_merged.closest.long.scale_0.center_0.byrow_0.H3K27me3.YyMerged.tab.txt.gz"
N.k27 <- fread(infN.k27)

Nsub.k27 <- subset(N.k27, select = c(Gene.ID, YyMerged)) %>%
  arrange(desc(YyMerged))

plot(density(Nsub.k27$YyMerged))

print(head(Nsub.k27))

# good regions
# chr2:172200000-172250000
# chr13:60050000-60100000




# Plot K9me3 --------------------------------------------------------------

infN.k9 <- "/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/mara_objs/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.H3K9me3.txt.gz"
N.k9 <- fread(infN.k9)

# Plzf == Zbtb16
Nsub.k9 <- subset(N.k9, select = c(Gene.ID, Zbtb16)) %>%
  arrange(desc(Zbtb16))
plot(density(Nsub.k9$Zbtb16))

print(head(Nsub.k9))

infE.k9 <- "/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/mara_objs/exprs_mats/ldaOut.count_mat_k9_dynamic_bins_50kb.H3K9me3.2021-01-28.K-30.keepNbins_0.txt.gz"
E.k9 <- as.data.frame(fread(infE.k9))
rownames(E.k9) <- E.k9$Gene.ID
E.k9$Gene.ID <- NULL

infE.k27 <- "/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/mara_objs/exprs_mats/ldaOut.count_tables_merged.H3K27me3.DE_bins_all_marks_top_6085_dists_to_TSS.annot_table.H3K27me3.2021-02-15.txt.K-30.keepNbins_0.withchr2.txt.gz"
E.k27 <- as.data.frame(fread(infE.k27))
rownames(E.k27) <- E.k27$Gene.ID
E.k27$Gene.ID <- NULL


# Get metadata  -----------------------------------------------------------

inf.meta.k9 <- "/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/mara_objs/metadata/metadata_batch_corrected.arranged_by_lineage.shuffled.H3K9me3.2021-02-19.txt"
dat.meta.k9 <- fread(inf.meta.k9)

inf.meta.k27 <- "/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/mara_objs/metadata/metadata_batch_corrected.arranged_by_lineage.shuffled.H3K27me3.2021-02-19.txt"
dat.meta.k27 <- fread(inf.meta.k27)

cnames.keep.k9.lst <- lapply(split(dat.meta.k9, dat.meta.k9$cluster), function(x) make.names(x$cell))
cnames.keep.k27.lst <- lapply(split(dat.meta.k27, dat.meta.k27$cluster), function(x) make.names(x$cell))

E.k9.mean <- MeanAcrossClusters(E.k9, cnames.keep.k9.lst) %>%
  do.call(cbind, .) 
E.k9.mean.long <- E.k9.mean %>%
  melt()
colnames(E.k9.mean.long) <- c("Gene.ID", "cluster", "chicsignal")
# E.k9.mean$Gene.ID <_ rownames(E.k9.mean)
E.k9.mean.long <- E.k9.mean.long %>%
  group_by(cluster) %>%
  mutate(zscore = scale(chicsignal, center = TRUE, scale = TRUE))

E.k27.mean <- MeanAcrossClusters(E.k27, cnames.keep.k27.lst) %>%
  do.call(cbind, .) 
E.k27.mean.long <- E.k27.mean %>%
  melt()
colnames(E.k27.mean.long) <- c("Gene.ID", "cluster", "chicsignal")

E.k27.mean.long <- E.k27.mean.long %>%
  group_by(cluster) %>%
  mutate(zscore = scale(chicsignal, center = TRUE, scale = TRUE))

# plot differences for differenet regions
jregions.k9 <- Nsub.k9$Gene.ID[1:100]
jregions.k27 <- Nsub.k27$Gene.ID[1:100]

ggplot(E.k9.mean.long %>% filter(Gene.ID %in% jregions.k9),  
       aes(x = cluster, y = zscore)) + 
  geom_boxplot () +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(E.k27.mean.long %>% filter(Gene.ID %in% jregions.k27),  
       aes(x = cluster, y = zscore)) + 
  geom_boxplot () +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())





