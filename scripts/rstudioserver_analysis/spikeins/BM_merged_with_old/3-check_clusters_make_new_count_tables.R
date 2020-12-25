# Jake Yeung
# Date of Creation: 2020-11-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/3-filter_H3K27me3.R
# rerun LDA

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

hubprefix <- "/home/jyeung/hub_oudenaarden"
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

indir.mat <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs.merge_with_new_BM"
outdir.mat <- file.path(indir.mat, paste0("clusterfilt.", Sys.Date()))
dir.create(outdir.mat)
outpdf <- file.path(outdir.mat, paste0("cluster_filt_UMAPs.", Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)

# inf <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs.merge_with_new_BM/count_mat_old_merged_with_new.H3K27me3.rds")
# outf <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs.merge_with_new_BM/clusterfilt/count_mat_old_merged_with_new.H3K27me3.rds")
# assertthat::assert_that(file.exists(inf))
# count.mat <- readRDS(inf)


# Load annots -------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

inf.annot1 <- lapply(jmarks, function(jmark){
  inf.tmp <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios/spikeins_mouse.BMround2_umaps_and_ratios.colfix.celltyping.2020-10-31.WithRelLevels.mark_", jmark, ".cell_cluster_tables.txt")
})
dat.annot1 <- lapply(inf.annot1, fread) %>%
  bind_rows()


inf.annot2 <- lapply(jmarks, function(jmark){
  inf.tmp <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2/spikeins_mouse.BMround1and2_umaps_and_ratios.colfix.celltyping.2020-11-01.WithRelLevels.mark_", jmark, ".cell_cluster_tables.txt")
})
dat.annot2 <- lapply(inf.annot2, fread) %>%
  bind_rows()


# 
# jmark <- "H3K27me3"
# jmark <- "H3K4me3"
# jmark <- "H3K4me1"
# 
# ggplot(dat.annot1 %>% filter(mark == jmark), aes(x = umap1, y = umap2, color = cluster)) + 
#   geom_point() + 
#   theme_bw() + 
#   scale_color_manual(values = cbPalette) + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# ggplot(dat.annot2 %>% filter(mark == jmark), aes(x = umap1, y = umap2, color = cluster)) + 
#   geom_point() + 
#   theme_bw() + 
#   scale_color_manual(values = cbPalette) + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# show louvain clusters onto first umap
dat.annot.merge <- left_join(dat.annot1, dat.annot2, by = c("cell"))

m.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  m1 <- ggplot(dat.annot.merge %>% filter(mark.x == jmark), aes(x = umap1.x, y = umap2.x, color = cluster.y)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_manual(values = cbPalette, na.value = "grey95") +  
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
  m2 <- ggplot(dat.annot.merge %>% filter(mark.x == jmark), aes(x = umap1.y, y = umap2.y, color = cluster.y)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_manual(values = cbPalette, na.value = "grey95") + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  JFuncs::multiplot(m1, m2, cols = 2)
  return(list(m1, m2))
})





# Clean up k4me1 ----------------------------------------------------------



bad.clsts.k4me1 <- c("Macrophages", NA)
merge.clsts.k4me1 <- list("BcellsNaive" = "Bcells",
                    "BcellsPlasma" = "Bcells")
merge.clsts.hash.k4me1 <- hash(merge.clsts.k4me1)

dat.k4me1.new <- subset(dat.annot2, (mark == "H3K4me1") & (!cluster %in% bad.clsts.k4me1)) %>%
  rowwise() %>%
  mutate(cluster = AssignHash(x = cluster, jhash = merge.clsts.hash.k4me1, null.fill = cluster))

ggplot(subset(dat.annot2, mark == "H3K4me1"), aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("H3K4me1", "before removing clusters") + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.k4me1.new, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("H3K4me1", paste(bad.clsts.k4me1, collapse = ",")) + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

JFuncs::multiplot(m.lst$H3K4me1[[1]], m.lst$H3K4me1[[2]], cols = 2)

# Clean up k4me3 ----------------------------------------------------------


bad.clsts.k4me3 <- c(NA)
merge.clsts.k4me3 <- list("BcellsNaive" = "Bcells",
                          "BcellsPlasma" = "Bcells")
merge.clsts.hash.k4me3 <- hash(merge.clsts.k4me3)

dat.k4me3.new <- subset(dat.annot2, (mark == "H3K4me3") & (!cluster %in% bad.clsts.k4me3)) %>%
  rowwise() %>%
  mutate(cluster = AssignHash(x = cluster, jhash = merge.clsts.hash.k4me3, null.fill = cluster))

ggplot(subset(dat.annot2, mark == "H3K4me3"), aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("H3K4me3", "before removing clusters") + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.k4me3.new, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("H3K4me3", paste(bad.clsts.k4me3, collapse = ",")) + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

JFuncs::multiplot(m.lst$H3K4me3[[1]], m.lst$H3K4me3[[2]], cols = 2)

# Clean up k27me3 ---------------------------------------------------------

bad.clsts.k27me3 <- c("BasophilsCheck", "louvain1", "louvain3", "louvain5", "louvain18", NA)
merge.clsts.k27me3 <- list("pDCsCheck" = "pDCs")
merge.clsts.hash.k27me3 <- hash("pDCsCheck" = "pDCs")

# dat.k27me3.new.check <- subset(dat.annot2, (mark == "H3K27me3") & (!cluster %in% bad.clsts.k27me3)) %>%
#   rowwise() %>%
#   mutate(cluster = AssignHash(x = cluster, jhash = merge.clsts.hash.k27me3, null.fill = cluster))

dat.k27me3.new <- subset(dat.annot2, (mark == "H3K27me3") & (!cluster %in% bad.clsts.k27me3)) %>%
  rowwise() %>%
  mutate(cluster = ifelse(cluster == "pDCsCheck", "pDCs", cluster))

# xtest <- sapply(dat.k27me3.new$cluster, function(x) AssignHash(x = x, jhash = merge.clsts.hash.k27me3, null.fill = x))
# xtest <- sapply(dat.k27me3.new$cluster, function(x) x)

ggplot(subset(dat.annot2, mark == "H3K27me3"), aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("H3K27me3", "before removing clusters") + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.k27me3.new, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("H3K27me3", paste(bad.clsts.k27me3, collapse = ",")) + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


JFuncs::multiplot(m.lst$H3K27me3[[1]], m.lst$H3K27me3[[2]], cols = 2)


# K9me3 no need to clean up  ----------------------------------------------

dat.k9me3.new <- subset(dat.annot2, mark == "H3K9me3")

ggplot(subset(dat.annot2, mark == "H3K9me3"), aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("H3K9me3", "before removing clusters") + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.k9me3.new, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("H3K9me3", "after:") + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

JFuncs::multiplot(m.lst$H3K9me3[[1]], m.lst$H3K9me3[[2]], cols = 2)

# Combine new dat annots --------------------------------------------------

dat.annots.new <- list("H3K4me1" = dat.k4me1.new, "H3K4me3" = dat.k4me3.new, "H3K27me3" = dat.k27me3.new, "H3K9me3" = dat.k9me3.new)

# Clean up count matrix ---------------------------------------------------


mats.unfilt <- lapply(jmarks, function(jmark){
  inf <- file.path(indir.mat, paste0("count_mat_old_merged_with_new.", jmark, ".rds"))
  readRDS(inf)
})

mats.filt <- lapply(jmarks, function(jmark){
  print(jmark)
  jmat.tmp <- mats.unfilt[[jmark]]
  cells.keep <- dat.annots.new[[jmark]]$cell
  jmat.filt.tmp <- jmat.tmp[, cells.keep]
  print(dim(jmat.tmp))
  print(dim(jmat.filt.tmp))
  return(jmat.filt.tmp)
})

# write output
outfs <- lapply(jmarks, function(jmark){
  print(jmark)
  outf <- file.path(outdir.mat, paste0("count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.", Sys.Date(), ".rds"))
  if (!file.exists(outf.table.tmp)){
    saveRDS(mats.filt[[jmark]], file = outf)
  } else {
    print(paste("Exists. Skipping", outf))
  }
  return(outf)
})

# write new cell table to output

lapply(jmarks, function(jmark){
  print(jmark)
  outf.table.tmp <- file.path(outdir.mat, paste0("cell_cluster_table.old_merged_with_new.", jmark, ".remove_bad_clusters.", Sys.Date(), ".txt"))
  print(dim(dat.annots.new[[jmark]]))
  if (!file.exists(outf.table.tmp)){
    fwrite(dat.annots.new[[jmark]], file = outf.table.tmp, quote = FALSE, sep = "\t", na = "NA")
  } else {
    print(paste("Exists. Skipping", outf.table.tmp))
  }
})

# split by batch
batch2cells <- unique(dat.annot1$cell)
batch1cells <- subset(dat.annot2, !cell %in% batch2cells)$cell

# write output
outfs <- lapply(jmarks, function(jmark){
  print(jmark)
  outf <- file.path(outdir.mat, "batch1", paste0("count_mat_BM_batch1.", jmark, ".remove_bad_clusters.", Sys.Date(), ".rds"))
  if (!file.exists(outf)){
    jmat <- mats.filt[[jmark]]
    cells.keep <- colnames(jmat) %in% batch1cells
    jmat.filt <- jmat[, cells.keep]
    print(dim(jmat))
    print(dim(jmat.filt))
    saveRDS(jmat.filt, file = outf)
  } else {
    print(paste("Exists. Skipping", outf))
  }
  return(outf)
})

# write output
outfs <- lapply(jmarks, function(jmark){
  print(jmark)
  outf <- file.path(outdir.mat, "batch2", paste0("count_mat_BM_batch2.", jmark, ".remove_bad_clusters.", Sys.Date(), ".rds"))
  if (!file.exists(outf)){
    jmat <- mats.filt[[jmark]]
    cells.keep <- colnames(jmat) %in% batch2cells
    jmat.filt <- jmat[, cells.keep]
    print(dim(jmat))
    print(dim(jmat.filt))
    saveRDS(jmat.filt, file = outf)
  } else {
    print(paste("Exists. Skipping", outf))
  }
  return(outf)
})


dev.off()