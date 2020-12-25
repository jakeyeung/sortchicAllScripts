# Jake Yeung
# Date of Creation: 2020-11-20
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/3-bone_marrow_celltypes.from_glmpca.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)
library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

library(glmpca)


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

# Load LDA output ---------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmarks.k9 <- "H3K9me3"  # add this separately

hubprefix <- "/home/jyeung/hub_oudenaarden"

# niter <- "1000"
# binskeep <- 0
niter <- "500"
binskeep <- 1000
jsuffix <- paste0("bincutoff_0.binskeep_", binskeep, ".byplate.szname_none.niter_", niter, ".reorder_rownames.dupfilt")
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2"
# outname <- paste0("bonemarrow_celltypes.", Sys.Date(), ".niter_", niter, ".pdf")
outname <- paste0("BM_celltypes.", jsuffix, ".", Sys.Date(), ".pdf")
outf <- file.path(outdir, outname)

infs.lst <- lapply(jmarks, function(jmark){
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/glmpca.", jmark, ".", jsuffix, ".RData"))
  assertthat::assert_that(file.exists(inf))
  return(inf)
})


pdf(outf, useDingbats = FALSE)

dat.umap.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- infs.lst[[jmark]]
  # inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/glmpca.", jmark, ".", jsuffix, ".RData"))
  load(inf, v=T)
  dat.umap <- DoUmapAndLouvain(glm.out$factors, jsettings = jsettings)
  return(dat.umap)
})

# inf.k9 <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/glmpca.", jmarks.k9, ".bincutoff_0.binskeep_0.byplate.szname_none.niter_1000.reorder_rownames.dupfilt.RData"))
# assertthat::assert_that(file.exists(inf.k9))
# load(inf.k9, v=T)
# dat.umap.k9 <- DoUmapAndLouvain(glm.out$factors, jsettings = jsettings)


# Plot inits --------------------------------------------------------------

m.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  m <- ggplot(dat.umap.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    ggtitle(jmark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
})
print(m.lst)




# Add celltuypes ----------------------------------------------------------

dat.annots <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins/cell_cluster_table_with_spikeins.", jmark, ".2020-11-18.dupfilt.txt"))
  dat.annot <- fread(inf)
  return(dat.annot)
})



# match colors to k4me1 motif analysis ------------------------------------

inf.colors <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/motif_activity_H3K4me1.niter_1000.2020-11-20.uniquecolor.txt"
dat.colors.all <- fread(inf.colors)
dat.colors <- subset(dat.colors.all, select = c(cell, cluster, batch, col))

dat.colors.sum <- dat.colors %>%
  group_by(cluster) %>%
  summarise(col = unique(col))

dat.colors.sum.k9me3 <- dat.colors.sum %>%
  rowwise() %>%
  mutate(cluster = gsub(pattern = "Bcells", replacement = "Lymphoid", x = cluster),
         cluster = gsub(pattern = "HSPCs", replacement = "HSPCs", x = cluster),
         cluster = gsub(pattern = "Eryths", replacement = "Eryth", x = cluster)) %>%
  ungroup() %>%
  filter(cluster %in% c("Lymphoid", "HSPCs", "Eryth", "Granulocytes"))

dat.colors.sum.merge <- rbind(dat.colors.sum, dat.colors.sum.k9me3)

colhash <- hash::hash(dat.colors.sum.merge$cluster, dat.colors.sum.merge$col)


# K9me3 has specific clusternames -----------------------------------------

dat.umap.annot.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jtmp1 <- left_join(dat.umap.lst[[jmark]], subset(dat.annots[[jmark]], select = c(cell, cluster, batch, plate)), by = c("cell"))
  jtmp1$col <- sapply(jtmp1$cluster, function(x) AssignHash(x, colhash, null.fill = x))
  return(jtmp1)
  # if (jmark != "H3K9me3"){
  #   print("Not H3K9me3")
  #   jtmp1 <- left_join(jtmp1, dat.colors.sum, by = c("cluster"))
  # } else {
  #   print("Is H3K9me3")
  #   jtmp1 <- left_join(jtmp1, dat.colors.sum.k9me3, by = c("cluster"))
  # }
})

m.annot.lst2 <- lapply(jmarks, function(jmark){
  print(jmark)
  m <- ggplot(dat.umap.annot.lst[[jmark]], aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + 
    ggtitle(jmark) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_manual(values = cbPalette) 
  return(m)
})
print(m.annot.lst2)

m.annot.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  m <- ggplot(dat.umap.annot.lst[[jmark]], aes(x = umap1, y = umap2, color = col)) + 
    geom_point() + 
    ggtitle(jmark) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_identity(guide = "legend", labels = unique(dat.umap.annot.lst[[jmark]]$cluster), breaks = unique(dat.umap.annot.lst[[jmark]]$col), na.value = "grey85")
  return(m)
})
print(m.annot.lst)

dev.off()


# Write tables to output --------------------------------------------------

for (jmark in jmarks){
  print(jmark)
  outname.mark <- paste0("BM_celltypes.", jsuffix, ".", Sys.Date(), ".", jmark, ".txt")
  fwrite(dat.umap.annot.lst[[jmark]], file = file.path(outdir, outname.mark), sep = "\t", quote = FALSE)
}





# 
# # add celltypes -----------------------------------------------------------
# 
# jmarks
# 
# 
# inf2 <- file.path(hubprefix, "jyeung/data/scChiC/glmpca_outputs/glmpca.H3K4me1.bincutoff_0.binskeep_0.byplate.szname_none.niter_500.reorder_rownames.dupfilt.RData")
# 
# load(inf1, v=T)
# 
# glm.out1 <- glm.out
# 
# load(inf2, v=T)
# 
# glm.out2 <- glm.out
# 
# dat.umap1 <- DoUmapAndLouvain(glm.out1$factors, jsettings = jsettings)
# dat.umap2 <- DoUmapAndLouvain(glm.out2$factors, jsettings = jsettings)
# 
# ggplot(dat.umap1, aes(x = umap1, y = umap2, color = louvain)) + 
#   geom_point() + 
#   theme_bw() + ggtitle("1000iter") + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggplot(dat.umap2, aes(x = umap1, y = umap2, color = louvain)) + 
#   geom_point() + 
#   theme_bw() +  ggtitle("500iter") + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 



# 
# 
# 
# indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04"
# 
# 
# outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/BM_UMAP_celltypes.", Sys.Date(), ".glmpca.pdf")
# 
# pdf(file = outpdf, useDingbats = FALSE)
# 
# out.lst <- lapply(jmarks, function(jmark){
#   print(jmark)
#   fname <- paste0("lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj")
#   inf <- file.path(indir, fname)
#   assertthat::assert_that(file.exists(inf))
#   load(inf, v=T)
#   tm.result <- posterior(out.lda)
#   dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)
#   return(list(dat.umap = dat.umap, tm.result = tm.result, count.mat = count.mat))
# })
# 
# 
# dat.umap.lda <- lapply(jmarks, function(jmark){
#   jout <- out.lst[[jmark]]
#   jdat <- jout$dat.umap
#   jdat$mark <- jmark
#   return(jdat)
# }) %>%
#   bind_rows()
# 
# ggplot(dat.umap.lda %>% filter(mark == "H3K4me1"), aes(x = umap1, y = umap2)) + 
#   geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# # Load cell cluster -------------------------------------------------------
# 
# inf.rdata <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/LDA_downstream.BMround2all.merge_with_old.cleaner/LDA_downstream_objects.2020-11-05.again.RData"
# load(inf.rdata, v=T)
# 
# dat.merge <- dat.merge %>%
#   rowwise() %>%
#   mutate(cluster = ifelse(cluster == "Eryth", "Eryths", cluster))
# # merge
# dat.umap.lda.annot <- left_join(dat.umap.lda, subset(dat.merge, select = c(cell, cluster), by = c("cell")))
# 
# ctypes <- sort(unique(subset(dat.umap.lda.annot, mark == "H3K4me1")$cluster))
# ctypes2 <- gsub(pattern = "Bcells", replacement = "Lymphoid", ctypes)
# # ctypes2 <- sort(unique(subset(dat.umap.lda.annot, mark == "H3K9me3")$cluster))
# 
# cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
# m.lst <- lapply(jmarks, function(jmark){
#   print(jmark)
#   jsub <- dat.umap.lda.annot %>% filter(mark == jmark) 
#   if (jmark != "H3K9me3"){
#     jsub$cluster <- factor(jsub$cluster, levels = ctypes)
#     jctypes <- ctypes
#   } else {
#     jsub$cluster <- factor(jsub$cluster, levels = ctypes2)
#     jctypes <- ctypes2
#   }
#   jpal <- cbPalette[sort(unique(as.numeric(jsub$cluster)))]
#   jsub$jcol <- cbPalette[as.numeric(jsub$cluster)]
#   m <- ggplot(jsub, aes(x = umap1, y = umap2, color = cluster)) + 
#     geom_point() + 
#     ggtitle(jmark) + 
#     theme_minimal() + 
#     scale_color_manual(values = jpal) + 
#     # scale_color_identity(labels = unique(jsub$cluster), guide = "legend") + 
#     theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#   return(m)
# })
# 
# print(m.lst)
# 
# dev.off()

