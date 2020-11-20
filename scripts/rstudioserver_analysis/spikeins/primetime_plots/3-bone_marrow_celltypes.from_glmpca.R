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

# Load LDA output ---------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmarks.k9 <- "H3K9me3"  # add this separately

hubprefix <- "/home/jyeung/hub_oudenaarden"

dat.umap.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/glmpca.", jmark, ".bincutoff_0.binskeep_0.byplate.szname_none.niter_1000.reorder_rownames.dupfilt.RData"))
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
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
})



# Add celltuypes ----------------------------------------------------------

dat.annots <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins/cell_cluster_table_with_spikeins.", jmark, ".2020-11-18.dupfilt.txt"))
  dat.annot <- fread(inf)
  return(dat.annot)
})



# Manually get colors -----------------------------------------------------





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

