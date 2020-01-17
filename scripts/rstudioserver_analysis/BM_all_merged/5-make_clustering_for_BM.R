# Jake Yeung
# Date of Creation: 2020-01-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/intestinal_all_merged/2-make_clustering.R
# Make clustering table so we can make bigwigs out of bams (pseudobulks)

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(hash)
library(igraph)
library(umap)

library(leiden)
library(topicmodels)

library(RColorBrewer)


# Functions ---------------------------------------------------------------

DoLeiden <- function(topics.mat, custom.settings, dat.umap.long = NULL, jres.param = 1){
  # Do Louvain for clustering
  
  dat.umap <- umap(topics.mat, config = custom.settings)
  dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout),
                                   stringsAsFactors = FALSE)
  cell.indx <- hash::hash(rownames(dat.umap$knn$indexes), dat.umap$knn$indexes[, 1])
  cell.indx.rev <- hash(dat.umap$knn$indexes[, 1], rownames(dat.umap$knn$indexes))
  nr <- nrow(dat.umap$knn$indexes)
  nc <- ncol(dat.umap$knn$indexes)
  edgelist <- matrix(NA, nrow = nr * nc, ncol = 2)
  colnames(edgelist) <- c("from", "to")
  for (vertex.i in seq(nr)){
    istart <- nc*(vertex.i - 1)+1
    iend <- nc*vertex.i
    edgelist[istart : iend, 1] <- cell.indx.rev[[as.character(vertex.i)]]
    edgelist[istart : iend, 2] <- sapply(dat.umap$knn$indexes[vertex.i, 1:nc], function(x) cell.indx.rev[[as.character(x)]])
    # edgelist[istart : iend, 3] <- 1 / (dat.umap$knn$distances[vertex.i, 1:nc] + 0.1)
  }
  g <- graph_from_data_frame(edgelist, directed=FALSE)
  g.out <- leiden(g, resolution_parameter = jres.param)
  V(g)$color <- g.out$membership
  clstr <- hash(g.out$names, g.out$membership)
  if (is.data.frame(dat.umap.long)){
    dat.umap.long$louvain <- as.character(sapply(dat.umap.long$cell, function(x) clstr[[as.character(x)]]))
  } else {
    dat.umap.long <- clstr
  }
  return(dat.umap.long)
  
}


# Settings ----------------------------------------------------------------


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")


# Analyze for 50kb 25kb ---------------------------------------------------

# jwins <- "20000"
jwins <- c("20000_10000", "50000_25000", "100000_20000")
names(jwins) <- jwins
# jmark <- "H3K9me3"
# jmarks <- c("H3K27me3", "H3K9me3")
# jwin <- "50000_25000"
jmarks <- c("H3K4me1")

outdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/pdfs_clustering_BM"

for (jmark in jmarks){
  
  infs <- list(paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2019-12-30.bsizestepsize_20000_10000/lda_outputs.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj"),
               paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2019-12-30.bsizestepsize_50000_25000/lda_outputs.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj"),
               paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2019-12-16/lda_outputs.B6BM_AllMerged_", jmark, ".TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-16.K-50.binarize.FALSE/ldaOut.B6BM_AllMerged_", jmark, ".TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-16.K-50.Robj"))
  names(infs) <- jwins
  
  for (jwin in jwins){
    inf <- infs[[jwin]]
    # inf <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2019-12-30.bsizestepsize_", jwin, "/lda_outputs.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj")
    outpdf <- file.path(outdir, paste0("plots_BM_All_merged_", jmark, "_", jwin, ".pdf"))
    outtxt <- file.path(outdir, paste0("clusters_BM_All_merged_", jmark, "_", jwin, ".txt"))
    if (file.exists(outpdf) & file.exists(outtxt)){
      print(outpdf)
      print(outtxt)
      print("Both output files exist, skipping...")
      next
    }
    assertthat::assert_that(file.exists(inf))
    load(inf, v=T)
    print(out.lda@terms[1:5])
    
    tm.result <- posterior(out.lda)
    topics.mat <- posterior(out.lda)$topics
    
    umap.out <- umap(topics.mat, config = jsettings)
    dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
    dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
    ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    # color by experiment 
    dat.umap.long <- dat.umap.long %>%
      rowwise() %>%
      mutate(experi = ClipLast(cell, jsep = "-"),
             plate = ClipLast(cell, jsep = "_"),
             prefix = paste(strsplit(cell, split = "-")[[1]][1:3], collapse = "-"))
    
    dat.umap.long <- DoLouvain(topics.mat, custom.settings.louv = jsettings, dat.umap.long)
    
    m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
      scale_color_manual(values = cbPalette) + 
      facet_wrap(~prefix) + 
      ggtitle("WinSize:", jwin)
    print(m)
    
    m.together <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
      scale_color_manual(values = cbPalette) + 
      ggtitle("WinSize:", jwin)
    print(m.together)
    
    m.louv <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
      scale_color_manual(values = cbPalette) + 
      facet_wrap(~prefix) +
      ggtitle("WinSize:", jwin)
    print(m.louv)
    
    m.louv.together <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
      scale_color_manual(values = cbPalette) + 
      ggtitle("WinSize:", jwin)
    print(m.louv.together)
    
    dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
    jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
    dat.var <- CalculateVarAll(dat.impute.log, jchromos)
    dat.merge <- left_join(dat.umap.long, dat.var)
    
    m.var <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
      scale_color_viridis_c(direction = -1) + 
      facet_wrap(~prefix) +
      ggtitle("WinSize:", jwin)
    
    m.var.together <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
      scale_color_viridis_c(direction = -1) + 
      ggtitle("WinSize:", jwin)
    
    m.var.plate <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
      scale_color_viridis_c(direction = -1) + 
      facet_wrap(~plate) +
      ggtitle("WinSize:", jwin)
    
    # write table output with pdf to accompany
    pdf(file = outpdf, useDingbats = FALSE)
    print(m)
    print(m.together)
    print(m.louv)
    print(m.louv.together)
    print(m.var)
    print(m.var.together)
    print(m.var.plate)
    dev.off()
    
    # write textfile, first column is cellname, second column is cluster name, other columns are metadata. First row is colnames
    out.dat <- subset(dat.merge, select = c(cell, louvain, umap1, umap2, cell.var.within.sum.norm)) %>%
      mutate(louvain = paste("louvain", louvain, sep = "_"))
    fwrite(x = out.dat, file = outtxt, col.names = TRUE, sep = "\t")
  }
  
}





  