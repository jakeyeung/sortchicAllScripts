# Jake Yeung
# Date of Creation: 2020-01-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/intestinal_all_merged/2-make_clustering.R
# Make clustering table so we can make bigwigs out of bams (pseudobulks)

rm(list=ls())

library(scchicFuncs)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(hash)
library(igraph)
library(umap)

library(topicmodels)

library(leidenbase)
library(jclustering)


# Settings ----------------------------------------------------------------


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")


# Analyze for 50kb 25kb ---------------------------------------------------

jmarks <- c("k36me3", "k4me1", "k4me3", "k27me3", "k9me3")

clstr.meths <- c("leiden", "louvain")
# clstr.meth <- "louvain"

for (clstr.meth in clstr.meths){
  
  # for (jmark in jmarks){
  parallel::mclapply(jmarks, function(jmark){
    
    # for (jmark in jmarks){
    inf <- paste0("/home/jyeung/hpc/intestinal_scchic/LDA_outputs/topicmodels/ldaAnalysisBins_intestines.2019-12-22/lda_outputs.mat.Scraped.AllMerged.", jmark, ".TAcutoff_0.5.countscutoff_500_1000.binfilt_cellfilt.2019-12-22.K-30.binarize.FALSE/ldaOut.mat.Scraped.AllMerged.", jmark, ".TAcutoff_0.5.countscutoff_500_1000.binfilt_cellfilt.2019-12-22.K-30.Robj")
    assertthat::assert_that(file.exists(inf))
    print(inf)
    
    outpdf <- paste0("/home/jyeung/hpc/intestinal_scchic/from_rstudioiserver/pdfs_clustering/plots_intestine_AllMerged_", jmark, ".clstmeth_", clstr.meth, ".pdf")
    outtxt <- paste0("/home/jyeung/hpc/intestinal_scchic/from_rstudioiserver/pdfs_clustering/clusters_intestine_AllMerged_", jmark, ".clstmeth_", clstr.meth, ".txt")
    
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
    
    if (clstr.meth == "louvain"){
      dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long, clstr.cname = "cluster")  # make clstr.cname match DoLeiden()
    } else if (clstr.meth == "leiden"){
      dat.umap.long <- DoLeiden(topics.mat, dat.umap.long = dat.umap.long, K=jsettings$n_neighbors, res_param=10^seq(-5, 0), random_seed = 123, weight=FALSE, verbose = TRUE)
    } else {
      stop(paste("clstr.methd must be louv or leiden: ", clstr.meth))
    }
    
    ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    # color by experiment 
    dat.umap.long <- dat.umap.long %>%
      rowwise() %>%
      mutate(experi = ClipLast(cell, jsep = "-"),
             plate = ClipLast(cell, jsep = "_"),
             prefix = paste(strsplit(cell, split = "-")[[1]][1:3], collapse = "-"))
    
    # dat.umap.long <- DoLouvain(topics.mat, custom.settings.louv = jsettings, dat.umap.long)
    
    m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
      scale_color_viridis_d() + 
      facet_wrap(~prefix) + 
      ggtitle(jmark)
    print(m)
    
    m.together <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
      scale_color_viridis_d() + 
      ggtitle(jmark)
    print(m.together)
    
    m.louv <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cluster)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
      scale_color_viridis_d() + 
      facet_wrap(~prefix) +
      ggtitle(jmark)
    print(m.louv)
    
    m.louv.together <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cluster)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
      scale_color_viridis_d() + 
      ggtitle(jmark)
    print(m.louv.together)
    
    dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
    jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
    dat.var <- CalculateVarAll(dat.impute.log, jchromos)
    dat.merge <- left_join(dat.umap.long, dat.var)
    
    m.var <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
      scale_color_viridis_c(direction = -1) + 
      facet_wrap(~prefix) +
      ggtitle(jmark)
    
    m.var.together <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
      scale_color_viridis_c(direction = -1) + 
      ggtitle(jmark)
    
    m.var.plate <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
      scale_color_viridis_c(direction = -1) + 
      facet_wrap(~plate) +
      ggtitle(jmark)
    
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
    out.dat <- subset(dat.merge, select = c(cell, cluster, umap1, umap2, cell.var.within.sum.norm))
    fwrite(x = out.dat, file = outtxt, col.names = TRUE, sep = "\t")
  },
  mc.cores = length(jmarks))
  
}







  
