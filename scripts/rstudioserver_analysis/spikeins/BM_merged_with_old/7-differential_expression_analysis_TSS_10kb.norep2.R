# Jake Yeung
# Date of Creation: 2020-11-14
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/7-differential_expression_analysis_TSS_10kb.norep2.R
# 


rm(list=ls())

library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)
library(hash)
library(igraph)
library(umap)

library(scchicFuncs)

jstart <- Sys.time()



jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Load LDA (contains countmat)  ---------------------------------------------------------------

ncores <- 8
hubprefix <- "/home/jyeung/hub_oudenaarden"
jdist <- 10000

# outdir <- "/home/jyeung/data/from_rstudioserver/spikein_fits_BM_poisson"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3"
dir.create(outdir)

# jmark <- "H3K4me1"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

for (jmark in jmarks){
  
  indir <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.from_TSS.dist_", jdist))
  
  fname <- paste0("lda_outputs.count_mat_from_TSS.", jmark, ".dist_", jdist, ".K-30.binarize.FALSE/ldaOut.count_mat_from_TSS.", jmark, ".dist_", jdist, ".K-30.Robj")
  
  load(file.path(indir, fname), v=T)
  
  tm.result <- posterior(out.lda)
  
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)
  
  ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # Load metadata -----------------------------------------------------------
  
  # indir.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2")
  # inf.annot <- file.path(indir.annot, paste0("spikeins_mouse.BMround1and2_umaps_and_ratios.colfix.celltyping.2020-11-01.WithRelLevels.mark_", jmark, ".cell_cluster_tables.txt"))
  indir.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs.merge_with_new_BM/clusterfilt.2020-11-04")
  inf.annot <- file.path(indir.annot, paste0("cell_cluster_table.old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.txt"))
  dat.annot <- fread(inf.annot)
  
  cells.keep <- colnames(count.mat)[!grepl("rep2", colnames(count.mat))]
  dat.annot.filt <- subset(dat.annot, cell %in% cells.keep)
  
  dat.umap.merge <- left_join(dat.umap, subset(dat.annot, select = c(cell, cluster)))
  
  cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
  ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() +  
    scale_color_manual(values = cbPalette) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
  # Load spikeins -----------------------------------------------------------
  
  
  
  
  # Run fits gene by gene ---------------------------------------------------
  
  
  print(jmark)
  jmat.mark <- count.mat
  dat.annots.filt.mark <- dat.annot.filt %>%
    mutate(Cluster = ifelse(cluster == "HSPCs", "aHSPCs", cluster)) %>%
    rowwise() %>%
    mutate(batch = IsRound1(cell, mark = jmark)) %>%
    mutate(Plate = ifelse(batch == "Round2", plate, "Round1"))
  
  print(unique(dat.annots.filt.mark$Plate))
  
  ncuts.for.fit.mark <- data.frame(cell = colnames(count.mat), ncuts.total = colSums(count.mat), stringsAsFactors = FALSE)
  cnames <- colnames(jmat.mark)
  
  jrow.names <- rownames(jmat.mark)
  names(jrow.names) <- jrow.names
  
  
  jfits.lst <- parallel::mclapply(jrow.names, function(jrow.name){
    jrow <- jmat.mark[jrow.name, ]
    jout <- FitGlmRowClustersPlate(jrow, cnames, dat.annots.filt.mark, ncuts.for.fit.mark, jbin = NULL, returnobj = FALSE)
    return(jout)
  }, mc.cores = ncores)
  
  
  # Ssave outputs -----------------------------------------------------------
  outf <- file.path(outdir, paste0("poisson_fit_TSS_", jdist, ".",  jmark, ".", Sys.Date(), ".newannot2.norep2.RData"))
  # saveRDS(jfits.lst, outf)
  save(jfits.lst, dat.annots.filt.mark, jmat.mark, ncuts.for.fit.mark, file = outf)
  
  print(Sys.time() - jstart)
}

print("Done")
print(Sys.time() - jstart)








