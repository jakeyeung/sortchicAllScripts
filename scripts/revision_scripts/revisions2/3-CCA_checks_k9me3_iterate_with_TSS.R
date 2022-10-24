# Jake Yeung
# Date of Creation: 2022-07-22
# File: ~/projects/scchic/scripts/revision_scripts/revisions2/3-CCA_checks_k9me3.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksold) <- jmarks

# Metas -------------------------------------------------------------------


inf.colors.fixed <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/ctypes_on_umap_batch_corrected_colors_fixed/dat_colors_DC_monocyte_fixed.2022-05-17.txt"
dat.colors.fixed <- fread(inf.colors.fixed)



dat.meta.lst <- lapply(jmarks, function(jmark){
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umaps_pcas_with_batch_corrections/umap_metadata_primetime.", jmark, ".2022-04-21.txt")
  dat.meta <- fread(inf.meta) %>%
    left_join(., dat.colors.fixed) %>%
    rowwise() %>%
    mutate(colcode = colcodenew)
  # replace colcode with colcodenew
})

dat.meta.colors <- subset(dat.meta.lst$k4me1, select = c(ctype.from.LL, colcode))
ctype2col <- hash::hash(dat.meta.colors$ctype.from.LL, dat.meta.colors$colcode)

dat.meta.merge <- lapply(jmarks, function(jmark){
  jdat <- dat.meta.lst[[jmark]]
  subset(jdat, select = c(cell, ctype.from.LL, colcode)) %>%
    rowwise() %>%
    mutate(mark = jmark)
}) %>%
  bind_rows()

# Check CCA outliers?  ----------------------------------------------------

# inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/with_k9me3_cluster_specific_bins/UMAP_of_CCA_repressive_k9me3_and_k4me1.2022-07-22.RData"
jmarkref <- "k9me3"
jmarkothers <- jmarks[jmarks != jmarkref]
dat.umap.cca.lst <- lapply(jmarkothers, function(jmarkother){
  print(jmarkother)
  # inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/with_k9me3_cluster_specific_bins_keeptop_500/UMAP_of_CCA_repressive_", jmarkref, "_and_", jmarkother, ".2022-07-22.RData")
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/with_k9me3_cluster_specific_bins_keeptop_500_bymark_factor_-1/UMAP_of_CCA_repressive_", jmarkref, "_and_", jmarkother, ".2022-07-22.RData")
  load(inf.meta, v=T)
  dat.umap.cca <- dat.umap.cca %>%
    left_join(., dat.meta.merge)
})

m.lst <- lapply(jmarkothers, function(jmarkother){
  ggplot(dat.umap.cca.lst[[jmarkother]], aes(x = umap1, y = umap2, color = colcode)) + 
    geom_point() + 
    ggtitle(paste(jmarkref, "in", jmarkother)) + 
    scale_color_identity() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
JFuncs::multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], cols = 3)

m.mark.lst <- lapply(jmarkothers, function(jmarkother){
  ggplot(dat.umap.cca.lst[[jmarkother]], aes(x = umap1, y = umap2, color = mark)) + 
    geom_point() + 
    ggtitle(paste(jmarkref, "in", jmarkother)) + 
    facet_wrap(~ctype.from.LL) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
print(m.mark.lst)
# JFuncs::multiplot(m.mark.lst[[1]], m.mark.lst[[2]], m.mark.lst[[3]])



# Get signal around TSS  --------------------------------------------------

jmarkold2new <- hash::hash(jmarksold, names(jmarksold))

counts.tss.lst <- lapply(jmarksold, function(jmark){
  print(jmark)
  # inf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_TSS_k9_cluster_specific_bins_keeptop_500/lda_outputs.count_mat_allmerged_for_LDA_k9_cluster_specific_bins_keeptop_500.", jmark, ".2022-07-22/ldaOut.count_mat_allmerged_for_LDA_k9_cluster_specific_bins_keeptop_500.", jmark, ".2022-07-22.Robj")
  # inf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_TSS_only/lda_outputs.count_mat_allmerged_for_LDA_TSS_only.", jmark, ".2022-07-20/ldaOut.count_mat_allmerged_for_LDA_TSS_only.", jmark, ".2022-07-20.Robj")
  inf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/count_mats_TSS/counts_tables_50000/BM_", jmark, "/BM_allmerged_", jmark, ".countTable.binsize_50000.csv")
  print(inf)
  count.mat <- ReadMatTSSFormat(inf, as.sparse = TRUE, add.coord = TRUE, sort.rnames = TRUE)
  # filter good cells
  cells.keep <- dat.meta.lst[[jmarkold2new[[jmark]]]]$cell
  print(paste("Cells keep:", length(cells.keep)))
  print("Dim before...")
  print(dim(count.mat))
  # cols.keep <- colnames(count.mat) %in% cells.keep
  count.mat <- count.mat[, cells.keep]
  print("Dim after...")
  print(dim(count.mat))
  return(count.mat)
})

# summarize tss 
summary.tss.lst <- lapply(jmarks, function(jmark){
  cuts.total <- colSums(counts.tss.lst[[jmark]])
  data.frame(samp = names(cuts.total), cuts.total = cuts.total, jtype = "TSS", nbins = nrow(counts.tss.lst[[jmark]]), stringsAsFactors = FALSE)
})




# Get signal around interesting bins ?  -----------------------------------

# different dynamic bins
indir.bins <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/dynamic_bins"
dat.bins.lst <- lapply(jmarks, function(jmark){
  inf <- file.path(indir.bins, paste0("dynamic_bins_50kb.", jmark, ".2022-07-24.txt"))
  bed.bins <- fread(inf)
  return(bed.bins$V4)
})


summary.dynbins.bybin.lst <- lapply(jmarksold, function(jmark){
  print(jmark)
  inmain <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/count_mats_bins/counts_tables_50000"
  indir <- file.path(inmain, paste0("BM_", jmark))
  inf <- file.path(indir, paste0("BM_allmerged_", jmark, ".countTable.binsize_50000.csv"))
  print(inf)
  count.mat <- ReadMatSlideWinFormat(inf, as.sparse = TRUE, sort.rnames = TRUE)
  # filter good cells
  jmarknew <- jmarkold2new[[jmark]]
  cells.keep <- dat.meta.lst[[jmarknew]]$cell
  print(paste("Cells keep:", length(cells.keep)))
  print("Dim before...")
  print(dim(count.mat))
  # cols.keep <- colnames(count.mat) %in% cells.keep
  count.mat <- count.mat[, cells.keep]
  print("Dim after...")
  print(dim(count.mat))
  # summarize: four dynamic bins per mat
  print("Getting four dynamic bins:")
  summary.dynbins.lst <- lapply(jmarks, function(jmarknew){
    print(jmarknew)
    dyn.bins <- dat.bins.lst[[jmarknew]]
    count.mat.filt <- count.mat[dyn.bins, ]
    cuts.total <- colSums(count.mat.filt)
    dat.summary <- data.frame(samp = names(cuts.total), 
                              cuts.total = cuts.total, 
                              jtype = jmarknew, 
                              nbins = nrow(count.mat.filt), 
                              stringsAsFactors = FALSE)
  })
  return(summary.dynbins.lst)
})

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/counts_different_bins"
save(summary.tss.lst, summary.dynbins.bybin.lst, file = file.path(outdir, "summary_bins_objs.", Sys.Date(), ".RData"))

# 
# # Sum across marks, normalize by depth  -----------------------------------
# 
# # get counts in k9me3 bins 
# dat.rz.with.counts <- lapply(jmarks, function(jmark){
#   print(jmark)
#   dat.rz.tmp <- dat.rz.lst[[jmark]] %>%
#     left_join(., cuts.at.k9.lst[[jmark]], by = "samp") %>%
#     left_join(., cuts.at.tss.lst[[jmark]], by = "samp") 
#     filter(!is.na(cut.k9)) %>%
#     rowwise() %>%
#     mutate(cut.k9.per.bin = cut.k9 / nbins.k9, 
#           cut.tss.per.bin = cut.tss / nbins.tss)
# })
# 
# 
# #   
# # # Can we get UMAP from k9me3 LDA ------------------------------------------
# # 
# # # inf.k9 <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_TSS_k9_cluster_specific_bins/lda_outputs.count_mat_allmerged_for_LDA_k9_cluster_specific_bins.H3K9me3.2022-07-21/ldaOut.count_mat_allmerged_for_LDA_k9_cluster_specific_bins.H3K9me3.2022-07-21.Robj"
# # inf.k9 <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_TSS_k9_cluster_specific_bins_keeptop_500/lda_outputs.count_mat_allmerged_for_LDA_k9_cluster_specific_bins_keeptop_500.H3K9me3.2022-07-22/ldaOut.count_mat_allmerged_for_LDA_k9_cluster_specific_bins_keeptop_500.H3K9me3.2022-07-22.Robj"
# # load(inf.k9, v=T)
# # 
# # tm <- out.lda %>%
# #   posterior() %>%
# #   scchicFuncs::AddTopicToTmResult()
# # 
# # 
# # library(hash)
# # library(igraph)
# # library(umap)
# # jsettings <- umap.defaults
# # jsettings[["n_neighbors"]] <- 30
# # jsettings[["min_dist"]] <- 0.1
# # jsettings[["random_state"]] <- 123
# # 
# # dat.umap.pca <- DoUmapAndLouvain(tm$topics, jsettings)
# # dat.umap.pca.annot <- dat.umap.pca %>%
# #   left_join(., dat.meta.merge)
# #   
# # ggplot(dat.umap.pca.annot, aes(x = umap1, y = umap2, color = colcode))  + 
# #   geom_point() + 
# #   theme_bw() + 
# #   scale_color_identity() + 
# #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# # 
