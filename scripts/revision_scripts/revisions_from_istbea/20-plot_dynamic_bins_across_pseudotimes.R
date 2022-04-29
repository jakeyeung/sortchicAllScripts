# Jake Yeung
# Date of Creation: 2022-04-05
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/20-plot_dynamic_bins_across_pseudotimes.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)

# jmarksnew <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarksnew) <- jmarksnew
jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksold) <- jmarksold

# jmarknew <- "k4me1"
# jmarksold <- "H3K4me1"; names(jmarksold) <- jmarksold
# jmarkold <- "H3K4me1"

# Load meta ----------------------------------------------------------------

dat.meta.lst <- lapply(jmarksnew, function(jmark){
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths/metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.", jmark, ".txt")
  dat.meta <- fread(inf.meta)
})


# Load LDAs: all  ---------------------------------------------------------


tm.lst <- lapply(jmarks, function(jmark){
  if (jmark == "k4me1"){
    inf.ldaout <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.k4me1.2022-01-28/ldaOut.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.k4me1.2022-01-28.Robj"
  } else if (jmark == "k4me3"){
    inf.ldaout <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k4me3_cleaned/lda_outputs.count_mat_cleaned_no3no7_dynbins_allcells.k4me3.2022-04-12/ldaOut.count_mat_cleaned_no3no7_dynbins_allcells.k4me3.2022-04-12.Robj"
  } else if (jmark == "k27me3"){
    inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k27me3_clean_eryths/lda_outputs.count_mat_merged_with_old_dynbins.k27me3.2022-04-15/ldaOut.count_mat_merged_with_old_dynbins.", jmark, ".2022-04-15.Robj")
  } else if (jmark == "k9me3"){
    inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_repressive_cleaned_from_jupyter/ldaAnalysis_fripfilt_varfilt_binfilt/lda_outputs.count_mat_cleaned_dynbins.", jmark, ".2022-02-16/ldaOut.count_mat_cleaned_dynbins.", jmark, ".2022-02-16.Robj")
  }
  
  load(inf.ldaout, v=T)
  tm <- posterior(out.lda)
  return(tm)
})

print(lapply(tm.lst, function(x) dim(x$topics)))

dat.impute.lst <- lapply(tm.lst, function(tm){
  tm <- tm.lst[[jmarknew]]
  dat.impute <- log2(t(tm$topics %*% tm$terms))
})


# Load trajs --------------------------------------------------------------

indir.traj <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/cleaned"
dat.trajs <- lapply(jmarks, function(jmark){
  inf.trajs <- file.path(indir.traj, paste0("trajs_outputs.", jmark, ".rds"))
  readRDS(inf.trajs)
})


# # Load dynamic genes  -----------------------------------------------------
# 
# inf.genes <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/genesets/geneset_on_umap.binskeep_0.niter_500.2020-12-09.from_LDA_topics.condensed.heatmap.famousgenes.keepn_400.refmark_H3K4me3.2020-12-09.txt"
# dat.genes <- fread(inf.genes)
# 
# 
# # Show pseudotime on a gset specific genes  --------------------------------------------------------
# 
# jgenes <- subset(dat.genes, jset == "Granulocytes")$gene


# Plot on UMAP  -----------------------------------------------------------

# for TSS
# 
# # jgene <- "Cebpe$"
# jgene <- "Ebf1"
# jgene <- "Tal1"
# jgene <- "S100a8"
# # jrow <- jgenes[grepl(jgene, jgenes)]
# (jsub <- subset(dat.genes, grepl(jgene, gene)))
# jrow <- jsub$gene[[1]]
# 
# dat.exprs <- data.frame(cell = colnames(dat.impute), signal = dat.impute[jrow, ], stringsAsFactors = FALSE) %>%
#   left_join(dat.meta, .)
# 
# ggplot(dat.exprs, aes(x = umap1, y = umap2, color = signal)) + 
#   geom_point() + 
#   ggtitle(jrow) + 
#   scale_color_viridis_c() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Plot along traj ---------------------------------------------------------

# show granu traj
# dat.traj.sub <- subset(dat.trajs$Granulocytes)



# Get dynamic bins that are up in traj ------------------------------------

inf.dynbins <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/count_tables.BM.dynamic_bins_TSS_TES_regions/dynamic_bins.50kb.corrected_DE_tables.", jmarkold, ".2021-04-07.txt")
dat.de.bins <- fread(inf.dynbins)
bins.filt <- dat.de.bins$V4

# inf.rds <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/dendrograms/variance_by_clusters_by_HighBins.2021-03-09.1000000.log10filt_-1.methods.check_other_norms.add_DE_bins.rds"
# dat.mat.lst <- readRDS(inf.rds)[jmarksold]
# dat.mat <- dat.mat.lst[[jmarkold]]


dat.mat.log2.diff.lst <- lapply(dat.mat.lst, function(dat.mat){
  dat.mat.log2 <- log2(dat.mat)[bins.filt, ]
  dat.mat.log2.diff <- sweep(dat.mat.log2, MARGIN = 1, STATS = dat.mat.log2[, "HSPCs"], FUN = "-")
  cnames.keep <- colnames(dat.mat.log2.diff) != "HSPCs"
  return(dat.mat.log2.diff[, cnames.keep])
})

ctypes.k4 <- colnames(dat.mat.log2.diff.lst$H3K4me1)[colnames(dat.mat.log2.diff.lst$H3K4me1) != "HSPCs"]
names(ctypes.k4) <- ctypes.k4

ctypes.k9 <- colnames(dat.mat.log2.diff.lst$H3K9me3)[colnames(dat.mat.log2.diff.lst$H3K9me3) != "HSPCs"]
names(ctypes.k9) <- ctypes.k9


ctype.spec.up <- lapply(jmarksold, function(jmark){
  if (jmark != "H3K9me3"){
    ctypes <- ctypes.k4
  } else {
    ctypes <- ctypes.k9
  }
  print(jmark)
  print(ctypes)
  lapply(ctypes, function(ctype){
    cnames.keep <- colnames(dat.mat.log2.diff.lst[[jmark]]) != ctype
    jmat <- dat.mat.log2.diff.lst[[jmark]][, cnames.keep]
    jref <- dat.mat.log2.diff.lst[[jmark]][, ctype]
    ctype.up.any <- jref > 0
    others.down.all <- apply(jmat, MARGIN = 1, FUN = function(jrow) all(jrow < 0))
    ctype.up.vec <- ctype.up.any & others.down.all
    return(names(which(ctype.up.vec)))
  })
})

ctype.spec.down <- lapply(jmarksold, function(jmark){
  if (jmark != "H3K9me3"){
    ctypes <- ctypes.k4
  } else {
    ctypes <- ctypes.k9
  }
  print(jmark)
  lapply(ctypes, function(ctype){
    cnames.keep <- colnames(dat.mat.log2.diff.lst[[jmark]]) != ctype
    jmat <- dat.mat.log2.diff.lst[[jmark]][, cnames.keep]
    jref <- dat.mat.log2.diff.lst[[jmark]][, ctype]
    ctype.down.any <- jref < 0
    others.up.all <- apply(jmat, MARGIN = 1, FUN = function(jrow) all(jrow > 0))
    ctype.down.vec <- ctype.down.any & others.up.all
    return(names(which(ctype.down.vec)))
  })
})

# check 

bins.check <- ctype.spec.down$H3K4me1$Granulocytes
bins.check <- ctype.spec.up$H3K4me1$Bcells
bins.check <- ctype.spec.up$H3K4me1$Eryths

rows.keep <- which(rownames(dat.impute) %in% bins.check)
bins.notin <- bins.check[which(!bins.check %in% rownames(dat.impute))]
print(length(rows.keep))
print(length(bins.notin))

dat.signal <- data.frame(cell = colnames(dat.impute), signal = colMeans(dat.impute[rows.keep, ]), stringsAsFactors = FALSE)
dat.exprs <- dat.signal %>%
  left_join(dat.meta, .)

ggplot(dat.exprs, aes(x = umap1, y = umap2, color = signal)) + 
  geom_point() + 
  scale_color_viridis_c() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Plot those along pseudotime along all trajs -----------------------------

jctypes <- names(dat.trajs[[jmarknew]]); names(jctypes) <- jctypes

dat.trajs.long <- lapply(jctypes, function(jctype){
  dat.trajs.sub <- dat.trajs[[jmarknew]][[jctype]] %>%
    filter(is.ctype) %>%
    dplyr::select(cell, ctype.ordered, ptime) %>%
    mutate(traj = jctype)
}) %>%
  bind_rows() %>%
  left_join(., dat.signal)

# add exprs
ggplot(dat.trajs.long, aes(x = ptime, y = signal, color = traj)) + 
  geom_point() + 
  facet_wrap(~traj) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


