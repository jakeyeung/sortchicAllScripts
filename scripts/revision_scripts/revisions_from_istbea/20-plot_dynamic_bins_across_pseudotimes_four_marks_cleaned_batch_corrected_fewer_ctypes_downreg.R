# Jake Yeung
# Date of Creation: 2022-04-05
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/20-plot_dynamic_bins_across_pseudotimes.R
# 


rm(list=ls())

library(scchicFuncs)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)

jratio <- 0.66

# jmarksnew <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarksnew) <- jmarksnew
jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksold) <- jmarks

# jmarknew <- "k4me1"
# jmarksold <- "H3K4me1"; names(jmarksold) <- jmarksold
# jmarkold <- "H3K4me1"

# Load meta ----------------------------------------------------------------

dat.meta.lst <- lapply(jmarks, function(jmark){
  # inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths2/metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.", jmark, ".txt")
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umaps_pcas_with_batch_corrections/umap_metadata_primetime.", jmark, ".2022-04-21.txt")
  dat.meta <- fread(inf.meta)
})

dat.meta.colors <- subset(dat.meta.lst$k4me1, select = c(ctype.from.LL, colcode))
ctype2col <- hash::hash(dat.meta.colors$ctype.from.LL, dat.meta.colors$colcode)

# Load LDAs: all  ---------------------------------------------------------


tm.lst <- lapply(jmarks, function(jmark){
  if (jmark == "k4me1"){
    inf.ldaout <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.k4me1.2022-01-28/ldaOut.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.k4me1.2022-01-28.Robj"
  } else if (jmark == "k4me3"){
    inf.ldaout <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k4me3_cleaned/lda_outputs.count_mat_cleaned_no3no7_dynbins_allcells.k4me3.2022-04-12/ldaOut.count_mat_cleaned_no3no7_dynbins_allcells.k4me3.2022-04-12.Robj"
  } else if (jmark == "k27me3"){
    # inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k27me3_clean_eryths/lda_outputs.count_mat_merged_with_old_dynbins.k27me3.2022-04-15/ldaOut.count_mat_merged_with_old_dynbins.", jmark, ".2022-04-15.Robj")
    # inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k27me3_clean_eryths/lda_outputs.count_mat_new_only_dynbins.k27me3.2022-04-15/ldaOut.count_mat_new_only.", jmark, ".2022-04-15.Robj")
    inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k27me3_clean_eryths/lda_outputs.count_mat_new_only_dynbins.k27me3.2022-04-15/ldaOut.count_mat_new_only_dynbins.", jmark, ".2022-04-15.Robj")
  } else if (jmark == "k9me3"){
    inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_repressive_cleaned_from_jupyter/ldaAnalysis_fripfilt_varfilt_binfilt/lda_outputs.count_mat_cleaned_dynbins.", jmark, ".2022-02-16/ldaOut.count_mat_cleaned_dynbins.", jmark, ".2022-02-16.Robj")
  }
  
  load(inf.ldaout, v=T)
  tm <- posterior(out.lda)
  return(tm)
})

# assertthat::assert_that(nrow(dat.meta.lst$k27me3) == nrow(tm.lst$k27me3$topics))

print(lapply(tm.lst, function(x) dim(x$topics)))

dat.impute.lst <- lapply(tm.lst, function(tm){
  dat.impute <- log2(t(tm$topics %*% tm$terms))
})


# Load batch corrected k27me3 ---------------------------------------------------

inf.impute.k27me3.bc <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/mat_wide_k27me3_batch_corrected.2022-04-19.rds"
dat.impute.k27me3.bc <- readRDS(inf.impute.k27me3.bc)
dat.impute.lst$k27me3 <- dat.impute.k27me3.bc

# Load trajs --------------------------------------------------------------

# indir.traj <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/cleaned2_batch_corrected_eryth_fix"
indir.traj <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/cleaned"
dat.trajs <- lapply(jmarks, function(jmark){
  inf.trajs <- file.path(indir.traj, paste0("trajs_outputs.", jmark, ".rds"))
  if (jmark == "k27me3"){
    inf.trajs <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/cleaned2_batch_corrected_eryth_fix/trajs_outputs_batch_corrected.k27me3.2022-04-21.rds")
  }
  print(inf.trajs)
  readRDS(inf.trajs)
})


# Get dynamic bins that are up in traj ------------------------------------

bins.filt.lst <- lapply(jmarksold, function(jmarkold){
  inf.dynbins <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/count_tables.BM.dynamic_bins_TSS_TES_regions/dynamic_bins.50kb.corrected_DE_tables.", jmarkold, ".2021-04-07.txt")
  dat.de.bins <- fread(inf.dynbins)
  bins.filt <- dat.de.bins$V4
  return(bins.filt)
})


# Mean across mat ---------------------------------------------------------


dat.pbulk.lst <- lapply(jmarks, function(jmark){
  dat.impute <- dat.impute.lst[[jmark]]
  dat.meta <- dat.meta.lst[[jmark]]
  cnames.keep.lst <- lapply(split(x = dat.meta, f = dat.meta$ctype.from.LL), function(jlst){
    jlst$cell
  })
  mean.lst <- ApplyAcrossClusters(count.mat = dat.impute, cnames.keep.lst = cnames.keep.lst, fn = mean)
  mat.mean <- do.call(cbind, mean.lst)
})

# Get diff relative to HSCs -----------------------------------------------

dat.mat.log2.diff.lst <- lapply(jmarks, function(jmark){
  dat.mat.log2 <- dat.pbulk.lst[[jmark]]
  bins.filt <- bins.filt.lst[[jmark]]
  dat.mat.log2.diff <- sweep(dat.mat.log2, MARGIN = 1, STATS = dat.mat.log2[, "HSCs"], FUN = "-")
  cnames.keep <- colnames(dat.mat.log2.diff) != "HSCs"
  return(dat.mat.log2.diff[, cnames.keep])
})

ctypes.lst <- lapply(dat.mat.log2.diff.lst, function(jdat){
  jnames <- colnames(jdat); names(jnames) <- jnames
  return(jnames)
})

ctype.spec.up <- lapply(jmarks, function(jmark){
  print(jmark)
  ctypes <- ctypes.lst[[jmark]]
  # print(ctypes)
  lapply(ctypes, function(ctype){
    # cnames.keep <- colnames(dat.mat.log2.diff.lst[[jmark]]) != ctype
    # jmat 
    # jmat <- dat.mat.log2.diff.lst[[jmark]][, cnames.keep]
    jref <- dat.mat.log2.diff.lst[[jmark]][, ctype]
    ctype.up.any <- jref > 0
    # others.down.all <- apply(jmat, MARGIN = 1, FUN = function(jrow) all(jrow < 0))
    # ctype.up.vec <- ctype.up.any & others.down.all
    return(names(which(ctype.up.any)))
  })
})

ctype.spec.down <- lapply(jmarks, function(jmark){
  print(jmark)
  ctypes <- ctypes.lst[[jmark]]
  lapply(ctypes, function(ctype){
    # cnames.keep <- colnames(dat.mat.log2.diff.lst[[jmark]]) != ctype
    # jmat <- dat.mat.log2.diff.lst[[jmark]][, cnames.keep]
    jref <- dat.mat.log2.diff.lst[[jmark]][, ctype]
    ctype.down.any <- jref < 0
    # others.up.all <- apply(jmat, MARGIN = 1, FUN = function(jrow) all(jrow > 0))
    # ctype.down.vec <- ctype.down.any & others.up.all
    return(names(which(ctype.down.any)))
  })
})

# check downs
cnames.keep <- which(colnames(dat.mat.log2.diff.lst$k4me1) %in% c("Granulocytes", "Eryths", "Bcells"))
dat.mat.check <- dat.mat.log2.diff.lst$k4me1[, c("Granulocytes", "Eryths", "Bcells")]
jbins <- ctype.spec.up$k4me1$Granulocytes
jbins <- ctype.spec.up$k4me1$Bcells
jbins <- ctype.spec.down$k4me1$Bcells
jbins <- ctype.spec.down$k4me1$Granulocytes
# jbins <- ctype.spec.down$k4me1$Granulocytes
jfilt <- as.data.frame(dat.mat.check[jbins, ])

colMeans(jfilt)
apply(jfilt, 2, median)

plot(density(jfilt[, "Eryths"]))
plot(density(jfilt[, "Bcells"]))
plot(density(jfilt[, "Granulocytes"]))

jfilt.ordered <- head(as.data.frame(jfilt) %>% arrange(desc(Eryths)))
jfilt.ordered <- head(as.data.frame(jfilt) %>% arrange(Eryths))
jfilt.ordered <- jfilt %>% filter(Granulocytes > -1.5)
colMeans(jfilt.ordered)
apply(jfilt.ordered, 2, median)

# check 
jmarktest <- "k9me3"
jmarktest <- "k27me3"
jctypetest <- "Granulocytes"

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umap_and_trajs/k27me3_batch_corrected_fewer_ctypes_downreg"
dir.create(outdir)
for (jmark in jmarks){
  print(jmark)
  
  
  # jctypes <- names(ctype.spec.up[[jmark]])
  jctypes <- c("Granulocytes", "Bcells", "Eryths")
  names(jctypes) <- jctypes
  for (jctype in jctypes){
    
    pdfout <- file.path(outdir, paste0("umap_on_trajs.", jmark, ".", jctype, ".", Sys.Date(), ".pdf"))
    pdf(pdfout, useDingbats = FALSE)
    print(jctype)
    bins.check <- ctype.spec.up[[jmark]][[jctype]]
    bins.check <- ctype.spec.down[[jmark]][[jctype]]
    
    # bins.check <- rownames(jfilt)
    
    # bins.check <- sample(bins.check, size = 10, replace = FALSE)
    # bins.check <- rownames(jfilt.ordered)
    
    # bins.check <- ctype.spec.up$k27me3$Granulocytes
    print("Nbins")
    print(length(bins.check))
    
    rows.keep <- which(rownames(dat.impute.lst[[jmark]]) %in% bins.check)
    # bins.notin <- bins.check[which(!bins.check %in% rownames(dat.impute))]
    print(length(rows.keep))
    jtitle <- paste(jmark, jctype, "Nbins:", length(rows.keep))
    # print(length(bins.notin))
    
    dat.signal <- data.frame(cell = colnames(dat.impute.lst[[jmark]]), signal = colMeans(dat.impute.lst[[jmark]][rows.keep, ]), stringsAsFactors = FALSE)
    # dat.signal <- data.frame(cell = colnames(dat.impute.lst[[jmark]]), 
    #                          signal = apply(dat.impute.lst[[jmark]][rows.keep, ], MARGIN = 2, median), 
    #                          stringsAsFactors = FALSE)
    dat.exprs <- dat.signal %>%
      left_join(dat.meta.lst[[jmark]], .)
    
    m.signal <- ggplot(dat.exprs, aes(x = umap1, y = umap2, color = signal)) + 
      geom_point() + 
      scale_color_viridis_c() + 
      ggtitle(jtitle) + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    m.signal.batch <- ggplot(dat.exprs, aes(x = umap1, y = umap2, color = signal)) + 
      geom_point() + 
      scale_color_viridis_c() + 
      facet_wrap(~batch) + 
      ggtitle(jtitle) + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    print(m.signal)
    print(m.signal.batch)
    
    # show trajectory
    
    traj.ctypes <- names(dat.trajs[[jmark]])
    names(traj.ctypes) <- traj.ctypes
    dat.trajs.long <- lapply(traj.ctypes, function(traj.ctype){
      dat.trajs.sub <- dat.trajs[[jmark]][[traj.ctype]] %>%
        filter(is.ctype) %>%
        dplyr::select(cell, ctype.ordered, ptime) %>%
        mutate(traj = traj.ctype, 
               colcode = ctype2col[[traj.ctype]])
    }) %>%
      bind_rows() %>%
      left_join(., dat.signal, by = "cell") %>%
      left_join(., subset(dat.meta.lst[[jmark]], select = c(cell, batch)), by = "cell")
    
    # order so Granulocyte is first? 
    dat.trajs.long <- dat.trajs.long %>%
      arrange(traj == jctype) %>%
      mutate(jalpha = ifelse(traj == jctype, 1, 0.25))
    
    m.traj <- ggplot(dat.trajs.long, aes(x = ptime, y = signal, color = colcode, alpha = jalpha)) + 
      geom_point() + 
      scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                            guide = "legend") + 
      theme_bw() + 
      ggtitle(paste(jmark, jctype)) + 
      theme(aspect.ratio=jratio, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m.traj)
    
    m.traj.batch <- ggplot(dat.trajs.long, aes(x = ptime, y = signal, color = colcode, alpha = jalpha)) + 
      geom_point() + 
      ggtitle(paste(jmark, jctype)) + 
      scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                            guide = "legend") + 
      theme_bw() + 
      facet_wrap(~batch) + 
      theme(aspect.ratio=jratio, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m.traj.batch)
    
    m.traj.batch2 <- ggplot(dat.trajs.long, aes(x = ptime, y = signal, color = batch, alpha = jalpha)) + 
      geom_point(alpha = 0.25) + 
      geom_point() + 
      theme_bw() + 
      ggtitle(paste(jmark, jctype)) + 
      facet_wrap(~traj) + 
      theme(aspect.ratio=jratio, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m.traj.batch2)
    
    dev.off()
  }
}



