# Jake Yeung
# Date of Creation: 2022-04-06
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/7-annotate_celltypes_from_reference_dynamic_bins_update_ctypes_from_LDA_all_marks.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(topicmodels)
library(AnnotateCelltypes)
library(parallel)
library(hash)

EuclideanDistance <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))


jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
jmark <- jmarks[[1]]

jmarks.active <- c("k4me1", "k4me3"); names(jmarks.active)
jmarks.repress <- c("k27me3", "k9me3"); names(jmarks.repress)

jsuffix <- "dynamicbins"
outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_newannot"
dir.create(outdir)
# for (jmark in jmarks){

ctypes.to.annotate <- c("Basophils", "Bcells", "CMP/GMP", "DCs", "Eryths", "Granulocytes", "HSCs", "MEP", "Monocytes", "MPPs", "NKs", "pDCs")
names(ctypes.to.annotate) <- ctypes.to.annotate

LL.all.lst <- mclapply(jmarks, function(jmark){
  print(jmark)
  
  outrdata <- file.path(outdir, paste0("LLmat_by_batch_", jsuffix, "_", jmark, ".", Sys.Date(), ".RData"))
  
  # Load metadata -----------------------------------------------------------
  
  if (jmark %in% jmarks.active){
    indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/varfilt"
    inf.meta <- file.path(indir.meta, paste0("metadata_celltyping_", jmark, ".allbins.2022-02-01.txt"))
    dat.meta <- fread(inf.meta)
  } else if (jmark %in% jmarks.repress){
    indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/repressive_cleaned_check_eryths"
    inf.meta <- file.path(indir.meta, paste0("metadata_celltyping_", jmark, ".dynamicbins.2022-04-06.txt"))
    dat.meta <- fread(inf.meta)
  }
  
  inf.meta.newannot <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata_new_prog_annots/metadata_newProgenitors.", jmark, ".2022-03-01.txt")
  dat.meta.newannot <- fread(inf.meta.newannot)
  print(unique(sort(dat.meta.newannot$ctype)))
  
  ggplot(dat.meta, aes(x = umap1, y = umap2, color = batch)) + 
    geom_point() + 
    facet_wrap(~ctype) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # reclassify celltypes for AllCells, HSPCs, IL7RLinNeg, LinNeg, LSK
  
  # Load rawcounts ----------------------------------------------------------
  
  if (jmark %in% jmarks.active){
    indir.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_", jsuffix, ".out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28")
    inf.ldaout <- file.path(indir.ldaout, paste0("ldaOut.count_mat_var_filt_", jsuffix, ".out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28.Robj"))
    load(inf.ldaout, v=T)
  } else if (jmark %in% jmarks.repress){
    indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_repressive_cleaned_check_eryths_keep_old/ldaAnalysis_fripfilt_varfilt_binfilt"
    dname <- paste0("count_mat_cleaned_dynbins.", jmark, ".2022-04-05")
    inf <- file.path(indir, paste0("lda_outputs.", dname), paste0("ldaOut.", dname, ".Robj"))
    assertthat::assert_that(file.exists(inf))
    load(inf, v=T)
  }
  tm <- posterior(out.lda)
  dat.impute.linear <- t(tm$topics %*% tm$terms)
  
  # Annotate celltypes ------------------------------------------------------
  
  # use Euclidean distance on the topics matrix
  
  # assume ctype from new data is "fixed", if NaN then reannotate to new ctype
  dat.meta.newannot.sub.nan <- subset(dat.meta.newannot, ctype == "NaN") %>%
    rowwise() %>%
    mutate(ctype = ifelse(ctype == "CMP/GMP", "CMP", ctype))
  
 cells.new.toreannotate <- dat.meta.newannot.sub.nan$cell
  
  # dangling cells
  dat.meta.new.sub <- subset(dat.meta, batch == "New")
  cells.toadd <- subset(dat.meta.new.sub, !cell %in% dat.meta.newannot$cell)$cell
  print(paste("Cells to add:", length(cells.toadd)))
  
  cells.new.toreannotate <- c(cells.new.toreannotate, cells.toadd)
  names(cells.new.toreannotate) <- cells.new.toreannotate
  
  print(length(cells.new.toreannotate))
  
  dat.meta.newannot.sub.notnan <- subset(dat.meta.newannot, ctype != "NaN")
  table(dat.meta.newannot.sub.notnan$ctype)
  new.groundtruth.cells <- dat.meta.newannot.sub.notnan$cell
  names(new.groundtruth.cells) <- new.groundtruth.cells
  
  dat.meta.sub.old.groundtruth <- subset(dat.meta, batch == "Old" & ctype %in% c("Basophils", "Bcells", "Eryths", "Granulocytes", "NKs", "pDCs"))
  old.groundtruth.cells <- dat.meta.sub.old.groundtruth$cell
  
  cell2label <- hash::hash(dat.meta.newannot.sub.notnan$cell, dat.meta.newannot.sub.notnan$ctype)
  cell2label.old <- hash::hash(dat.meta.sub.old.groundtruth$cell, dat.meta.sub.old.groundtruth$ctype)
  cell2label.merge <- hash::hash(c(keys(cell2label), keys(cell2label.old)), c(values(cell2label), values(cell2label.old)))
  
  # for old data: reannotate all cells
  cells.old.toreannotate <- subset(dat.meta, batch == "Old")$cell
  names(cells.old.toreannotate) <- cells.old.toreannotate
  
  topics.mat.new.groundtruth <- tm$topics[new.groundtruth.cells, ]
  topics.mat.old.groundtruth <- tm$topics[c(old.groundtruth.cells, new.groundtruth.cells), ]
  
  N <- 15
  
  # reannotate new cells
  closest.cells.lst.new <- lapply(cells.new.toreannotate, function(jcell){
    topics.vec <- tm$topics[jcell, ]
    euc.dists <- apply(topics.mat.new.groundtruth, 1, function(ref.row){
      EuclideanDistance(topics.vec, ref.row)
    }) %>%
      sort(decreasing = FALSE)
    # get labels for top N
    closest.cells <- names(euc.dists[1:N])
    closest.cells.labels <- sapply(closest.cells, function(x) AssignHash(x = x, jhash = cell2label, null.fill = NA))
    assertthat::assert_that(all(!is.na(closest.cells.labels)))
    return(closest.cells.labels)
  })
  
  closest.cells.best.new <- lapply(closest.cells.lst.new, function(closests){
    tabsort <- sort(table(closests), decreasing = TRUE)
    return(names(tabsort)[[1]])
  }) %>%
    unlist()
  
  # reannotate old cells
  closest.cells.lst.old <- lapply(cells.old.toreannotate, function(jcell){
    # print(jcell)
    topics.vec <- tm$topics[jcell, ]
    euc.dists <- apply(topics.mat.old.groundtruth, 1, function(ref.row){
      EuclideanDistance(topics.vec, ref.row)
    }) %>%
      sort(decreasing = FALSE)
    # get labels for top N
    closest.cells <- names(euc.dists[1:N])
    closest.cells.labels <- sapply(closest.cells, function(x) AssignHash(x = x, jhash = cell2label.merge, null.fill = NA))
    assertthat::assert_that(all(!is.na(closest.cells.labels)))
    return(closest.cells.labels)
  })
  
  closest.cells.best.old <- lapply(closest.cells.lst.old, function(closests){
    tabsort <- sort(table(closests), decreasing = TRUE)
    return(names(tabsort)[[1]])
  }) %>%
    unlist()
  
  # combine
  # closest.cells.out <- list("New" = closest.cells.lst, "Old" = closest.cells.lst.old)
  # closest.cells.best.merge <- 
  save(closest.cells.lst.new, closest.cells.lst.old, closest.cells.best.new, closest.cells.best.old, file = outrdata)
  print(paste("Done for", jmark))
  return(list(closest.cells.best.new, closest.cells.best.old))
  # }
}, mc.cores = length(jmarks))
  
print(LL.all.lst)  
  


