# Jake Yeung
# Date of Creation: 2022-02-17
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/15-split_mat_by_trajs.R
# Split mat by trajs, then re-run LDA to get a better pseudotime 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks

# Load metas --------------------------------------------------------------

indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_cleaned_up"
infs.meta <- lapply(jmarks, function(jmark){
  file.path(indir.meta, paste0("metadata.", jmark, ".txt"))
})
dats.meta <- lapply(infs.meta, fread)


# Load trajs --------------------------------------------------------------

trajs.lst <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/trajs_dims_ctypes_list_output.2022-02-08.rds")

# # add "GMP" to "DCs" traj
# for (jmark in jmarks){
#   
# }

# Load mats ---------------------------------------------------------------

countmat.dir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/countmats_cleaned_up"
mats.tss <- lapply(jmarks, function(jmark){
  readRDS(file.path(countmat.dir, paste0("count_mat_cleaned_tss.", jmark, ".rds")))
})

mats.dynbins <- lapply(jmarks, function(jmark){
  readRDS(file.path(countmat.dir, paste0("count_mat_cleaned_dynbins.", jmark, ".rds")))
})
lapply(mats.dynbins, dim)

# Split mats by trajs  ----------------------------------------------------

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/countmats_cleaned_up/filter_by_trajs"
for (jmark in jmarks){
  print(jmark)
  print("Orig size")
  print(dim(mats.tss[[jmark]]))
  ctypes <- trajs.lst[[jmark]]$ctypes.names
  for (ctype in ctypes){
    print(ctype)
    ctypes.filt <- trajs.lst[[jmark]]$ctypes.lst[[ctype]]
    cells.keep <- subset(dats.meta[[jmark]], ctype.from.LL %in% ctypes.filt)$cell
    mat.tss.filt <- mats.tss[[jmark]][, cells.keep]
    print(dim(mat.tss.filt))
    # save output
    outf.tmp <- file.path(outdir, paste0("countmat_tss.", ctype, ".", jmark, ".", Sys.Date(), ".rds"))
    saveRDS(mat.tss.filt, file = outf.tmp)
    
  }
}

outdir.dynbins <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/countmats_cleaned_up/filter_by_trajs_dynbins"
# save allbins
for (jmark in jmarks){
  print(jmark)
  print("Orig size")
  print(dim(mats.dynbins[[jmark]]))
  ctypes <- trajs.lst[[jmark]]$ctypes.names
  for (ctype in ctypes){
    print(ctype)
    ctypes.filt <- trajs.lst[[jmark]]$ctypes.lst[[ctype]]
    cells.keep <- subset(dats.meta[[jmark]], ctype.from.LL %in% ctypes.filt)$cell
    mat.dynbins.filt <- mats.dynbins[[jmark]][, cells.keep]
    print(dim(mat.dynbins.filt))
    # save output
    outf.tmp <- file.path(outdir, paste0("countmat_dynbins.", ctype, ".", jmark, ".", Sys.Date(), ".rds"))
    saveRDS(mat.dynbins.filt, file = outf.tmp)
  }
}


