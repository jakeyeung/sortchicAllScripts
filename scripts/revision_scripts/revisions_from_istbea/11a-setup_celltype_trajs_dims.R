# Jake Yeung
# Date of Creation: 2022-02-08
# File: ~/projects/scchic/scripts/scripts_analysis/11a-setup_celltype_trajs_dims.R
# description

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)

library(topicmodels)
library(DescTools)


# GetTrajs ----------------------------------------------------------------

GetTrajs <- function(jmark = "k4me1"){
  
  if (jmark != "k9me3"){
    hpcs <- c("HSCs", "LT", "ST", "MPPs")
    # granu tracj
    ctypes.granus <- c(hpcs, "CMP", "GMP", "Granulocytes")
    ctypes.bcells <- c(hpcs, "Bcells")
    
    ctypes.eryths <- c(hpcs, "MEP", "Eryths")
    
    ctypes.pdcs <- c(hpcs, "pDCs")
    
    ctypes.dcs <- c(hpcs, "DCs")
    
    ctypes.basos <- c(hpcs, "CMP", "GMP", "Basophils")
    
    ctypes.nks <- c(hpcs, "NKs")
    
    ctypes.monos <- c(hpcs, "CMP", "GMP", "Monocytes")
    
    ctypes.lst <- list(ctypes.granus, ctypes.bcells, ctypes.eryths, ctypes.pdcs, ctypes.dcs, ctypes.basos, ctypes.nks, ctypes.monos)
    names(ctypes.lst) <- sapply(ctypes.lst, function(x) x[length(x)])
    
    ctypes.names <- names(ctypes.lst); names(ctypes.names) <- ctypes.names
    
    if (jmark == "k4me1"){
      dims.granus <- c("dim1", "dim2")
      dims.bcells <- c("dim7", "dim8")
      dims.eryths <- c("dim2", "dim3")
      dims.pdcs <- c("dim8", "dim9")
      dims.dcs <- c("dim1", "dim2")
      dims.basos <- c("dim1", "dim2")
      dims.nks <- c("dim4", "dim5")
      dims.monos <- c("dim1", "dim2")
    } else if (jmark == "k4me3"){
      dims.granus <- c("dim3", "dim4")
      dims.bcells <- c("dim2", "dim3")
      dims.eryths <- c("dim3", "dim4")
      dims.pdcs <- c("dim14", "dim15")
      dims.dcs <- c("dim4", "dim5")
      dims.basos <- c("dim1", "dim2")
      dims.nks <- c("dim6", "dim7")
      dims.monos <- c("dim4", "dim5")
    } else if (jmark == "k27me3"){
      dims.granus <- c("dim3", "dim4")
      dims.bcells <- c("dim3", "dim4")
      dims.eryths <- c("dim3", "dim4")
      dims.pdcs <- c("dim3", "dim4")
      dims.dcs <- c("dim3", "dim4")
      dims.basos <- c("dim3", "dim4")
      dims.nks <- c("dim3", "dim4")
      dims.monos <- c("dim3", "dim4")
    }
    dims.lst  <- list(dims.granus, dims.bcells, dims.eryths, dims.pdcs, dims.dcs, dims.basos, dims.nks, dims.monos)
    names(dims.lst) <- ctypes.names
  } else if (jmark == "k9me3"){
    hpcs <- c("HSCs", "LT", "ST", "MPPs")
    # granu tracj
    ctypes.granus <- c(hpcs, "CMP", "GMP", "Granulocytes")
    dims.granus <- c("dim1", "dim2")
    
    ctypes.bcells <- c(hpcs, "Bcells")
    dims.bcells <- c("dim1", "dim2")
    
    ctypes.eryths <- c(hpcs, "MEP", "Eryths")
    dims.eryths <- c("dim5", "dim6")
    
    ctypes.pdcs <- c(hpcs, "pDCs")
    dims.pdcs <- c("dim1", "dim2")
    
    ctypes.monos <- c(hpcs, "CMP", "GMP", "Monocytes")
    dims.monos <- c("dim1", "dim2")
    
    ctypes.lst <- list(ctypes.granus, ctypes.bcells, ctypes.eryths, ctypes.pdcs, ctypes.monos)
    names(ctypes.lst) <- sapply(ctypes.lst, function(x) x[length(x)])
    
    ctypes.names <- names(ctypes.lst); names(ctypes.names) <- ctypes.names
    
    dims.lst  <- list(dims.granus, dims.bcells, dims.eryths, dims.pdcs, dims.monos)
    names(dims.lst) <- ctypes.names
  } 
  return(list(dims.lst = dims.lst, ctypes.lst = ctypes.lst, ctypes.names = ctypes.names))
}

# Define marks ------------------------------------------------------------


# jmark <- "k4me1"
jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks

# Define trajs ------------------------------------------------------------

trajs.out <- lapply(jmarks, function(jmark){
  GetTrajs(jmark)
})


# Save output ------------------------------------------------------------

outf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/trajs_dims_ctypes_list_output.", Sys.Date(), ".rds")
saveRDS(trajs.out, file = outf)
