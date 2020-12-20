# Jake Yeung
# Date of Creation: 2020-11-30
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/3-load_count_tables_fill_zeros_for_LDA.same_annot_file.R
# 


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

cbind.fill<-function(..., fill = 0){
  nm <- list(...) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(data = fill, n-nrow(x), ncol(x))))) 
}

cbind.fill.lst <- function(mats.lst, all.rnames, fill = 0){
  mats.lst.filled <- lapply(mats.lst, function(mat.tmp){
    missing.rnames <- all.rnames[!all.rnames %in% rownames(mat.tmp)]
    mat.tmp.to.fill <- matrix(data = fill, nrow = length(missing.rnames), ncol = ncol(mat.tmp), dimnames = list(missing.rnames, colnames(mat.tmp)))
    mat.tmp.bind <- rbind(mat.tmp, mat.tmp.to.fill)
    mat.tmp.bind <- mat.tmp.bind[all.rnames, ]
    return(mat.tmp.bind)
  })
  return(do.call(cbind, mats.lst.filled))
}

hubprefix <- "/home/jyeung/hub_oudenaarden"
outdir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BMround2.from_peaks.same_annot_file")
dir.create(outdir)

indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/count_tables_from_peaks.same_annot_file")
assertthat::assert_that(dir.exists(indir))

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

mats.bymark.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  infs.bymark <- list.files(path = indir, pattern = paste0("BM_round1_round2_merged_", jmark, "z*.txt"), full.names = TRUE)
  names(infs.bymark) <- sapply(infs.bymark, basename)
  
  mats.bymark <- lapply(infs.bymark, function(inf){
    ReadMatTSSFormat(inf, add.coord = TRUE, sort.rnames = TRUE)
  })
  
  jall.rnames <- sort(unique(unlist(lapply(mats.bymark, rownames))))
  mats.bymark.fill <- cbind.fill.lst(mats.bymark, all.rnames = jall.rnames, fill = 0)
  return(mats.bymark.fill)
})


# Write to output ---------------------------------------------------------

for (jmark in jmarks){
  print(jmark)
  outf <- file.path(outdir, paste0("count_mat_from_hiddendomains.", jmark, ".rds"))
  saveRDS(mats.bymark.lst[[jmark]], file = outf)
}




