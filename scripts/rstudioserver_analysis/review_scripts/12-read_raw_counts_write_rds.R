# Jake Yeung
# Date of Creation: 2021-08-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/12-read_raw_counts_write_rds.R
# 

# Load scATAC-seq 10kb, 50kb, peaks. Load sortChIC 10kb, 50kb, peaks?


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# Load  -------------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"


mapq <- 30
bnames <- c("BoneMarrow_62016", "BoneMarrow_62216")
names(bnames) <- bnames

jsuffs <- c("peaks_filt", "10kb_TSS", "50kb_TSS")
names(jsuffs) <- jsuffs

outdir <- file.path(hubprefix, paste0("jyeung/data/scChiC/public_data/Cusanovich_2018/rds_objs"))
dat.mat.lst <- lapply(jsuffs, function(jsuff){
  infs.lst <- lapply(bnames, function(bname){
    # inf.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/public_data/Cusanovich_2018/count_tables/", bname, ".peaks_filt.mapq_", mapq, ".count_table.txt"))
    # inf.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/public_data/Cusanovich_2018/count_tables/", bname, ".10kb_TSS.mapq_", mapq, ".count_table.txt"))
    # inf.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/public_data/Cusanovich_2018/count_tables/", bname, ".50kb_TSS.mapq_", mapq, ".count_table.txt"))
    inf.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/public_data/Cusanovich_2018/count_tables/", bname, ".", jsuff, ".mapq_", mapq, ".count_table.txt"))
    assertthat::assert_that(file.exists(inf.tmp))
    return(inf.tmp)
  })
  
  outf.tmp <- file.path(outdir, paste0("countmat_", jsuff, ".scATACseq.rds"))
  
  dat.mat <- lapply(infs.lst, function(jinf){
    mat.tmp <- fread(jinf)
    mat.sparse <- Matrix(as.matrix(mat.tmp[, -1]), sparse = TRUE)
    rownames(mat.sparse) <- mat.tmp$coordname
    return(mat.sparse)
  }) 
  dat.mat <- do.call(cbind, dat.mat)
  print(dim(dat.mat))
  saveRDS(dat.mat, file = outf.tmp)
  return(dat.mat)
})


# Load bone marrow  ------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K9me3"); names(jmarks) <- jmarks

inmain <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/count_tables_from_Cusanovich_atac_peaks")

mats.lst <- lapply(jmarks, function(jmark){
  infs.tmp <- list.files(inmain, pattern = paste0(jmark, ".*.txt"), full.names = TRUE)
  mat.lst.tmp <- lapply(infs.tmp, function(inf.tmp){
    scchicFuncs::ReadMatTSSFormat(inf.tmp)
  })
  all.rnames <- unique(unlist(lapply(mat.lst.tmp, function(jmat) rownames(jmat))))
  JFuncs::cbind.fill.lst(mat.lst.tmp, all.rnames = all.rnames)
})
# print(infs.lst)

outdir <- file.path(hubprefix, paste0("jyeung/data/scChiC/public_data/Cusanovich_2018/rds_objs"))
for (jmark in jmarks){
  outrds <- file.path(outdir, paste0("countmat_CusanovichPeaks.sortChIC.", jmark, ".rds"))
  mat.tmp <- mats.lst[[jmark]]
  print(dim(mat.tmp))
  saveRDS(mat.tmp, file = outrds)
}
