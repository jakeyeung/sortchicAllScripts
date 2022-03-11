# Jake Yeung
# Date of Creation: 2022-01-10
# File: ~/projects/scchic/scripts/macbook_analysis_2021/new_experiments/qc_K562_create_metadata_for_spliitting_bams.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
jdate <- "2022-01-10"
indir <- "/Users/yeung/data/scchic/from_cluster_2021/new_experiments/filtered_count_tables_for_LDA/K562"
outdir <- "/Users/yeung/data/scchic/from_cluster_2021/new_experiments/filtered_count_tables_for_LDA/K562/metadata_for_splitting"

inf.meta.lst <- lapply(jmarks, function(jmark){
  file.path(indir, paste0("K562_", jmark), paste0("qc_metadata_new_only.K562_", jmark, ".0.8_0.5_3000.", jdate, ".txt"))
})

outf.meta.lst <- lapply(jmarks, function(jmark){
  outf <- file.path(outdir, paste0("metadata_for_splitting.K562_", jmark, ".0.8_0.5_3000.", Sys.Date() , ".txt"))
  assertthat::assert_that(!file.exists(outf))
  return(outf)
})

dat.meta.lst <- lapply(jmarks, function(jmark){
  jinf <- inf.meta.lst[[jmark]]
  fread(jinf)  %>%
    rowwise() %>%
    dplyr::mutate(cluster = ifelse(as.character(is.good), "GoodCells", "BadCells")) %>%
    dplyr::select(samp, cluster) %>%
    dplyr::rename(cell = samp) %>%
    mutate(mark = jmark)
})

# row 1: header
# column 1: cell 
# column 2: "cluster"
# column 3: mark

for (jmark in jmarks){
  fwrite(dat.meta.lst[[jmark]], file = outf.meta.lst[[jmark]], sep = "\t")
}



