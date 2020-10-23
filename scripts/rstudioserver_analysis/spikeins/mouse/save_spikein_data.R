# Jake Yeung
# Date of Creation: 2020-09-13
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/mouse/save_spikein_data.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



jspikeinchromo <- "J02459.1"
hubprefix <- "/home/jyeung/hub_oudenaarden"
indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5046_BM/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters")
assertthat::assert_that(dir.exists(indir))


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
jchromos <- paste("", seq(19), sep = "")

dat.out.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  infs.rz <- list.files(indir, pattern = paste0("PZ-ChIC-mouse-BM-", jmark, ".*.RZ.csv"), full.names = TRUE)
  infs.chromo <- list.files(indir, pattern = paste0("PZ-ChIC-mouse-BM-", jmark, ".*.NoChromo.csv"), full.names = TRUE)
  infs.bins <- list.files(indir, pattern = paste0("PZ-ChIC-mouse-BM-", jmark, ".*.binsize_50000.csv"), full.names = TRUE)
  dat.rz <- lapply(infs.rz, function(inf.rz){
    dat.tmp <- ReadLH.SummarizeTA(inf.rz)
  }) %>%
    bind_rows()
  
  dat.chromo <- lapply(infs.chromo, function(inf.bychromo){
    dat.tmp <- GetChromoCounts(inf.bychromo, spikeinchromo = jspikeinchromo, chromos.keep = jchromos) %>%
      filter(chromo == "1")
    return(dat.tmp)
  }) %>%
    bind_rows()
  
  dat.rzchromo <- left_join(dat.rz, dat.chromo, by = "samp")
  dat.rzchromo$mark <- jmark
  
  
  dat.mat <- lapply(infs.bins, function(inf.bin){
    dat.tmp <- ReadMatSlideWinFormat(inf.bin)
  })
  rows.common <- Reduce(intersect, x = lapply(dat.mat, rownames))
  
  dat.mat.filt.lst <- lapply(dat.mat, function(jdat){
    jdat[rows.common, ]
  }) 
  dat.mat.filt <- do.call(cbind, dat.mat.filt.lst)
  
  # filter good cells
  return(list(dat.rz = dat.rzchromo, dat.mat = dat.mat.filt))
})


dat.spikein.mat <- lapply(jmarks, function(jmark){
  dat.out.lst[[jmark]]$dat.rz
}) %>%
  bind_rows() %>%
  as.data.frame()
rownames(dat.spikein.mat) <- dat.spikein.mat$samp

fwrite(dat.spikein.mat, file = "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda.chromo2spikeinfilt.VAN5046/spikein_info.txt")


