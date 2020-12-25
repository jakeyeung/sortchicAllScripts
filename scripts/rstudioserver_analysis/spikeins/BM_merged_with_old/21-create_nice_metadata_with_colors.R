# Jake Yeung
# Date of Creation: 2020-12-13
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/21-create_nice_metadata_with_colors.R
# Create metadata tables for all 4 marks after cleaning up 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"


inf.k4me1 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.H3K4me1.txt")
inf.k4me3 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.H3K4me3.txt")
inf.k27me3 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.H3K27me3_techrepmerged/BM_rep2_rep3reseq_H3K27me3.2020-12-10.txt")
inf.k9me3 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/BM_celltypes.bincutoff_0.binskeep_1000.byplate.szname_none.niter_500.reorder_rownames.dupfilt.2020-11-23.H3K9me3.txt")
assertthat::assert_that(file.exists(inf.k4me1))
assertthat::assert_that(file.exists(inf.k4me3))
assertthat::assert_that(file.exists(inf.k27me3))
assertthat::assert_that(file.exists(inf.k9me3))

inf.lst <- list(H3K4me1 = inf.k4me1,
                H3K4me3 = inf.k4me3,
                H3K27me3 = inf.k27me3,
                H3K9me3 = inf.k9me3)

dat.meta.lst <- lapply(jmarks, function(jmark){
  inf <- inf.lst[[jmark]]
  dat.annot <- fread(inf)
  if (jmark == "H3K27me3"){
    dat.annot <- dat.annot %>%
      dplyr::rename(cluster.more = cluster,
                    cluster = cluster.fewer)
    
  }
  if (jmark == "H3K9me3"){
    dat.annot <- dat.annot %>%
      rowwise() %>%
      mutate(cluster = gsub("Eryth", "Eryths", cluster),
             cluster = gsub("Lymphoid", "Bcells", cluster))
  }
  
  dat.annot <- dat.annot %>%
    rowwise() %>%
    mutate(plate = ClipLast(x = cell,jsep = "_"))
  
  dat.annot <- dat.annot %>%
    arrange(cluster, jrep)
})

clstrs.uniq <- unique(dat.meta.lst$H3K4me1$cluster)

