# Jake Yeung
# Date of Creation: 2021-08-28
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/13-LDA_atacseq_downstream.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)
library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

hubprefix <- "/home/jyeung/hub_oudenaarden"


# Function ----------------------------------------------------------------

FixNames <- function(jpeak){
  ending <- strsplit(jpeak, "-")[[1]][[2]]
  endinghalf <- substr(ending, 1, nchar(ending) / 2)
  jstart <- strsplit(jpeak, "-")[[1]][[1]]
  jpeak.wrangled <- paste(jstart, endinghalf, sep = "-")
  jpeak.atac <- paste(jpeak.wrangled, jpeak.wrangled, sep = ";")
  return(jpeak.atac)
}

# Load dat ----------------------------------------------------------------


inf.atac <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisRevisions_countmat_10kb_TSS.scATACseq.K-30/lda_outputs.countmat_10kb_TSS.scATACseq.K-30.binarize.FALSE/ldaOut.countmat_10kb_TSS.scATACseq.K-30.Robj"))
load(inf.atac, v=T)

tm.result <- posterior(out.lda)
dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)


# Add meta  ---------------------------------------------------------------

inf.meta <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Cusanovich_2018/metadata/cell_metadata.bonemarrow_only.withheader.txt")
dat.meta <- fread(inf.meta)
dat.meta.filt <- subset(dat.meta, select = c(cell, cell_label, tissue.replicate))
dat.meta.filt$barcode <- paste(dat.meta.filt$tissue.replicate, dat.meta.filt$cell, sep = ".")



# Annotate ----------------------------------------------------------------

dat.umap.annot <- left_join(dat.umap, dat.meta.filt)

m.atac <- ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = cell_label)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("scATACseq") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")



# Check sortChIC ----------------------------------------------------------

jmarks <- c("H3K4me1"); names(jmarks) <- jmarks

infs.chic.lst <- lapply(jmarks, function(jmark){
  # inf.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisRevisions_countmat_CusanovichPeaks.sortChIC.", jmark, ".K-30/lda_outputs.countmat_CusanovichPeaks.sortChIC.", jmark, ".K-30.binarize.FALSE/ldaOut.countmat_CusanovichPeaks.sortChIC.", jmark, ".K-30.Robj"))
  inf.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.from_TSS.dist_10000.RemoveDupRows/lda_outputs.count_mat_from_TSS.H3K4me1.dist_10000.K-30.binarize.FALSE/ldaOut.count_mat_from_TSS.H3K4me1.dist_10000.K-30.Robj"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

out.lda.lst <- lapply(infs.chic.lst, function(jinf){
  load(jinf, v=T)
  return(out.lda)
})

count.mat.sortchic.lst <- lapply(infs.chic.lst, function(jinf){
  load(jinf, v=T)
  # rownames(count.mat) <- sapply(rownames(count.mat), FixNames)
  return(count.mat)
})


tm.result.lst <- lapply(out.lda.lst, function(jout){
  posterior(jout)
})


dat.umap.lst <- lapply(tm.result.lst, function(jtm){
  DoUmapAndLouvain(jtm$topics, jsettings)
})


cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

m.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.umap.lst[[jmark]]
  ggplot(jdat, aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    scale_color_manual(values = cbPalette) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})



# Get meta ----------------------------------------------------------------

infs.meta <- lapply(jmarks, function(jmark){
  inf.meta <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt"))
  assertthat::assert_that(file.exists(inf.meta))
  return(inf.meta)
})

dat.meta.chic.lst <- lapply(infs.meta, function(jinf){
  dat.meta.tmp <- fread(jinf)
  
  dat.meta.tmp <- dat.meta.tmp %>%
    dplyr::select(c(cell, cluster, stypecol, clustercol, batch))
})

dat.umap.chic.annot.lst <- lapply(jmarks, function(jmark){
  left_join(dat.umap.lst[[jmark]], dat.meta.chic.lst[[jmark]])
})

m.celltypes.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.umap.chic.annot.lst[[jmark]]
  ggplot(jdat, aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + 
    scale_color_manual(values = cbPalette) + 
    theme_bw() + 
    ggtitle(paste("sortChIC", jmark)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})


# Chekc commons -----------------------------------------------------------

commons <- intersect(rownames(count.mat), rownames(count.mat.sortchic.lst[[1]]))

assertthat::assert_that(length(commons) > 0)

# Write outputs -----------------------------------------------------------

count.mat.atac <- count.mat
count.mat.chic <- count.mat.sortchic.lst$H3K4me1
tm.atac <- tm.result
tm.chic <- tm.result.lst$H3K4me1
meta.atac <- dat.meta.filt
meta.chic <- dat.meta.chic.lst$H3K4me1

outdir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/post_submission/analysis_scrnaseq_atacseq")
outmatatac <- file.path(outdir, paste0("scATACseq_Cusanovich_10kbTSS_countmat.rds"))
outmatchic <- file.path(outdir, paste0("H3K4me1-sortChIC_Zeller_10kbTSS_countmat.rds"))

outldaatac <- file.path(outdir, paste0("scATACseq_Cusanovich_10kbTSS_LDA.rds"))
outldachic <- file.path(outdir, paste0("H3K4me1-sortChIC_Zeller_10kbTSS_LDA.rds"))

metaatac <- file.path(outdir, paste0("scATACseq_Cusanovich_meta.txt"))
metachic <- file.path(outdir, paste0("H3K4me1-sortChIC_Zeller_meta.txt"))


# write count mats
saveRDS(count.mat.atac, file = outmatatac)
saveRDS(count.mat.chic, file = outmatchic)
saveRDS(tm.atac, file = outldaatac)
saveRDS(tm.chic, file = outldachic)
fwrite(meta.atac, file = metaatac)
fwrite(meta.chic, file = metachic)

