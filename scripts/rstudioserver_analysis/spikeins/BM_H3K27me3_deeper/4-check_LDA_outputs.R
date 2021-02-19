# Jake Yeung
# Date of Creation: 2021-02-14
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_deeper/4-check_LDA_outputs.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/variance_clusters"

outpdf <- file.path(outdir, paste0("check_variance_across_cells_and_clusters.", Sys.Date(), ".pdf"))

pdf(file = outpdf)

# Metadata ----------------------------------------------------------------

indir.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage"


dat.metas <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.meta <- file.path(indir.meta, paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt"))
  fread(inf.meta)
})


# Load data  --------------------------------------------------------------

inf.lda.lst <- lapply(jmarks, function(jmark){
  if (jmark != "H3K27me3"){
    inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj"))
  } else {
    inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_", jmark, ".cleaned.varfilt_2.K-30.Robj"))
  }
  assertthat::assert_that(file.exists(inf.lda.tmp))
  return(inf.lda.tmp)
})

out.lda.lst <- lapply(inf.lda.lst, function(inf){
  load(inf, v=T)
  return(out.lda)
})


# Get imputed ------------------------------------------------------------

tm.result.lst <- lapply(out.lda.lst, function(jout){
  posterior(jout)
})


# Get var -----------------------------------------------------------------

dat.imputed.lst <- lapply(tm.result.lst, function(tm.result){
  t(log2(tm.result$topics %*% tm.result$terms))
})



# Get var -----------------------------------------------------------------

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

dat.var.lst <- lapply(dat.imputed.lst, function(dat.imputed){
  CalculateVarAll(dat.imputed, jchromos)
})



# Add to umap and plot  ---------------------------------------------------


jmark <- "H3K27me3"

dat.merge.lst <- lapply(jmarks, function(jmark){
  left_join(dat.metas[[jmark]], dat.var.lst[[jmark]])
})



m.lst <- lapply(jmarks, function(jmark){
  ggplot(dat.merge.lst[[jmark]], aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + 
    scale_color_viridis_c() + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

print(m.lst)


m.bar <- lapply(jmarks, function(jmark){
  ggplot(dat.merge.lst[[jmark]], aes(x = forcats::fct_reorder(.f = cluster, .x = cell.var.within.sum.norm), y = cell.var.within.sum.norm)) + 
    geom_point() + 
    geom_boxplot() + 
    ylab("Intrachromosomal Variance") + 
    xlab("") + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
print(m.bar)


dev.off()
