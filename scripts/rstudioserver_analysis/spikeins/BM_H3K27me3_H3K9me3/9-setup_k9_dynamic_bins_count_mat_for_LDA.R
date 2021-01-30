# Jake Yeung
# Date of Creation: 2021-01-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_H3K9me3/9-setup_k9_dynamic_bins_count_mat_for_LDA.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"

# Load count mats ---------------------------------------------------------


# load LDA 
outs.all.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  if (jmark != "H3K27me3"){
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj"))
  } else {
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.Robj"))
  }
  load(inf.lda, v=T)
  return(list(out.lda = out.lda, count.mat = count.mat))
})

out.lda.lst <- lapply(outs.all.lst, function(jout) jout$out.lda)


count.mat.lst <- lapply(outs.all.lst, function(jout){
  jmat <- jout$count.mat
  # jmat <- sweep(jmat, MARGIN = 2, STATS = colSums(jmat), FUN = "/")
  # jmat <- BinarizeMatrix(jmat)
})



# Load bins ---------------------------------------------------------------

inf.bed <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/tables_H3K9me3_dynamic_bins/coords_H3K9me3_dynamic_bins.bed"
dat.bed <- fread(inf.bed)

bins.keep <- dat.bed$V4



# Write filltered bins  ---------------------------------------------------


count.mat.filt.lst <- lapply(count.mat.lst, function(jmat.tmp){
  print(dim(jmat.tmp))
  rnames.keep <- rownames(jmat.tmp) %in% bins.keep
  jmat.tmp.filt <- jmat.tmp[rnames.keep, ]
  print(dim(jmat.tmp.filt))
  return(jmat.tmp.filt)
})

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k9_dynamic_50kb_bins"
for (jmark in jmarks){
  print(jmark)
  outrds <- file.path(outdir, paste0("count_mat_k9_dynamic_bins_50kb.", jmark, ".", Sys.Date(), ".rds"))
  saveRDS(count.mat.filt.lst[[jmark]], file = outrds)
}



