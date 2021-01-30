# Jake Yeung
# Date of Creation: 2020-12-28
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/correct_batch_effects/1-correct_batch_effects_from_peaks_from_LDA.R
# description
rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(topicmodels)

library(JFuncs)

library(hash)
library(igraph)
library(umap)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(DescTools)

library(topicmodels)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

ncores <- 4
hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_batch_correction_output"
dir.create(outdir)

# Load LDA outputs --------------------------------------------------------

# for K27me3 
inf.lda.lst <- lapply(jmarks, function(jmark){
  if (jmark != "H3K9me3"){
    # inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/tm_results_from_projs/tm_result_old_to_new.", jmark, ".2020-12-28.RData"))
    inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins.from_sitecount_mat.from_same_annot_file/lda_outputs.count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.K-30.binarize.FALSE/ldaOut.count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.K-30.Robj"))
  } else {
    # inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/tm_results_from_projs/tm_result_new_to_old.", jmark, ".2020-12-28.RData"))
    inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins/lda_outputs.count_mat_from_hiddendomains.", jmark, ".filtNAcells_allbins.K-30.binarize.FALSE/ldaOut.count_mat_from_hiddendomains.", jmark, ".filtNAcells_allbins.K-30.Robj"))
  }
  assertthat::assert_that(file.exists(inf.lda.tmp))
  
  return(inf.lda.tmp)
})

out.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.lda <- inf.lda.lst[[jmark]]
  load(inf.lda, v=T)  # out.lda, count.mat
  tm.result <- posterior(out.lda)
  return(list(tm.result = tm.result, count.mat = count.mat))
})

count.mat.lst <- lapply(out.lst, function(out){
  out$count.mat
})



# Save output -------------------------------------------------------------

print("Saving output")
print(jmarks)

system.time(
  for (jmark in jmarks){
    print(jmark)
    outrds <- file.path(outdir, paste0("count_mat.", jmark, ".from_LDA.rds"))
    saveRDS(count.mat.lst[[jmark]], file = outrds)
  }
)






