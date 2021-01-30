# Jake Yeung
# Date of Creation: 2020-12-28
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/correct_batch_effects/1-correct_batch_effects_from_peaks.from_project.R
# 

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
    inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/tm_results_from_projs/tm_result_old_to_new.", jmark, ".2020-12-28.RData"))
  } else {
    inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/tm_results_from_projs/tm_result_new_to_old.", jmark, ".2020-12-28.RData"))
  }
  assertthat::assert_that(file.exists(inf.lda.tmp))
  
  return(inf.lda.tmp)
})

out.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.lda <- inf.lda.lst[[jmark]]
  load(inf.lda, v=T)  # tm.result and count.mat
  return(list(tm.result = tm.result, count.mat = count.mat))
})

count.mat.lst <- lapply(out.lst, function(out){
  out$count.mat
})

# Load meta data  ---------------------------------------------------------

indir.metas <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BMround2.from_peaks.sitecount_mat.split_old_and_new/add_experi")
dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.metadata.2020-12-28.with_experi.txt")
  fread(file.path(indir.metas, fname))
})

# add jrep2 for batch correction?
dat.metas <- lapply(jmarks, function(jmark){
  dat.metas.tmp <- dat.metas[[jmark]]
  if (jmark != "H3K9me3"){
    dat.metas.tmp$jrep2 <- sapply(dat.metas.tmp$jrep, function(x) ifelse(x == "rep1old", "zold", "anew"))
  } else {
    dat.metas.tmp$jrep2 <- sapply(dat.metas.tmp$jrep, function(x) ifelse(x != "rep1old", "zold", "anew"))
  }
  return(dat.metas.tmp)
})


# Select bins and correct -------------------------------------------------

imputed.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  log2(t(out.lst[[jmark]]$tm.result$topics %*% out.lst[[jmark]]$tm.result$terms))
})


imputed.long.lst <- lapply(jmarks, function(jmark){
  rnames.keep <- rownames(imputed.lst[[jmark]])  # keep all
  jmat.filt <- imputed.lst[[jmark]][rnames.keep, ] %>%
    data.table::melt()
  colnames(jmat.filt) <- c("rname", "cell", "log2exprs")
  jmat.filt <- jmat.filt %>%
    left_join(., dat.metas[[jmark]])
  return(jmat.filt)
})



# Save output -------------------------------------------------------------

print("Saving output...")
system.time(
  for (jmark in jmarks){
    print(jmark)
    outrds <- file.path(outdir, paste0("imputed_long.", jmark, ".from_project.rds"))
    saveRDS(imputed.long.lst[[jmark]], file = outrds)
  }
)





