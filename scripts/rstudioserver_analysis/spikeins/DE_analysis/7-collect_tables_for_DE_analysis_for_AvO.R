# Jake Yeung
# Date of Creation: 2021-02-02
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/DE_analysis/7-collect_tables_for_DE_analysis_for_AvO.R
# 


rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)

hubprefix <- "/home/jyeung/hub_oudenaarden"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


# Load metas  -------------------------------------------------------------


repremove <- "rep1old"
indir.metas <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage"

dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  inf <- file.path(indir.metas, fname)
  fread(inf)
})

dat.metas.filt <- lapply(jmarks, function(jmark){
  dat.metas.filt.tmp <- subset(dat.metas[[jmark]], jrep != repremove)
})

norm.factors.lst <- lapply(jmarks, function(jmark){
  dat.metas.filt.tmp <- dat.metas.filt[[jmark]]
  dat.metas.filt.tmp.sum <- dat.metas.filt.tmp %>%
    group_by(cluster) %>%
    summarise(spikein_cuts = sum(spikein_cuts)) %>%
    arrange(cluster)
  return(dat.metas.filt.tmp.sum)
})

# Load for TSS ------------------------------------------------------------

jtype <- "TSS"
jdist <- 10000
indir <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.from_", jtype, ".dist_", jdist))
assertthat::assert_that(dir.exists(indir))

out.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0("lda_outputs.count_mat_from_", jtype, ".", jmark, ".dist_", jdist, ".K-30.binarize.FALSE/ldaOut.count_mat_from_", jtype, ".", jmark, ".dist_", jdist, ".K-30.Robj")
  inf.tmp <- file.path(indir, fname)
  assertthat::assert_that(file.exists(inf.tmp))
  load(inf.tmp, v=T)
  tm.result <- posterior(out.lda)
  return(list(tm.result = tm.result, count.mat = count.mat))
})


# Get pseudobulks  --------------------------------------------------------

pbulk.tss.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  mat.tmp <- out.lst[[jmark]]$count.mat
  dat.metas.filt.tmp <- dat.metas.filt[[jmark]]
  norm.factors <- norm.factors.lst[[jmark]]
  cells.lst <- split(x = dat.metas.filt.tmp$cell, f = dat.metas.filt.tmp$cluster)
  pbulk.tmp.lst <- SumAcrossClusters(mat.tmp, cells.lst)
  pbulk.mat <- do.call(cbind, pbulk.tmp.lst)
  
  # make colnames alphabetiacl
  pbulk.mat <- pbulk.mat[, sort(colnames(pbulk.mat))]
  
  # get norm factors from spikeins 
  
  pbulk.mat.norm <- sweep(pbulk.mat, MARGIN = 2, STATS = norm.factors$spikein_cuts, FUN = "/")
  return(list(pbulk.mat = pbulk.mat, pbulk.mat.norm = pbulk.mat.norm, norm.factors = norm.factors))
})


# Load for bins  ----------------------------------------------------------


indir.cuts.other <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/count_tables_from_bins/binsize_10000")
indir.cuts.k27me3 <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/bams_split_by_cluster/bams_merged_by_cluster/counts_tables_10000")

pbulk.bin.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  indir.cuts <- ifelse(jmark == "H3K27me3", indir.cuts.k27me3, indir.cuts.other)
  infs.cuts <- list.files(path = indir.cuts, pattern = paste0(jmark, ".*.txt"), all.files = TRUE, full.names = TRUE)
  
  dat.csv.lst <- lapply(infs.cuts, function(inf.tmp){
    print(inf.tmp)
    # mat.tmp <- ReadMatTSSFormat(inf.tmp, as.sparse = TRUE, add.coord = TRUE, sort.rnames = FALSE)
    mat.tmp <- ReadMatSlideWinFormat(inf.tmp, as.sparse = TRUE, sort.rnames = FALSE, add.chromo = TRUE)
    rownames(mat.tmp) <- paste("chr", rownames(mat.tmp), sep = "")
    return(mat.tmp)
  })
  rnames.all <- sort(unique(unlist(lapply(dat.csv.lst, function(jmat) rownames(jmat)))))
  mat.tmp <- cbind.fill.lst(dat.csv.lst, all.rnames = rnames.all)
  
  dat.metas.filt.tmp <- dat.metas.filt[[jmark]]
  norm.factors <- norm.factors.lst[[jmark]]
  
  cells.lst <- split(x = dat.metas.filt.tmp$cell, f = dat.metas.filt.tmp$cluster)
  
  pbulk.tmp.lst <- SumAcrossClusters(mat.tmp, cells.lst)
  pbulk.mat <- do.call(cbind, pbulk.tmp.lst)
  
  
  # make colnames alphabetiacl
  pbulk.mat <- pbulk.mat[, sort(colnames(pbulk.mat))]
  
  pbulk.mat.norm <- sweep(pbulk.mat, MARGIN = 2, STATS = norm.factors$spikein_cuts, FUN = "/")
  return(list(pbulk.mat = pbulk.mat, pbulk.mat.norm = pbulk.mat.norm, norm.factors = norm.factors))
  
  return(pbulk.mat)
})


# Do different normalizations  --------------------------------------------




# Write outputs -----------------------------------------------------------



# write raw 
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pseudobulk_tables_TSS_and_bins"

for (jmark in jmarks){
  print(jmark)
  fname1 <- paste0("raw_cuts_pseudobulk_TSS.new_batch_only.", jmark, ".", Sys.Date(), ".txt")
  fname2 <- paste0("spikeinnorm_cuts_pseudobulk_TSS.new_batch_only.", jmark, Sys.Date(), ".txt")
  outf1 <- file.path(outdir, fname1)
  outf2 <- file.path(outdir, fname2)
  mat.tmp1 <- as.data.frame(pbulk.tss.lst[[jmark]]$pbulk.mat)
  mat.tmp2 <- as.data.frame(pbulk.tss.lst[[jmark]]$pbulk.mat.norm)
  fwrite(x = mat.tmp1, file = outf1, sep = "\t", quote = FALSE, row.names = TRUE)  
  fwrite(x = mat.tmp2, file = outf2, sep = "\t", quote = FALSE, row.names = TRUE)  
}

for (jmark in jmarks){
  print(jmark)
  fname1 <- paste0("raw_cuts_pseudobulk_10000_bins.new_batch_only.", jmark, ".", Sys.Date(), ".txt")
  fname2 <- paste0("spikeinnorm_cuts_pseudobulk_10000_bins.new_batch_only.", jmark, ".", Sys.Date(), ".txt")
  outf1 <- file.path(outdir, fname1)
  outf2 <- file.path(outdir, fname2)
  mat.tmp1 <- as.data.frame(pbulk.bin.lst[[jmark]]$pbulk.mat)
  mat.tmp2 <- as.data.frame(pbulk.bin.lst[[jmark]]$pbulk.mat.norm)
  fwrite(x = mat.tmp1, file = outf1, sep = "\t", quote = FALSE, row.names = TRUE)  
  fwrite(x = mat.tmp2, file = outf2, sep = "\t", quote = FALSE, row.names = TRUE)  
}

# write norm factors
for (jmark in jmarks){
  fname1 <- paste0("spikein_norm_factors.", jmark, ".txt")
  outf1 <- file.path(outdir, fname1)
  fwrite(norm.factors.lst[[jmark]], file = outf1, sep = "\t")
}


