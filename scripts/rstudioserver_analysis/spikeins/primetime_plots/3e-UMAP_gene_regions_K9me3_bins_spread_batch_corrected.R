# Jake Yeung
# Date of Creation: 2020-12-28
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/3e-UMAP_gene_regions_K9me3_bins_spread_batch_corrected.R
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

outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/batch_corrected/umap_gene_sets.H3K9me3_bins.", Sys.Date(), ".again.pdf")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


# Load meta data  ---------------------------------------------------------


indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28")

dat.metas <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.meta <- file.path(indir.meta, paste0("metadata_batch_corrected.", jmark, ".2020-12-28.txt"))
  print(inf.meta)
  dat.meta <- fread(inf.meta)
})


# Load mats  --------------------------------------------------------------


inf.mat.adj <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_H3K4me1_H3K9me3_analysis/batch_corrected_imputed_values.bins.all_marks.mat.namesfix.2020-12-20.H3K27me3rep2rep3reseq.RData"
load(inf.mat.adj, v=T)

mat.adj.lst <- lapply(mat.adj.lst, function(jmat){
  rownames(jmat) <- jmat$rname
  jmat$rname <- NULL
  return(jmat)
})


# Get raw cuts  -----------------------------------------------------------


out.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  if (jmark != "H3K27me3"){
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj"))
  } else {
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.Robj"))
  }
  load(inf.lda, v=T)
  tm.result <- posterior(out.lda)
  return(list(tm.result = tm.result, count.mat = count.mat))
})

count.mat.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jout <- out.lst[[jmark]]
  sweep(jout$count.mat, MARGIN = 2, STATS = colSums(jout$count.mat), FUN = "/")
})


# Get DE outputs ----------------------------------------------------------


jlow.in.k9 <- TRUE
jkeeptop <- 150

jfits.lst.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jinf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.RData"))
  load(jinf, v=T)
  return(jfits.lst)
})

params.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jfits.lst <- jfits.lst.lst[[jmark]]
  jnames <- names(jfits.lst); names(jnames) <- jnames
  params.dat.all <- lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xkeep <- grepl("^Cluster.*.Estimate$", x = names(x))
    jparams <- x[xkeep]
    data.frame(bin = jname, param = names(jparams), estimate = unlist(jparams), mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  if (jmark == "H3K9me3"){
    params.dat.all <- params.dat.all %>% 
      mutate(param = gsub("Eryth", "Eryths", param),
             param = gsub("Lymphoid", "Bcells", param))
  }
  # make params more readable
  params.dat.all$ctype <- params.dat.all$param
  params.dat.all$ctype <- gsub("Cluster", "", params.dat.all$ctype)
  params.dat.all$ctype <- gsub(".Estimate", "", params.dat.all$ctype)
  return(params.dat.all)
})

pvals.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jfits.lst <- jfits.lst.lst[[jmark]]
  jnames <- names(jfits.lst); names(jnames) <- jnames
  lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xvec <- x$pval
    data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
})

pvals.lst2 <- lapply(jfits.lst.lst$H3K9me3, function(x) x$pval)
k9.bins <- which(pvals.lst2 < 1e-10)


pval.k9.sub <- subset(pvals.lst$H3K9me3, pval < 1e-10) %>%
  arrange(desc(pval))

k9.bins.names <- unique(pval.k9.sub$bin)
ctypes.keep <- c("Eryths", "Bcells", "Granulocytes")
params.keep <- paste("Cluster", ctypes.keep, ".Estimate", sep = "")


params.dat.wide.lst <- lapply(jmarks, function(jmark){
  jsub <- subset(params.lst[[jmark]], bin %in% k9.bins.names & param %in% params.keep) %>%
    group_by(bin) %>% filter(max(abs(estimate)) < 5)
  jdat <- GetParamsWideFormat(jsub, jvalue.var = "estimate")
  # # keep only effect cnames ( do this later maybe ? )
  # cnames.keep.i <- grep("effect$", colnames(jdat))
  # cnames.new <- paste(colnames(jdat)[cnames.keep.i], jmark, sep = "_")
  # colnames(jdat)[cnames.keep.i] <- cnames.new
  # cnames.keep.bin.i <- grep("bin", colnames(jdat))
  # cnames.keep.merged.i <- c(cnames.keep.bin.i, cnames.keep.i)
  # jdat.filt <- jdat[, cnames.keep.merged.i]
  return(jdat)
})

bins.keep.lst <- GetK9CelltypeBins(params.dat.wide.lst$H3K9me3, low.in.k9 = jlow.in.k9, keeptop = jkeeptop)
bnames <- names(bins.keep.lst); names(bnames) <- bnames

jbins.eryth <- bins.keep.lst[["Eryths"]]
jbins.bcell <- bins.keep.lst[["Bcells"]]
jbins.granu <- bins.keep.lst[["Granulocytes"]]
jbins.hspcs <- bins.keep.lst[["HSPCs"]]

bins.keep <- c(jbins.eryth, jbins.bcell, jbins.granu, jbins.hspcs)


bins.common <- Reduce(f = intersect, x = lapply(mat.adj.lst, rownames))
bins.keep.common <- bins.keep[bins.keep %in% bins.common]




# Get mean across HSPC-specific bins and others  -------------------------------------

# jset <- "HSPCs"
# # jmarktest <- "H3K4me1"
# jmarktest <- "H3K9me3"


pdf(outpdf, useDingbats = FALSE)

jmarktest.vec <- c("H3K4me1", "H3K9me3", "H3K4me3", "H3K27me3")
names(jmarktest.vec) <- jmarktest.vec

for (jset in bnames){
  for (jmarktest in jmarktest.vec){
    jmat.raw <- count.mat.lst[[jmarktest]]
    jbins.jset <- bins.keep.common[bins.keep.common %in% bins.keep.lst[[jset]]]
    jbins.jset2 <- rownames(jmat.raw)[rownames(count.mat.lst[[jmarktest]]) %in% jbins.jset]
    dat.imputed.bins <- data.frame(cell = colnames(mat.adj.lst[[jmarktest]]), exprs.bins = colMeans(mat.adj.lst[[jmarktest]][jbins.jset, ]), stringsAsFactors = FALSE) %>%
      left_join(., dat.metas[[jmarktest]])
    dat.raw.bins <- data.frame(cell = colnames(count.mat.lst[[jmarktest]]), log2cuts = log2(colMeans(count.mat.lst[[jmarktest]][jbins.jset2, ]) + 1), stringsAsFactors = FALSE) %>%
      left_join(., dat.metas[[jmarktest]])
    
    m <- ggplot(dat.imputed.bins, aes(x = umap1, y = umap2, color = exprs.bins)) + 
      geom_point() + 
      theme_minimal(2) + 
      ggtitle(paste(jmarktest, jset, "imputed")) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
      scale_color_viridis_c() 
    print(m)
    m <- ggplot(dat.raw.bins, aes(x = umap1, y = umap2, color = Winsorize(log2cuts, probs = c(0.05, 0.95)))) + 
    # m <- ggplot(dat.raw.bins, aes(x = umap1, y = umap2, color = aes(x = umap1, y = umap2, color = log2cuts))) + 
      geom_point() + 
      theme_minimal(2) + 
      ggtitle(paste(jmarktest, jset, "raw")) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
      scale_color_viridis_c() 
    print(m)
  }
}

dev.off()
