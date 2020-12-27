# Jake Yeung
# Date of Creation: 2020-12-18
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/9-check_k9_bins_on_other_marks.R
# Check k9 bins on other marks


rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(JFuncs)
library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

library(topicmodels)
library(DescTools)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

hubprefix <- "/home/jyeung/hub_oudenaarden"

# jmarks <- c("H3K4me1", "H3K9me3"); names(jmarks) <- jmarks
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks



keeptop <- 150
low.in.k9 <- FALSE
# low.in.k9 <- TRUE
# outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_differential_expression_outputs/heatmap_k9me3_k4me1_signif_bins_k9.highink9_", low.in.k9, ".", Sys.Date(), ".WithLogFCmaps.pdf")

# Load LDA outputs --------------------------------------------------------

out.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj"))
  load(inf.lda, v=T)
  tm.result <- posterior(out.lda)
  return(list(tm.result = tm.result, count.mat = count.mat))
})


# Load meta data  ---------------------------------------------------------

dat.metas <- lapply(jmarks, function(jmark){
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/umaps_final/umaps_final_get_marks.", jmark, ".metadata.2020-12-17.txt")
  fread(inf)
})

# # add jrep2 for batch correction?
# dat.metas$H3K4me1$jrep2 <- sapply(dat.metas$H3K4me1$jrep, function(x) ifelse(x == "rep1old", "zold", "anew"))  # batch2 is better than batch1
# dat.metas$H3K9me3$jrep2 <- sapply(dat.metas$H3K9me3$jrep, function(x) ifelse(x != "rep1old", "zold", "anew"))  # batch1 is better than batch2

# Select bins  ------------------------------------------------------------

load("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_H3K4me1_H3K9me3_analysis/de_analysis_H3K4me1_H3K9me3.RData", v=T)

k9.bins <- names(which(pvals.lst2 < 1e-10))



# Pick bins ---------------------------------------------------------------

# check 


k9.bins <- which(pvals.lst2 < 1e-10)

k9.bins.names <- names(k9.bins)

params.dat2.wide <- data.table::dcast(subset(params.dat2.all, bin %in% k9.bins.names), formula = bin ~ param, value.var = "estimate2") %>%
  rowwise() %>%
  mutate(ClusterBcells.Estimate = ClusterBcells.Estimate / log(2),
         ClusterGranulocytes.Estimate = ClusterGranulocytes.Estimate / log(2),
         ClusterEryths.Estimate = ClusterEryths.Estimate / log(2),
         ClusterBcells.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterBcells.Estimate, ClusterGranulocytes.Estimate)),
         ClusterEryths.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate)),
         ClusterBcells.Estimate_ClusterEryths.Estimate =  mean(c(ClusterBcells.Estimate, ClusterEryths.Estimate)),
         bcell.effect = ClusterBcells.Estimate - ClusterEryths.Estimate_ClusterGranulocytes.Estimate,
         eryth.effect = ClusterEryths.Estimate - ClusterBcells.Estimate_ClusterGranulocytes.Estimate,
         granu.effect = ClusterGranulocytes.Estimate - ClusterBcells.Estimate_ClusterEryths.Estimate,
         mean.effect = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate, ClusterBcells.Estimate)))


if (low.in.k9){
  jsort.hspcs <- params.dat2.wide %>%
    group_by(bin) %>%
    # arrange(mean.effect)
    arrange(desc(mean.effect))
  jbins.hspcs <- jsort.hspcs$bin[1:keeptop]
  
  jsort.bcell <- params.dat2.wide %>%
    group_by(bin) %>%
    # arrange(desc(bcell.effect)) 
    arrange(bcell.effect)
  jbins.bcell <- jsort.bcell$bin[1:keeptop]
  
  jsort.granu <- params.dat2.wide %>%
    group_by(bin) %>%
    # arrange(desc(granu.effect))
    arrange(granu.effect)
  jbins.granu <- jsort.granu$bin[1:keeptop]
  
  jsort.eryth <- params.dat2.wide %>%
    group_by(bin) %>%
    # arrange(desceryth.effect)) 
    arrange(eryth.effect)
  jbins.eryth <- jsort.eryth$bin[1:keeptop]
} else {
  jsort.hspcs <- params.dat2.wide %>%
    group_by(bin) %>%
    arrange(mean.effect)
  jbins.hspcs <- jsort.hspcs$bin[1:keeptop]
  
  jsort.bcell <- params.dat2.wide %>%
    group_by(bin) %>%
    arrange(desc(bcell.effect))
  jbins.bcell <- jsort.bcell$bin[1:keeptop]
  
  jsort.granu <- params.dat2.wide %>%
    group_by(bin) %>%
    arrange(desc(granu.effect))
  jbins.granu <- jsort.granu$bin[1:keeptop]
  
  jsort.eryth <- params.dat2.wide %>%
    group_by(bin) %>%
    arrange(desc(eryth.effect))
  jbins.eryth <- jsort.eryth$bin[1:keeptop]
}


# Check raw cuts on umap  -------------------------------------------------

count.mat.lst <- lapply(jmarks, function(jmark){
  jmat <- out.lst[[jmark]]$count.mat
  jmat <- sweep(x = jmat, MARGIN = 2, STATS = colSums(jmat), FUN = "/")
})

# take mean across rows, log2 transform


bins.keep <- c(jbins.eryth, jbins.bcell, jbins.granu, jbins.hspcs)

bins.keep.lst <- list("Eryths" = jbins.eryth,
                      "Bcells" = jbins.bcell,
                      "Granulocytes" = jbins.granu,
                      "HSPCs" = jbins.hspcs)
bnames <- names(bins.keep.lst); names(bnames) <- bnames

dat.cuts.lst.lst <- lapply(bnames, function(bname){
  jbins.keep <- bins.keep.lst[[bname]]
  dat.cuts.lst <- lapply(jmarks, function(jmark){
    x <- count.mat.lst[[jmark]]
    rkeep <- rownames(x) %in% jbins.keep
    dat.cuts <- data.frame(log2cuts = log2(colSums(x[rkeep, ]) + 1), cell = colnames(x)) %>%
      left_join(., dat.metas[[jmark]])
    dat.cuts$bname <- bname
    return(dat.cuts)
  })
})


# pdf("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_differential_expression_outputs/UMAPs_raw_cuts_K9_bins")
outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_differential_expression_outputs/UMAPs_raw_cuts_K9_bins_", low.in.k9, ".", Sys.Date(), ".WithLogFCmaps.pdf")

# plot reference 
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

pdf(file = outpdf, useDingbats = FALSE)
mlst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.metas[[jmark]], aes(x = umap1, y = umap2, color = clustercol)) + 
    geom_point(size = 0.5, alpha = 0.2) + 
    # scale_color_manual(values = cbPalette) + 
    scale_color_identity() + 
    theme_bw() + 
    ggtitle(jmark) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m)
})
JFuncs::multiplot(mlst[[1]], mlst[[2]], mlst[[3]], mlst[[4]], cols = 4)
for (bname in bnames){
  mlst <- lapply(jmarks, function(jmark){
  # for (jmark in jmarks){
    m <- ggplot(dat.cuts.lst.lst[[bname]][[jmark]], aes(x = umap1, y = umap2, color = Winsorize(log2cuts, probs = c(0.05, 0.95)))) + 
      geom_point(size = 0.5) + 
      scale_color_viridis_c() + 
      theme_bw() + 
      # ggtitle(paste(bname, jmark, paste0("Low in K9:", low.in.k9))) + 
      ggtitle(bname, paste(jmark, paste0(low.in.k9))) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    # print(m)
    return(m)
  })
  JFuncs::multiplot(mlst[[1]], mlst[[2]], mlst[[3]], mlst[[4]], cols = 4)
  # }
}
dev.off()
