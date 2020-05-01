# Jake Yeung
# Date of Creation: 2020-04-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/1-plot_TSS_signal_onto_UMAP.R
# Plot TSS signal onto UMAP 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(topicmodels)

library(hash)
library(igraph)
library(umap)


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jmark <- "H3K27me3"
jmark <- "H3K9me3"
jmark <- "H3K4me3"
jmark <- "H3K4me1"

jdist <- 10000
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


# Load TSS signal ---------------------------------------------------------

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/TSS_signal"
for (jmark in jmarks){
  
  pdfout <- file.path(outdir, paste0("TSS_signal_on_stringent_umap.", jmark, ".pdf"))
  
  pdf(pdfout, width = 1020/72, height = 815/72, useDingbats = FALSE)
  
  
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.TSS.winsize_", jdist, "/PZ-ChIC-ZF_", jmark, "_2020-04-07.countTable.TSS.csv")
  assertthat::assert_that(file.exists(inf))
  
  
  # indir.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZF_AllMerged.winsize_50000"
  # inf.lda <- file.path(indir.lda, paste0("lda_outputs.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.Robj"))
  # inf.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZF_AllMerged.winsize_50000.imputevarfilt.lessstringent.remove_eryth/lda_outputs.counts_table_var_filt.H3K4me3.imputevar_2.K-30.binarize.FALSE/ldaOut.counts_table_var_filt.H3K4me3.imputevar_2.K-30.Robj"
  inmain.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZF_AllMerged.winsize_50000.imputevarfilt.lessstringent"
  inf.lda <- list.files(inmain.lda, pattern = jmark, all.files = TRUE, full.names = TRUE, recursive = TRUE)
  
  count.tss <- ReadMatTSSFormat(inf)
  count.tss.filt <- CollapseRowsByGene(count.tss, as.long = FALSE, track.kept.gene = TRUE)
  
  
  # count.tss.filt <- BinarizeMatrix(count.tss.filt)
      
  # Load LDA ----------------------------------------------------------------
  
  load(inf.lda, v=T)
  
  
  
  # Do UMAP  ----------------------------------------------------------------
  
  jchromos <- paste("chr", seq(25), sep = "")
  
  tm.result <- posterior(out.lda)
  topics.mat <- tm.result$topics
  
  dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
  
  dat.umap <- DoUmapAndLouvain(topics.mat, jsettings)
  
  dat.var <- CalculateVarAll(dat.impute.log, jchromos)
  
  
  
  
  # Plot TSS signal onto UMAP  ----------------------------------------------
  
  # tss.cuts are not filtered by cells... so do left join 
  tss.cuts.dat <- data.frame(cell = colnames(count.tss.filt), tss.cuts = colSums(count.tss.filt), stringsAsFactors = FALSE)
  total.cuts.dat <- data.frame(cell = colnames(count.mat), total.cuts = colSums(count.mat), stringsAsFactors = FALSE)
  
  merge.cuts.dat <- left_join(total.cuts.dat, tss.cuts.dat) %>%
    mutate(frac.tss = tss.cuts / total.cuts) 
  
  
  dat.umap.merge <- left_join(dat.umap, merge.cuts.dat) %>%
    left_join(., dat.var) %>%
    rowwise() %>%
    mutate(plate = ClipLast(cell, jsep = "_"))
  
  m <- ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = frac.tss)) + geom_point() +  
    scale_color_viridis_c(direction = 1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
    facet_wrap(~plate) + ggtitle(jmark)
  print(m)
  
  m <- ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = frac.tss)) + geom_point() +  
    scale_color_viridis_c(direction = 1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
    ggtitle(jmark)
  print(m)
  
  ngenes <- nrow(count.tss.filt)
  m <- ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = tss.cuts / ngenes)) + geom_point() +  
    scale_color_viridis_c(direction = 1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
    facet_wrap(~plate) + ggtitle(jmark)
  print(m)
  
  
  m <- ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() +  
    scale_color_viridis_c(direction = 1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~plate) + ggtitle(jmark)
  print(m)
  
  m <- ggplot(dat.umap.merge, aes(x = tss.cuts, y = cell.var.within.sum.norm)) + geom_point() +  
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~plate) + ggtitle(jmark)
  print(m)
  dev.off()
  
}


