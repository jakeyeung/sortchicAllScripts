# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/5-LDA_downstream.R
# 


rm(list=ls())

library(heatmap3)

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

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K9me3"); names(jmarks) <- jmarks

keeptop <- 150
high.in.k9 <- FALSE

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

# add jrep2 for batch correction?
dat.metas$H3K4me1$jrep2 <- sapply(dat.metas$H3K4me1$jrep, function(x) ifelse(x == "rep1old", "zold", "anew"))  # batch2 is better than batch1
dat.metas$H3K9me3$jrep2 <- sapply(dat.metas$H3K9me3$jrep, function(x) ifelse(x != "rep1old", "zold", "anew"))  # batch1 is better than batch2

# Select bins  ------------------------------------------------------------

load("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_H3K4me1_H3K9me3_analysis/de_analysis_H3K4me1_H3K9me3.RData", v=T)

k9.bins <- names(which(pvals.lst2 < 1e-10))
k4.bins <- names(which(pvals.lst1 < 1e-100))


# Load RData  -------------------------------------------------------------

inf.mat.adj <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_H3K4me1_H3K9me3_analysis/batch_corrected_imputed_values_H3K4me1_H3K9me3.mat.namesfix.RData"
load(inf.mat.adj, v=T)


mat.adj.lst <- lapply(mat.adj.lst, function(jmat){
  rownames(jmat) <- jmat$rname
  jmat$rname <- NULL
  return(jmat)
})


# Select bins  ------------------------------------------------------------

k9.bins <- which(pvals.lst2 < 1e-10)
k4.bins <- which(pvals.lst1 < 1e-100)

k9.bins.names <- names(k9.bins)
k4.bins.names <- names(k4.bins)

params.dat2.wide <- reshape2::dcast(subset(params.dat2.all, bin %in% k9.bins.names), formula = bin ~ param, value.var = "estimate2") %>%
  rowwise() %>%
  mutate(bcell.effect = ClusterBcells.Estimate - mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate)),
         eryth.effect = ClusterEryths.Estimate - mean(c(ClusterBcells.Estimate, ClusterGranulocytes.Estimate)),
         granu.effect = ClusterGranulocytes.Estimate - mean(c(ClusterBcells.Estimate, ClusterEryths.Estimate)),
         mean.effect = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate, ClusterBcells.Estimate)))


if (high.in.k9){
  
  
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

# Plot heatmaps  ----------------------------------------------------------


bins.keep <- c(jbins.eryth, jbins.bcell, jbins.granu, jbins.hspcs)

bins.keep.lst <- list("Eryths" = jbins.eryth,
                      "Bcells" = jbins.bcell,
                      "Granulocytes" = jbins.granu,
                      "HSPCs" = jbins.hspcs)
bnames <- names(bins.keep.lst); names(bnames) <- bnames


bins.common <- intersect(rownames(mat.adj.lst$H3K4me1), rownames(mat.adj.lst$H3K9me3))

bins.keep.common <- bins.keep[bins.keep %in% bins.common]

ctypes <- c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "pDCs", "DCs", "HSPCs")

jmeta1 <- dat.metas$H3K9me3
jmeta1$cluster <- factor(jmeta1$cluster, levels = c("Eryths", "Bcells", "Granulocytes", "HSPCs"))
jmeta1 <- jmeta1 %>% arrange(cluster, jrep)
cells.keep1 <- jmeta1$cell


jmeta2 <- dat.metas$H3K4me1
jmeta2$cluster <- factor(jmeta2$cluster, levels = ctypes)
jmeta2 <- jmeta2 %>% arrange(cluster, jrep)
cells.keep2 <- jmeta2$cell

# labrows = dat.genes.sub.join$genesymbol
# colsidecolors = dat.meta.sub$colorcode
# rowsidecolors = dat.genes.sub.join$colorcode
# heatmap3::heatmap3(mat.adj.tmp, cexRow = 0.08, labRow = labrows, labCol = "", Rowv = NA, Colv = NA, ColSideColors = colsidecolors, scale = "row", RowSideColors = rowsidecolors, revC = TRUE, main = paste0(jmark.test))

jmat1 <- mat.adj.lst$H3K9me3[bins.keep.common, cells.keep1]
jmat2 <- mat.adj.lst$H3K4me1[bins.keep.common, cells.keep2]

jsum <- data.frame(cell = cells.keep1, log2exprsadj.mean = colMeans(mat.adj.lst$H3K9me3[jbins.eryth, cells.keep1]), stringsAsFactors = FALSE) %>%
  left_join(., jmeta1)

# ggplot(jsum, aes(x = cluster, y = log2exprsadj.mean)) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   geom_point() +
#   geom_boxplot()
# 
# plot(mat.adj.lst$H3K9me3[jsort.eryth$bin, cells.keep1])

ctype2col <- hash::hash(jmeta2$cluster, jmeta2$colorcode)

# annotate genes
names(bins.keep) <- c(rep("Eryths", 150), rep("Bcells", 150), rep("Granulocytes", 150), rep("HSPCs", 150))
# bins.keep <- c(jsort.eryth$bin, jsort.bcell$bin, jsort.granu$bin, jsort.hspcs$bin)
colsvec <- sapply(names(bins.keep), function(x) AssignHash(x, jhash = ctype2col, null.fill = NA))

bin2col <- hash::hash(bins.keep, colsvec)

# color rows

jmat2 <- t(apply(jmat2, 1, function(jrow) Winsorize(jrow, probs = c(0.01, 0.99))))
jmat2 <- apply(jmat2, 2, function(jcol) Winsorize(jcol, probs = c(0.01, 0.99)))

jmat1 <- t(apply(jmat1, 1, function(jrow) Winsorize(jrow, probs = c(0.01, 0.99))))
jmat1 <- apply(jmat1, 2, function(jcol) Winsorize(jcol, probs = c(0.01, 0.99)))

outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_differential_expression_outputs/heatmap_k9me3_k4me1_signif_bins_k9.highink9_", high.in.k9, ".", Sys.Date(), ".pdf")
print("Making heatmaps")
pdf(outpdf, useDingbats = FALSE)




  heatmap3::heatmap3(jmat1, Rowv = NA, Colv = NA, scale = "row", ColSideColors = jmeta1$clustercol, RowSideColors = sapply(bins.keep.common, AssignHash, jhash = bin2col),  revC = TRUE, main = "H3K9me3 50kb bins", margins = c(5, 8))
  heatmap3::heatmap3(jmat2, Rowv = NA, Colv = NA, scale = "row",  ColSideColors = jmeta2$clustercol, RowSideColors = sapply(bins.keep.common, AssignHash, jhash = bin2col), revC = TRUE, main = "H3K4me1 50kb bins", margins = c(5, 8))
dev.off()

