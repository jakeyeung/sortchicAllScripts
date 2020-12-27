# Jake Yeung
# Date of Creation: 2020-12-18
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/8-summarize_DE_analysis_K4me1_K9me3.R
# 


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

jmarks <- c("H3K4me1", "H3K9me3"); names(jmarks) <- jmarks


keeptop <- 150
# low.in.k9 <- TRUE
low.in.k9 <- FALSE
outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_differential_expression_outputs/heatmap_k9me3_k4me1_signif_bins_k9.lowink9_", low.in.k9, ".", Sys.Date(), ".WithLogFCmaps.pdf")

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



# Show how we picked bins  ------------------------------------------------

print(head(params.dat2.all))

params.keep <- c("ClusterBcells.Estimate", "ClusterEryths.Estimate", "ClusterGranulocytes.Estimate")
dat.merge <- left_join(params.dat1.all %>% filter(param %in% params.keep), 
                       params.dat2.all %>% filter(param %in% params.keep), 
                       by = c("bin", "param"))
dat.merge <- dat.merge %>%
  rowwise() %>%
  mutate(celltype = gsub("^Cluster", "", x = param),
         celltype = gsub(".Estimate$", "", x = celltype))
  

ggplot(dat.merge %>% filter(abs(estimate2) < 5 & bin %in% k9.bins), aes(x = estimate1, y = estimate2)) + 
  geom_point(alpha = 0.25) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_density_2d() + 
  theme_bw() + 
  facet_wrap(~param) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Pick bins ---------------------------------------------------------------

# check 


k9.bins <- which(pvals.lst2 < 1e-10)
k4.bins <- which(pvals.lst1 < 1e-100)

k9.bins.names <- names(k9.bins)
k4.bins.names <- names(k4.bins)

params.dat2.wide <- GetParamsWideFormat(subset(params.dat2.all, bin %in% k9.bins.names), jvalue.var = "estimate2")
# params.dat2.wide <- data.table::dcast(subset(params.dat2.all, bin %in% k9.bins.names), formula = bin ~ param, value.var = "estimate2") %>%
#   rowwise() %>%
#   mutate(ClusterBcells.Estimate = ClusterBcells.Estimate / log(2),
#          ClusterGranulocytes.Estimate = ClusterGranulocytes.Estimate / log(2),
#          ClusterEryths.Estimate = ClusterEryths.Estimate / log(2),
#          ClusterBcells.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterBcells.Estimate, ClusterGranulocytes.Estimate)),
#          ClusterEryths.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate)),
#          ClusterBcells.Estimate_ClusterEryths.Estimate =  mean(c(ClusterBcells.Estimate, ClusterEryths.Estimate)),
#          Bcells.effect = ClusterBcells.Estimate - ClusterEryths.Estimate_ClusterGranulocytes.Estimate,
#          Eryths.effect = ClusterEryths.Estimate - ClusterBcells.Estimate_ClusterGranulocytes.Estimate,
#          Granulocytes.effect = ClusterGranulocytes.Estimate - ClusterBcells.Estimate_ClusterEryths.Estimate,
#          HSPCs.effect = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate, ClusterBcells.Estimate)))

bins.keep.lst <- GetK9CelltypeBins(params.dat2.wide, low.in.k9 = low.in.k9, keeptop = keeptop)

# if (low.in.k9){
#   jsort.hspcs <- params.dat.wide %>%
#     group_by(bin) %>%
#     # arrange(HSPCs.effect)
#     arrange(desc(HSPCs.effect))
#   jbins.hspcs <- jsort.hspcs$bin[1:keeptop]
#   
#   jsort.bcell <- params.dat.wide %>%
#     group_by(bin) %>%
#     # arrange(desc(Bcells.effect)) 
#     arrange(Bcells.effect)
#   jbins.bcell <- jsort.bcell$bin[1:keeptop]
#   
#   jsort.granu <- params.dat.wide %>%
#     group_by(bin) %>%
#     # arrange(desc(Granulocytes.effect))
#     arrange(Granulocytes.effect)
#   jbins.granu <- jsort.granu$bin[1:keeptop]
#   
#   jsort.eryth <- params.dat.wide %>%
#     group_by(bin) %>%
#     # arrange(descEryths.effect)) 
#     arrange(Eryths.effect)
#   jbins.eryth <- jsort.eryth$bin[1:keeptop]
# } else {
#   jsort.hspcs <- params.dat.wide %>%
#     group_by(bin) %>%
#     arrange(HSPCs.effect)
#   jbins.hspcs <- jsort.hspcs$bin[1:keeptop]
#   
#   jsort.bcell <- params.dat.wide %>%
#     group_by(bin) %>%
#     arrange(desc(Bcells.effect))
#   jbins.bcell <- jsort.bcell$bin[1:keeptop]
#   
#   jsort.granu <- params.dat.wide %>%
#     group_by(bin) %>%
#     arrange(desc(Granulocytes.effect))
#   jbins.granu <- jsort.granu$bin[1:keeptop]
#   
#   jsort.eryth <- params.dat.wide %>%
#     group_by(bin) %>%
#     arrange(desc(Eryths.effect))
#   jbins.eryth <- jsort.eryth$bin[1:keeptop]
# }




# Check raw cuts  ---------------------------------------------------------


print("Making heatmaps")
pdf(outpdf, useDingbats = FALSE)


jbins.eryth <- bins.keep.lst[["Eryths"]]
jbins.bcell <- bins.keep.lst[["Bcells"]]
jbins.granu <- bins.keep.lst[["Granulocytes"]]
jbins.hspcs <- bins.keep.lst[["HSPCs"]]

bins.keep <- c(jbins.eryth, jbins.bcell, jbins.granu, jbins.hspcs)

# bins.keep.lst <- list("Eryths" = jbins.eryth,
#                       "Bcells" = jbins.bcell,
#                       "Granulocytes" = jbins.granu,
#                       "HSPCs" = jbins.hspcs)
bnames <- names(bins.keep.lst); names(bnames) <- bnames

bins.common <- intersect(rownames(mat.adj.lst$H3K4me1), rownames(mat.adj.lst$H3K9me3))

bins.keep.common <- bins.keep[bins.keep %in% bins.common]


dat.merge.filt <- dat.merge %>% filter(abs(estimate2) < 5 & bin %in% k9.bins.names)
params.dat2.wide.filt <- subset(params.dat2.wide, abs(ClusterEryths.Estimate) < 5)
ctypes.sub <- c("Bcells", "Eryths", "Granulocytes")
cnames.all <- paste("Cluster", ctypes.sub, ".Estimate", sep = "")
assertthat::assert_that(all(cnames.all %in% colnames(params.dat2.wide.filt)))

m <- ggplot(dat.merge.filt,
            aes(x = estimate1 / log(2), y = estimate2 / log(2))) + 
  geom_point(alpha = 0.25) + 
  ggtitle("50kb bins, K9me3 pval < 10^-10")  + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  xlab("log2FC H3K4me1 relative to HSPCs") + 
  ylab("log2FC H3K9me3 relative to HSPCs") + 
  geom_density_2d(inherit.aes = FALSE, mapping = aes(x = estimate1, y = estimate2)) + 
  theme_bw() + 
  facet_wrap(~param) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

# show bins in 2d plot 
m.lst <- lapply(bnames, function(bname){
  m1 <- ggplot(dat.merge.filt %>% 
                mutate(tag = bin %in% bins.keep.lst[[bname]]) %>% 
                arrange(tag), 
              aes(x = estimate1 / log(2), y = estimate2 / log(2), color = tag)) + 
    geom_point(alpha = 0.25) + 
    ggtitle(paste("50kb bins, K9me3 pval < 10^-10, colored for celltypespecific bins", bname))  + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    xlab("log2FC H3K4me1 relative to HSPCs") + 
    ylab("log2FC H3K9me3 relative to HSPCs") + 
    # geom_density_2d(inherit.aes = FALSE, mapping = aes(x = estimate1, y = estimate2), color = "black") + 
    theme_bw() + 
    facet_wrap(~param) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  # show k9me3 specific estimates
  if (bname != "HSPCs"){
    cname.keep <- paste0("Cluster", bname, ".Estimate")
  } else {
    cname.keep <- "ClusterEryths.Estimate"
  }
  print(m1)
  # plot combined estimate
  cnames.sub <- cnames.all[!cnames.all %in% cname.keep]
  cnames.sub.paste <- paste(cnames.sub, collapse = "_")
  m2.combined <- ggplot(params.dat2.wide.filt %>% mutate(tag = bin %in% bins.keep.lst[[bname]]), 
               aes_string(x = cname.keep, y = cnames.sub.paste, color = "tag")) + 
    geom_point(alpha = 0.25) + 
    ggtitle(paste("H3K9me3 log2FC estimates", bname)) +
    geom_hline(yintercept = 0, linetype = "dotted") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m2.combined)
  
  for (cname.other in cnames.sub){
    m2 <- ggplot(params.dat2.wide.filt %>% mutate(tag = bin %in% bins.keep.lst[[bname]]), 
                 aes_string(x = cname.keep, y = cname.other, color = "tag")) + 
      geom_point(alpha = 0.25) + 
      ggtitle(paste("H3K9me3 log2FC estimates", bname)) +
      geom_hline(yintercept = 0, linetype = "dotted") + 
      geom_vline(xintercept = 0, linetype = "dotted") + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m2)
    
  }
  return(m1)
})

# print(m.lst#)



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

# jsum <- data.frame(cell = cells.keep1, log2exprsadj.mean = colMeans(mat.adj.lst$H3K9me3[jbins.eryth, cells.keep1]), stringsAsFactors = FALSE) %>%
#   left_join(., jmeta1)
# 
# ggplot(jsum, aes(x = cluster, y = log2exprsadj.mean)) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   geom_point() +
#   geom_boxplot()

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

heatmap3::heatmap3(jmat1, Rowv = NA, Colv = NA, scale = "row", ColSideColors = jmeta1$clustercol, RowSideColors = sapply(bins.keep.common, AssignHash, jhash = bin2col),  revC = TRUE, main = "H3K9me3 50kb bins")
heatmap3::heatmap3(jmat2, Rowv = NA, Colv = NA, scale = "row",  ColSideColors = jmeta2$clustercol, RowSideColors = sapply(bins.keep.common, AssignHash, jhash = bin2col), revC = TRUE, main = "H3K4me1 50kb bins")

dev.off()

