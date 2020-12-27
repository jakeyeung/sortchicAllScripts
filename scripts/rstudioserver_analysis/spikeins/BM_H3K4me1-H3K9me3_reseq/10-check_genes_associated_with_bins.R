# Jake Yeung
# Date of Creation: 2020-12-19
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/10-check_genes_associated_with_bins.R
# 

rm(list=ls())

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

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

hubprefix <- "/home/jyeung/hub_oudenaarden"

# jmarks <- c("H3K4me1", "H3K9me3"); names(jmarks) <- jmarks
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks



keeptop <- 150
low.in.k9 <- TRUE
# low.in.k9 <- FALSE
# outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_differential_expression_outputs/heatmap_k9me3_k4me1_signif_bins_k9.highink9_", low.in.k9, ".", Sys.Date(), ".WithLogFCmaps.pdf")

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_differential_expression_outputs"
outpdf <- file.path(outdir, paste0("correlate_k9_specific_bins.with_giladi.topn_", keeptop, ".LowInK9_", low.in.k9, ".", Sys.Date(), ".pdf"))
outtxt <- file.path(outdir, paste0("correlate_k9_specific_bins.topn_", keeptop, ".LowInK9_", low.in.k9, ".", Sys.Date(), ".txt"))


# Select bins  ------------------------------------------------------------

load("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_H3K4me1_H3K9me3_analysis/de_analysis_H3K4me1_H3K9me3.RData", v=T)

k9.bins <- names(which(pvals.lst2 < 1e-10))



# Pick bins ---------------------------------------------------------------

# check 

k9.bins <- which(pvals.lst2 < 1e-10)

k9.bins.names <- names(k9.bins)

params.dat2.wide <- GetParamsWideFormat(params.long.filt = subset(params.dat2.all, bin %in% k9.bins.names), jvalue.var = "estimate2")

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


bins.keep.lst <- GetK9CelltypeBins(params.dat.wide = params.dat2.wide, low.in.k9 = low.in.k9, keeptop = keeptop)

# add pval
dat.pvals.k9 <- data.frame(pvals = unlist(pvals.lst2), bin = names(unlist(pvals.lst2)), stringsAsFactors = FALSE)
params.k9.wide <- left_join(params.dat2.wide, dat.pvals.k9)

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
# 
# 
# bins.keep <- c(jbins.eryth, jbins.bcell, jbins.granu, jbins.hspcs)
# 
# bins.keep.lst <- list("Eryths" = jbins.eryth,
#                       "Bcells" = jbins.bcell,
#                       "Granulocytes" = jbins.granu,
#                       "HSPCs" = jbins.hspcs)

jbins.eryth <- bins.keep.lst[["Eryths"]]
jbins.bcell <- bins.keep.lst[["Bcells"]]
jbins.granu <- bins.keep.lst[["Granulocytes"]]
jbins.hspcs <- bins.keep.lst[["HSPCs"]]
bins.keep <- c(jbins.eryth, jbins.bcell, jbins.granu, jbins.hspcs)

bnames <- names(bins.keep.lst); names(bnames) <- bnames


jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
bins.uniq <- unique(bins.keep)
dat.annot <- AnnotateCoordsFromList(coords.vec = bins.uniq, inf.tss = "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.100000.bed", txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)

# assign closest gene and TSS?
params.k9.wide.annot <- left_join(params.k9.wide, subset(dat.annot$regions.annotated, select = c(region_coord, SYMBOL, distanceToTSS)), by = c("bin" = "region_coord"))

# Check giladi  -----------------------------------------------------------

inf.public <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/giladi_pseudobulk_exprs_data.rds"
dat.public <- readRDS(inf.public)

m2c <- MarkerToCelltype()
dat.public$celltype2 <- sapply(as.character(dat.public$celltype), function(x) m2c[[x]])

# check eryth bins

jgenes <- subset(dat.annot$out2.df.closest, region_coord %in% jbins.hspcs)$gene
jgenes <- subset(dat.annot$out2.df.closest, region_coord %in% jbins.bcell)$gene
jgenes <- subset(dat.annot$out2.df.closest, region_coord %in% jbins.granu)$gene
jgenes <- subset(dat.annot$out2.df.closest, region_coord %in% jbins.eryth)$gene
# eryth.genes <- subset(dat.annot$regions.annotated, region_coord %in% jbins.eryth)$SYMBOL

pdf(outpdf, useDingbats = FALSE)
dat.params <- lapply(bnames, function(bname){
  jbins <- bins.keep.lst[[bname]]
  jgenes <- subset(dat.annot$out2.df.closest, region_coord %in% jbins)$gene
  dat.public.sub <- subset(dat.public, gene %in% jgenes) %>%
    ungroup() %>%
    mutate(celltype2 = forcats::fct_reorder(.f = celltype2, .x = zscore, .fun = median, .desc = TRUE),
           celltype = forcats::fct_reorder(.f = celltype, .x = zscore, .fun = median, .desc = TRUE))
  # m <- ggplot(dat.public.sub, aes(x = celltype2, y = zscore)) + 
  #   ggtitle(paste(bname, "Ngenes:", length(jgenes))) +  
  #   geom_boxplot(outlier.shape = NA) + 
  #   geom_point(alpha = 0.25) + 
  #   theme_bw() + 
  #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  # print(m)
  m <- ggplot(dat.public.sub, aes(x = interaction(celltype, celltype2), y = zscore)) + 
    ggtitle(paste(bname, "Ngenes:", length(jgenes))) +  
    geom_boxplot(outlier.shape = NA) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
  m <- ggplot(dat.public.sub, aes(x = interaction(celltype, celltype2), y = exprs)) + 
    ggtitle(paste(bname, "Ngenes:", length(jgenes))) +  
    geom_boxplot(outlier.shape = NA) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
  
  # get annotated params
  params.sub <- subset(params.k9.wide, bin %in% jbins) %>%
    mutate(label = bname, "LowInK9me3" = low.in.k9)
  return(params.sub)
}) %>%
  bind_rows()
dev.off()

# Write DE outputs for k9me3  ---------------------------------------------

fwrite(dat.params, file = outtxt, na = "NA")




