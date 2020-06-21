# Jake Yeung
# Date of Creation: 2020-06-18
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/check_celltype_spec_genes_using_model.R
# Check that the genes im using at celltype sepcific zebrafish to double check the heatmap 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)


jdate <- "2020-06-08"
jdate2 <- "2020-06-09"

gsets.filt <- c("lymphocytes", "erythrocytes", "HSPCs", "granulocytes")
jgenes <- c("meis1b", "tal1", "pax5", "s100z")

ctypes.end.old <- list("Clustergranu", "Clusterlymph", "Clustereryth"); names(ctypes.end.old) <- ctypes.end.old
ctypes.end <- list("ClusterGranulocytes", "ClusterBcells", "ClusterErythroblasts"); names(ctypes.end) <- ctypes.end
ctypes.hash <- hash::hash(ctypes.end.old, ctypes.end)

# gset.specs <- list("granulocytes", "lymphocytes", "erythrocytes"); names(gset.specs) <- ctypes.end
gset.specs <- list("Neutrophil", "Bcell", "Erythro"); names(gset.specs) <- ctypes.end
# gset.others <- list(c("lymphocytes", "erythrocytes"), c("erythrocytes", "granulocytes"), c("lymphocytes", "granulocytes")); names(gset.others) <- ctypes.end
gset.others <- list(c("Bcell", "Erythro"), c("Erythro", "Neutrophil"), c("Bcell", "Neutrophil")); names(gset.others) <- ctypes.end
gsets.differentiated <- c("Neutrophil", "Bcell", "Erythro")
gset.hsc <- "HSPCs"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/setup_matrix_for_poisson_regression.likeBM.redo_count_tables"
assertthat::assert_that(dir.exists(indir))

infits.wrangled <- file.path(indir, paste0("fit_poisson_model_on_TSS_ZF.DownstreamWrangled.",jdate2, ".ClusterRenamed.RData"))

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/integrated_analysis_poisson_and_2D_clouds"
assertthat::assert_that(dir.exists(outdir))

inf.ci <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits/fit_poisson_model_on_TSS.ZF.DownstreamWrangled.2020-06-10.ConfidenceIntervals.smaller.Wrangled.RData"


load(infits.wrangled, v=T)
load(inf.ci, v=T)

# add arrows? check for Bcell-specific genes vs Eryth+Granu specific genes
print(unique(fits.bygenesets.long$geneset))
print(unique(fits.bygenesets.long$cluster))


# Load sets of genes ------------------------------------------------------


ctypes <- c("HSPCs", "Granus", "Bcells", "Eryths")
names(ctypes) <- ctypes


beddir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bedannotations/ZebrafishWKMFromTopics.2000"

bed.lst <- lapply(ctypes, function(ctype){
  inf <- file.path(beddir, paste0("ZebrafishWKM_TSS_FromTopics.", ctype, ".bsize_2.bed"))
  bed <- fread(inf, col.names = c("chromo", "start", "end", "gene"))
})

genes.lst <- lapply(bed.lst, function(jbed){
  jbed$gene
})

m.lst <- lapply(ctypes, function(ctype){
  m <- ggplot(subset(fits.bygenesets.long, gene %in% genes.lst[[ctype]]) %>% filter(abs(logLambda) < 5), aes(x = mark, y = logLambda, fill = cluster)) + geom_boxplot() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
    geom_hline(yintercept = 0) + ggtitle(ctype)
  return(m)
})

print(m.lst)
