t# Jake Yeung
# Date of Creation: 2020-04-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/5-make_pseudobulk_downsample_save_outputs.R
# 

rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(DropletUtils)


library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ChIPseeker)
library(GenomicRanges)
library(forcats)



# Load annots  ------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/rdata_robjs"
# stringent filter louvain? 

for (jmark in jmarks){
  print(jmark)
  
  
  # filter by previously defined celltypes? 
  inf.annot.louv <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/LDA_downstream/LDA_downstream_ZF.2020-04-23.imputevarfilt.lessstringent/ZF_LDA_output.", jmark, ".keepn_150.final.ClusterTables.txt")
  assertthat::assert_that(file.exists(inf.annot.louv))
  
  
  inf.annot.glmpca <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/celltyping/from_louvain/cell_to_cluster_table.", jmark, ".2020-04-13.txt")
  assertthat::assert_that(file.exists(inf.annot.glmpca))
  
  annot.louv <- fread(inf.annot.louv)
  annot.louv$clusterplate <- paste(annot.louv$cluster, annot.louv$plate, "_")
  
  annot.glmpca <- fread(inf.annot.glmpca)
  annot.glmpca.filt <- subset(annot.glmpca, cell %in% annot.louv$cell) %>%
    rowwise() %>%
    mutate(clusterplate = paste(cluster, plate, sep = "_")) %>%
    mutate(cluster = ifelse(cluster %in% c("lymph1", "lymph2"), "lymph", cluster)) %>%   # rename lymph1 and lymph2 into lymph
    ungroup() %>%
    filter(cluster != "Unknown")  # remove the small cluster Unknown
  
  # plot uMAP 
  
  cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
  m <- ggplot(annot.glmpca.filt, aes(x = umap1, y = umap2, color = cluster)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette)
  print(m)
  
  
  
  # Load count TSS data  ----------------------------------------------------
  
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.TSS.CodingOnly.imputevarfilt.lessstringent.mapq_40.winsize_10000/", jmark, ".imputevarfilt.lessstringent.mapq_40.remerged.countTable.TSS.csv")
  mat <- ReadMatTSSFormat(inf)
  
  
  # Pseudobulk --------------------------------------------------------------
  
  cnames.keep.lst <- lapply(split(annot.glmpca.filt, annot.glmpca.filt$cluster), function(jdat) jdat$cell)
  
  pbulk <- SumAcrossClusters(mat, cnames.keep.lst)
  
  pbulk <- do.call(what = cbind, args = pbulk)
  
  pbulk.filt <- CollapseRowsByGene(pbulk, as.long = FALSE, track.kept.gene = TRUE)
  
  
  boxplot(pbulk)
  boxplot(pbulk.filt)
  
  psums <- colSums(pbulk.filt)
  jprop <- min(psums) / psums
  
  # Down sample  ------------------------------------------------------------
  
  pbulk.filt.ds <- DropletUtils::downsampleMatrix(pbulk.filt, prop = jprop)
  
  # get ensembl genes
  genes.orig <- sapply(rownames(pbulk.filt.ds), function(x) strsplit(x, ";")[[1]][[2]])
  genes.ens <- JFuncs::Gene2Ensembl.ZF(genes.orig, return.original = TRUE, species = "drerio")
  
  pbulk.long <- data.frame(genefull = rownames(pbulk.filt.ds), gene = genes.orig, ens = genes.ens, pbulk.filt.ds, stringsAsFactors = FALSE) %>%
    reshape2::melt(., id.vars = c("genefull", "gene", "ens"), variable.name = "pbulk", value.name = "ncuts") %>%
    rowwise() %>%
    mutate(log2cuts = log2(ncuts + 1)) %>%
    group_by(genefull) %>%
    mutate(log2FC = log2cuts - mean(log2cuts),
           log2zscore = log2FC / sd(log2cuts))
  
  # save to output
  outf <- file.path(outdir, paste0("WKM_scChICseq_pseudobulk_downsampled.ens.", jmark, ".RData"))
  save(pbulk.filt.ds, pbulk.filt, annot.louv, annot.glmpca.filt, pbulk.long, file = outf)
  
}