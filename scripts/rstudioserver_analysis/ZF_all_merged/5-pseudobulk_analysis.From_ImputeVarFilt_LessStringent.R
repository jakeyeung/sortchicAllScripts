# Jake Yeung
# Date of Creation: 2020-04-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/5-pseudobulk_analysis.From_ImputeVarFilt_LessStringent.R
# Redo pseudobulk analysis with impute var filt. Maybe try a GLM model across plates? Downsampling? 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(DropletUtils)


# Functions ---------------------------------------------------------------


Gene2Ensembl.ZF <- function(gene.list, return.original=TRUE, species = "drerio"){
  library("biomaRt")
  mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = paste0(species, "_gene_ensembl"), 
                      host = "www.ensembl.org")
  gos <- getBM(gene.list, attributes = c("external_gene_name", 
                                         "ensembl_gene_id"), filters = c("external_gene_name"), 
               mart = mart.obj)
  gl <- gos[match(gene.list, gos[, 1]), 2]
  print(paste0("Could not match ", length(gl[is.na(gl)]), " genes."))
  if (return.original) {
    gl[is.na(gl)] <- gene.list[is.na(gl)]
  }
  return(gl)
}
# Load pbulk scrnaseq -----------------------------------------------------


inf.pbulk <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/rdata_robjs/WKM_pseudobulk_scrnaseq_downsampled.2020-04-26.EosinophilsKeep.RData"
load(inf.pbulk, v=T)

# rename objects to prevent clashing
pbulk.long.gexprs <- pbulk.long; rm(pbulk.long)
pbulk.ctypefilt.long.gexprs <- pbulk.ctypefilt.long; rm(pbulk.ctypefilt.long)


inf.de <- "/home/jyeung/data/from_rstudioserver/zebrafish.2020-04-26/diff_exprs_Chloe_seurat.full.rds"
de.out <- readRDS(inf.de)

plot(density(log10(abs(de.out$avg_logFC) + 1)))


# Load annots  ------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


# stringent filter louvain? 

for (jmark in jmarks){
  pdfout <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis.downsample/pseudobulk_on_DE_genesets.", jmark, ".pdf")
  pdf(pdfout, useDingbats = FALSE)
  
  
  
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
  
  jthres <- 3
  plot(density(log2(pbulk.filt.ds * 1 + 1)))
  abline(v = jthres)
  
  # filter genes where all genes are less than threshold
  
  genes.exprs <- which(apply(log2(pbulk.filt.ds + 1), MARGIN = 1, FUN = function(jrow) max(jrow)) > jthres)
  
  pbulk.filt2.ds <- pbulk.filt.ds[genes.exprs, ]
  
  # Calculate log2 fold changes?  -------------------------------------------
  
  pbulk.long <- data.frame(genefull = rownames(pbulk.filt2.ds), pbulk.filt2.ds, stringsAsFactors = FALSE) %>%
    reshape2::melt(., id.vars = "genefull", variable.name = "pbulk", value.name = "ncuts") %>%
    rowwise() %>%
    mutate(log2cuts = log2(ncuts + 1)) %>%
    group_by(genefull) %>%
    mutate(log2FC = log2cuts - mean(log2cuts),
           log2zscore = log2FC / sd(log2cuts))
    
  m <- ggplot(pbulk.long, aes(x = pbulk, y = log2zscore)) + geom_boxplot() + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
  # Find tissue-specific genes ----------------------------------------------
  
  # pick random genes
  jgenes <- sample(x = rownames(pbulk.filt2.ds), size = 100)
  
  m <- ggplot(pbulk.long %>% filter(genefull %in% jgenes), aes(x = pbulk, y = log2zscore)) + geom_boxplot() + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
  
  
  # Define gene sets from pseudobulk scrnaseq  ------------------------------
  
  pbulk.ctypefilt.long.gexprs <- pbulk.ctypefilt.long.gexprs %>%
    rowwise() %>%
    mutate(ens = strsplit(as.character(gene), "_")[[1]][[1]])
  
  jmean.min <- 4
  fc.min <- 1
  zscore.min <- 0
  
  jctypes <- unique(as.character(pbulk.ctypefilt.long.gexprs$pbulk))
  # jctypes <- unique(as.character(pbulk.long.gexprs$pbulk))
  names(jctypes) <- jctypes
  
  print(jctypes)
  
  jctype <- "lymphocytes"
  jsub <- subset(pbulk.ctypefilt.long.gexprs, pbulk == jctype & log2p1counts >= jmean.min & log2fc >= fc.min & zscore >= zscore.min) %>%
    arrange(desc(abs(log2fc)))
  
  # do for all ctypes
  de.genes.lst <- lapply(jctypes, function(jctype){
    jsub <- subset(pbulk.ctypefilt.long.gexprs, pbulk == jctype & log2p1counts >= jmean.min & log2fc >= fc.min & zscore >= zscore.min) %>%
    # jsub <- subset(pbulk.long.gexprs, pbulk == jctype & log2p1counts >= jmean.min & log2fc >= fc.min & zscore >= zscore.min) %>%
      arrange(desc(abs(log2fc)))
    return(jsub$gene)
  })
  
  # check by plotting zcsores
  
  print(names(de.genes.lst))
  
  mlst <- lapply(names(de.genes.lst), function(jname){
    glst <- de.genes.lst[[jname]]
    m <- ggplot(subset(pbulk.ctypefilt.long.gexprs, gene %in% glst), aes(x = forcats::fct_reorder(pbulk, zscore, median, .desc = TRUE), y = zscore)) + geom_boxplot() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle(paste(jname, "DE genes from scRNAseq pseudobulk"))
    print(m)
  })
  
  print(mlst)
  
  # mlst[[6]]
  # mlst[[5]]
  # mlst[[3]]
  
  # Check gene list in the scChIC-seq  --------------------------------------
  
  # get gene names right?
  glst.scchicseq <- sapply(rownames(pbulk.filt2.ds), function(x) strsplit(x, ";")[[1]][[2]])
  
  glst.ensembl <- Gene2Ensembl.ZF(glst.scchicseq, return.original = TRUE, species = "drerio")
  
  print(unique(length(glst.scchicseq)))
  print(unique(length(glst.ensembl)))
  
  g2e <- hash(unlist(glst.scchicseq), unlist(glst.ensembl))
  
  # convert degenes to ensembl
  de.genes.lst.ens <- lapply(de.genes.lst, function(jgenes){
    sapply(as.character(jgenes), function(x) strsplit(x, "_")[[1]][[1]], USE.NAMES = FALSE)
  })
  
  # assign ensembl name to pbulk
  pbulk.long <- pbulk.long %>%
    rowwise() %>%
    mutate(ens = g2e[[strsplit(genefull, ";")[[1]][[2]]]])
  
  # Check expression of each DE gene list with scChIC-seq  pseudobul --------
  
  print(colnames(pbulk.long))
  
  # jctype <- "granulocytes"
  
  jgenes <- de.genes.lst.ens[[jctype]]
  
  jctypes <- names(de.genes.lst.ens)
  names(jctypes) <- jctypes
  mlst.chic <- lapply(jctypes, function(jctype){
    jgenes <- de.genes.lst.ens[[jctype]]
    m <- ggplot(subset(pbulk.long, ens %in% jgenes), aes(x = forcats::fct_reorder(pbulk, log2FC, median, .desc = TRUE), y = log2zscore)) + 
      geom_boxplot() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle(paste("ctype:", jctype, "Ngenes:", length(jgenes)))
    return(m) 
  })
  
  print(mlst.chic)
  
  
  dev.off()
}



# # Find HSPCs in chic, ask what is its expression in scRNAseq  -------------
# 
# print(unique(pbulk.long$pbulk))
# jmean.min.chic <- 3
# plot(density(pbulk.long$log2cuts))
# plot(density(pbulk.long$log2FC))
# plot(density(pbulk.long$log2zscore, na.rm = TRUE))
# 
# jctype <- "eryth"
# jctype <- "eryth"
# jgenes.bychic <- subset(pbulk.long, pbulk == jctype & log2cuts >= jmean.min.chic & log2FC  >= 0.5 & log2zscore >= 0.25)$ens
# 
# m <- ggplot(subset(pbulk.ctypefilt.long.gexprs, ens %in% jgenes.bychic), aes(x = forcats::fct_reorder(pbulk, zscore, median, .desc = TRUE), y = zscore)) + geom_boxplot() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
# print(m)
# 
# 
# 
# print(mlst)

# Explore -----------------------------------------------------------------
# 
# 
# jgene <- "anxa5b"
# jgene <- "fabp3"
# jgene <- "pax5"
# jgene <- "mmp13"
# jgene <- "gata1a"
# jgene <- "viml"
# jgene <- "lyz"
# jgene <- "meis1b"
# jsub <- subset(pbulk.ctypefilt.long, grepl(jgene, gene)) %>%
#   arrange(desc(abs(log2fc)))
# print(jsub)
# 
# print(subset(de.out, grepl(jgene, gene)))
# 
# ggplot(jsub, aes(x = pbulk, y = log2p1counts)) + geom_col() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# print(jsub)
