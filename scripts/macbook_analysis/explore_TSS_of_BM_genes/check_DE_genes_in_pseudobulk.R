# Jake Yeung
# Date of Creation: 2020-03-13
# File: ~/projects/scchic/scripts/macbook_analysis/explore_TSS_of_BM_genes/check_DE_genes_in_pseudobulk.R
# Check DE genes in pseudobulk. Lowly expressed genes are not so meaningful??? 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(forcats)

# Load annots -------------------------------------------------------------

inf.annot <- "/Users/yeung/data/scchic/public_data/Giladi_et_al_2018/diff_exprs_Giladi_seurat.celltypes_filt.rds"
dat.annots <- readRDS(inf.annot)
jmark <- "H3K4me1"
jdir <- "lessthan"
logfcmin <- -0.5
pvalmax <- 0.01

MarkerToCelltype <- function(){
  # get correspondance between marker to celltype
  x <- list("Car1" = "Erythroblast",
            "core" = "HSCs",
            "Vpreb1" = "Bcell",
            "Siglech" = "pDendritic",
            "Prg2" = "Eosinophil",
            "Gstm1" = "Neutrophil",
            "Ly86" = "Monocyte",
            "Ccl5" = "NKcell",
            "Prss34" = "Basophil",
            "Cd74" = "cDendritic",
            "Pf4" = "Megakaryocyte",
            "Fcrla" = "Bcell",
            "Fcnb" = "Neutrophil",
            "Hba.a2" = "Erythroblast",
            "Ltf" = "Neutrophil")
  return(x)
}

# Load pseudobulkl  -------------------------------------------------------

inf.pseudo <- "/Users/yeung/data/scchic/public_data/Giladi_et_al_2018/giladi_pseudobulk_datsum_and_ncells.DESeq_and_QuantNorm.RData"
load(inf.pseudo, v=T)

dat.sum.long <- data.table::melt(dat.sum.norm.quantnorm)
colnames(dat.sum.long) <- c("gene", "celltype", "exprs")

dat.sum.long <- dat.sum.long %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE)) %>%
  filter(!is.nan(zscore))

ggplot(dat.sum.long, aes(x = celltype, y = zscore)) + geom_boxplot() + 
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.sum.long, aes(x = celltype, y = exprs)) + geom_boxplot() + 
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
# Find top genes in neutro, plot it ---------------------------------------
annot.lst <- MarkerToCelltype()
jclusts <- as.character(unique(dat.annots$cluster))
outpdf <- paste0("/Users/yeung/data/scchic/pdfs/marker_genes_Giladi_TSS_signal/AllClusts_zscore_across_scrnaseq_pseudobulk.pvalmax_", pvalmax, ".logfcmin_", jdir, "_", logfcmin, ".pdf")

pdf(outpdf, useDingbats = FALSE)
for (jclust in jclusts){
  jsub <- subset(dat.annots, cluster == jclust & p_val_adj <= pvalmax & avg_logFC <= logfcmin)
  jclust.annot <- annot.lst[[jclust]]
  top.genes <- jsub$gene
  jtitle <- paste0(jclust, " (", jclust.annot, "). ",  "N=", length(top.genes), " genes.\nPval<", pvalmax, " logFC ", jdir, " ", logfcmin)
  # rname?
  plot(density(jsub$avg_logFC), main = jtitle)
  # hist(jsub$avg_logFC, main = jtitle, col = 'blue')
    print(jclust)
  jsub$cluster <- sapply(jsub$cluster, function(x) paste0(x, "(", annot.lst[[x]], ")"))
    m <- ggplot(subset(dat.sum.long, gene %in% top.genes) %>%
                  rowwise() %>%
                  mutate(celltype = paste0(celltype, "(", annot.lst[[celltype]], ")")), 
                aes(x = forcats::fct_reorder(celltype, zscore, median, .desc = TRUE), y = zscore)) + 
      geom_boxplot() + geom_jitter(width = 0.1) + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
      xlab("") + 
      ggtitle(jtitle)
    print(m)
    m.exprs <- ggplot(subset(dat.sum.long, gene %in% top.genes) %>%
                  rowwise() %>%
                  mutate(celltype = paste0(celltype, "(", annot.lst[[celltype]], ")")), 
                aes(x = forcats::fct_reorder(celltype, exprs, median, .desc = TRUE), y = exprs)) + 
      geom_boxplot() + geom_jitter(width = 0.1) + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
      xlab("") + 
      ggtitle(jtitle)
    print(m.exprs)
}
dev.off()


